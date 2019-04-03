#!/usr/bin/env python
desc="""Generate violin plot from REDiscover.diff output.

TBA:
- half-violin plot for sense/antisense? https://stackoverflow.com/questions/29776114/half-violin-plot
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Fribourg, 14/08/2017
"""

import os, sys, gzip
import numpy as np
from itertools import izip
from datetime import datetime
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt
#from plot_hist import bases, base2rc

bases = "ACGTid"
base2index = {b: i for i, b in enumerate(bases)}
index2base = {i: b for i, b in enumerate(bases)}
strands = "+-"
base2rc= {"A": "T", "T": "A", "C": "G", "G": "C",
          ">": ">", "+": "-", "-": "+", "i": "i", "d": "d"}

def load_snps(fname, snp2id, id2snp, eid=-2, dbSNP={}, minDepth=10, minFreq=0.00,
              minAltReads=1, minSamples=1, log=sys.stderr):
    """Load SNP into list of lists:
    - snp type
    - freqs of each sample for give snp type
    """
    # process snps
    snps, names = [], []
    #name2id = {}
    for i, l in enumerate(gzip.open(fname)): 
        ldata = l[:-1].split('\t')
        if l.startswith('#'):
            if l.startswith('#chr'):
                names = [os.path.basename(fn).split()[0].split('.')[eid] for fn in l[:-1].split('\t')[3:]]
                sys.stderr.write("Processing variations from %s files: %s...\n"%(len(names), names))
                snps = [[[] for n in names] for _i in range(len(id2snp))]
            continue
        # print info
        if log and not i%10000: log.write(" %s \r"%i)
        if l[-1]!='\n':
            log.write("[WARNING] Interrupting due to corrupted line: %s\n"%l)
            break    
        # unload pos info
        chrom, pos, snp = ldata[:3]
        # unstranded
        if snp[-1] == ".":
            snp = snp.replace(".", "+")
        # there may be multiple changes on single position at given strand ie. T>CG-
        refbase = snp[0]
        refbasei = base2index[refbase]
        snptmp = "%s>%s%s"%(refbase, '%s', snp[-1])
        #if snp not in snp2id:
        #   continue
        if not chrom.startswith('chr'):
            chrom = "chr%s"%chrom
        # skip if present in dbSNP
        if dbSNP and chrom in dbSNP and int(pos) in dbSNP[chrom]:
            continue
        # unload base counts for pos # A,C,G,T,i,d 0,0,0,0,0,0     0,1,0,0,0,0
        sampledata = np.array([map(int, sc.split(',')) for sc in ldata[3:]])#; print(refbase, refbasei, sampledata)
        for ii in range(len(names)):
            # check depth
            cov = sampledata[ii].sum()
            if cov<minDepth: continue
            # get alt reads
            for bi, b in enumerate(sampledata[ii]):
                # skip reference or sites with too little alt reads
                if bi==refbasei or b<minAltReads: continue
                # check freq
                freq = 1.*b/cov
                if freq<minFreq: continue
                # store
                snps[snp2id[snptmp%index2base[bi]]][ii].append(freq)

    return snps, names

def violin_plot(outfn, snps, fnames, id2snp, xmax=1.0):
    """Generate violin plot"""
    fig, axes = plt.subplots(nrows=1, ncols=len(id2snp), figsize=(20, .5*len(fnames)+2))

    # process all files
    pos = range(1, len(fnames)+1)
    for ax, _snps, name in izip(axes, snps, id2snp):
        # get ax and set title
        name = "%s\n%s"%(name[:-1], max([len(x) for x in _snps]))
        ax.set_title(name, fontsize=10)
        ax.set_xlim(0, xmax)
        # violin plot after http://matplotlib.org/examples/statistics/violinplot_demo.html
        ax.violinplot(_snps, pos, points=1000, widths=.5, vert=False,
                      showmeans=True, showextrema=True, showmedians=True, bw_method='silverman')

    # show labels only of first plot
    for i, ax in enumerate(axes.flatten()):
        # strip leading 0 on x-axis tick labels
        #ax.set_xticklabels([0, xmax])#[item.get_text().lstrip('0') for item in ax.get_xticklabels()])
        # print tick labels only on the left-most y-axis
        if i:
            ax.set_yticklabels([])  
        else:
            ax.set_yticks(pos)
            ax.set_yticklabels(fnames)
            
    #fig.suptitle("RNA editing violin plot", fontsize=16)
    fig.subplots_adjust(hspace=.2, left=.15, right=.99, bottom=.05)
    plt.savefig(outfn, dpi=300, transparent=False)
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="file(s) to process")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType('w'), help="output stream [stdout]")
    parser.add_argument("-e", "--eid", default=0, type=int,
                        help="element of bam name (after . splitting) [%(default)s]")
    parser.add_argument("-s", "--selected", nargs="+", default=[], help="select samples [all]")
    parser.add_argument("-s2", "--startswith", default='', help="select samples that start with string [all]")
    parser.add_argument("-x", "--xmax", default=1.0, type=float, help="limit x axis [%(default)s]")
    parser.add_argument("--ext", default="png", choices=['png', 'svg', 'pdf', 'jpg'], 
                        help="figure extension [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    # get snp2id
    snp2id, id2snp = {}, []
    s = "+"
    # skip indels in ref
    for a in bases[:4]: 
        for b in bases[:4]:
            if a==b: continue
            snp = "%s>%s%s"%(a, b, s)
            snprc = "".join(base2rc[b] for b in snp)
            snp2id[snp] = len(id2snp)
            snp2id[snprc] = len(id2snp)
            id2snp.append(snp)
    # store indels at the end
    for a in bases[:4]: 
        for b in bases[4:]:
            if a==b: continue
            snp = "%s>%s%s"%(a, b, s)
            snprc = "".join(base2rc[b] for b in snp)
            snp2id[snp] = len(id2snp)
            snp2id[snprc] = len(id2snp)
            id2snp.append(snp)
            
    # process input files
    out, log = o.out, sys.stderr
    for fn in o.fnames:
        snps, names = load_snps(fn, snp2id, id2snp, o.eid)
        # report stats
        out.write("## %s\n#snp\t%s\n"%(fn, "\t".join(names)))
        for i, snp in enumerate(id2snp):
            lens = map(len, snps[i]) #map(len, [filter(lambda x: 0.02<x<0.98, _snps) for _snps in snps[i]])
            out.write("%s\t%s\n"%(snp[:-1], "\t".join(map(str, lens))))
            
        outfn = "%s.violin_plot.%s"%(fn, o.ext)
        violin_plot(outfn, snps, names, id2snp, o.xmax)
        sys.stderr.write("[INFO] Violin plot saved as: %s\n"%outfn)

if __name__=="__main__":
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)