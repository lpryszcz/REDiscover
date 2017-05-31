#!/usr/bin/env python
desc="""Plot histogram from REDiscover.diff output.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 231/05/2017
"""

import os, sys, gzip
import numpy as np
from collections import Counter
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt

bases = "ACGT"
strands = "+-"
base2rc= {"A": "T", "T": "A", "C": "G", "G": "C", ">": ">", "+": "-", "-": "+"}

def load_snps(fname, snp2id, eid=-2, dbSNP={}, minDepth=5, minFreq=0.05, minAltReads=3, minSamples=1):
    """Load SNP into list of lists:
    - snp type
    - freqs of each sample for give snp type
    """
    # process snps
    snps, names = [], []
    #name2id = {}
    for i, l in enumerate(gzip.open(fname)): 
        ldata = l[:-1].split('\t')
        if l.startswith('#') or len(ldata)<3:
            if i==3:
                #names = [os.path.basename(fn).split('.')[eid] for fn in l[:-1].split()[1:]]; print len(names), names
                names = [os.path.basename(fn).split()[0].split('.')[eid] for fn in l[:-1].split('\t')[3::2]]; print len(names), names
                snps = [[[] for n in names] for _i in range(12)]
                '''for n in names:
                    if n not in name2id:
                        name2id[n] = len(name2id)
                snps = [[[] for n in name2id.keys()] for _i in range(12)]'''
            continue
        chrom, pos, snp = ldata[:3]
        # unstranded
        if snp[-1] == ".":
            snp = snp.replace(".", "+")
        if snp not in snp2id:
           continue
        if not chrom.startswith('chr'):
            chrom = "chr%s"%chrom
        # skip if present in dbSNP
        if dbSNP and chrom in dbSNP and int(pos) in dbSNP[chrom]:
            continue
        sampledata = np.array(map(float, ldata[3:])).reshape(len(ldata[3:])/2, 2)
        # enough depth, frequency and more than 3 reads in alt allele    
        passed = sum((sampledata[:, 0]>=minDepth) & (sampledata[:, 1]>=minFreq) \
                     & (sampledata[:, 0]*sampledata[:, 1]>=minAltReads))
        # store
        if passed<minSamples:
            continue
        #print sampledata.shape
        for ii, (cov, freq) in enumerate(sampledata):
            if cov: #name2id[names[ii]]
                snps[snp2id[snp]][ii].append(freq)
    return snps, names

def _sorter(x):
    xlist = x.split('_')
    name2id = {'solid':0, 'total': 1, 'bound': 2, 'unbound': 3}
    data = [name2id[xlist[0]]]
    if xlist[1].lower()=="egg": data.append(0)
    elif xlist[1].startswith('128'): data.append(6)
    elif xlist[1].startswith('16'): data.append(4)
    elif xlist[1].startswith('1'): data.append(2)
    elif xlist[1].startswith('3'): data.append(8)
    elif xlist[1].startswith('5'): data.append(10)
    else: data.append(xlist[1])
    return data 
    
def plot_hist(bins, orgnames, outfn, snps, title, selected=[], startswith=''): #, colors, log=0):
    """Plot hist"""
    name2id = {}
    # subset of samples
    if selected or startswith:
        names = filter(lambda x: x in selected or x.startswith(startswith), orgnames)
    else:
        names = orgnames
    # get number of rows/columns
    ncol = nrow = int(np.sqrt(len(names)))
    if ncol*nrow < len(names):
        ncol += 1
    if ncol*nrow < len(names):
        nrow += 1
    #nrow, ncol = 1, 3
    # sort names
    #names = sorted(names, key=lambda x: _sorter(x)); print names; nrow, ncol = 4, 6; name2id = {n: i for i, n in enumerate(names)}
    
    # plot
    fig, subplots = plt.subplots(nrow, ncol, sharey='all', sharex='all', figsize=(ncol*3, nrow*3+2))
    fig.suptitle(title, size=20)
    i = 0
    for ii, freq in enumerate(snps):
        if selected and orgnames[ii] not in names or startswith and not orgnames[ii].startswith(startswith):
            continue
        if name2id:
            i = name2id[orgnames[ii]]
        if nrow>1:
            ax = subplots[i/ncol][i%ncol]
        else:
            ax = subplots[i%ncol] 
        n, bins, patches = ax.hist(freq, bins, normed=0)#, stacked=True, color=colors, label=["+", "-"])
        ax.set_title(orgnames[ii])
        ax.grid(True)
        i += 1
    # make y axis log
    #if log: plt.yscale('log')
    plt.savefig(outfn, dpi=100, transparent=False)#, orientation='landscape')
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="file(s) to process")
    parser.add_argument("-b", "--bins", default=50, type=int,
                        help="number of bins in histogram [%(default)s]")
    parser.add_argument("-e", "--eid", default=0, type=int,
                        help="element of bam name (after . splitting) [%(default)s]")
    parser.add_argument("-s", "--selected", nargs="+", default=[], help="select samples [all]")
    parser.add_argument("-s2", "--startswith", default='', help="select samples that start with string [all]")
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
    for a in bases:
        for b in bases:
            if a==b: continue
            snp = "%s>%s%s"%(a, b, s)
            snprc = "".join(base2rc[b] for b in snp)
            snp2id[snp] = len(id2snp)
            snp2id[snprc] = len(id2snp)
            id2snp.append(snp)

    for fn in o.fnames:
        print fn
        snps, names = load_snps(fn, snp2id, o.eid)

        for i, snp in enumerate(id2snp):
            #if not snp.startswith('A>G'): continue
            outfn = "%s.%s.%s"%(fn, snp[:-1].replace('>','_'), o.ext)
            print i, snp, outfn
            plot_hist(o.bins, names, outfn, snps[i], snp[:-1], o.selected, o.startswith)

if __name__=="__main__":
    main()
