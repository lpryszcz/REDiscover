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
    for i, l in enumerate(gzip.open(fname)): 
        ldata = l[:-1].split('\t')
        if l.startswith('#') or len(ldata)<3:
            if i==3:
                #names = [os.path.basename(fn).split('.')[eid] for fn in l[:-1].split()[1:]]; print len(names), names
                names = [os.path.basename(fn).split()[0].split('.')[eid] for fn in l[:-1].split('\t')[3::2]]; print len(names), names
                snps = [[[] for n in names] for _i in range(12)]
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
        for ii, (cov, freq) in enumerate(sampledata):
            if cov:
                snps[snp2id[snp]][ii].append(freq)
    return snps, names

def plot_hist(bins, names, outfn, snps, title): #, colors, log=0):
    """Plot hist"""
    x, y = 5, 4
    fig, subplots = plt.subplots(x, y, sharey='all', sharex='all', figsize=(15, 15))
    fig.suptitle(title, size=20)
    for ii, freq in enumerate(snps):
        ax = subplots[ii/y][ii%y] 
        n, bins, patches = ax.hist(freq, bins, normed=0)#, stacked=True, color=colors, label=["+", "-"])
        ax.set_title(names[ii])
        ax.grid(True)
    # make y axis log
    #if log: plt.yscale('log')
    plt.savefig(outfn, dpi=300, transparent=False) #orientation='landscape', 
    
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
    parser.add_argument("-e", "--eid", default=-2, type=int,
                        help="element of bam name (after . splitting) [%(default)s]")
    
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
            outfn = "%s.%s.png"%(fn, snp[:-1].replace('>','_'))
            print i, snp, outfn
            plot_hist(o.bins, names, outfn, snps[i], snp[:-1])

if __name__=="__main__":
    main()
