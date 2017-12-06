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

def load_gtf(fname, window=5000):
    """Return chr2intervals from gtf"""
    gtf = {}
    for l in open(fname):
        chrom, provider, ftype, s, e = l[:-1].split('\t')[:5]
        s, e = int(s)-window, int(e)+window
        if not chrom.startswith('chr'):
            chrom='chr%s'%chrom
        if chrom not in gtf:
            gtf[chrom] = []
        gtf[chrom].append((s, e))
    return gtf

def load_snps(fname, snp2id, eid=-2, gtf={}, minDepth=5, minFreq=0.05, minAltReads=3, minSamples=1):
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
                names = [os.path.basename(fn).split()[0].split('.')[eid] for fn in l[:-1].split('\t')[3::2]]; print len(names), names
                snps = [np.zeros(len(names), dtype=int) for _i in range(12)]
            continue
        chrom, pos, snp = ldata[:3]
        # unstranded
        if snp[-1] == ".":
            snp = snp.replace(".", "+")
        if snp not in snp2id:
            continue
        if not chrom.startswith('chr'):
            chrom = "chr%s"%chrom
        #print chrom, pos, snp
        # skip if SNP not overlapping with GTF
        if chrom not in gtf or not filter(lambda x: x[0]<=int(pos)<=x[1], gtf[chrom]):
            continue
        sampledata = np.array(map(float, ldata[3:])).reshape(len(ldata[3:])/2, 2)
        # enough depth, frequency and more than 3 reads in alt allele
        passed = (sampledata[:, 0]>=minDepth) & (sampledata[:, 1]>=minFreq) \
                 & (sampledata[:, 0]*sampledata[:, 1]>=minAltReads)
        # store
        if sum(passed)<minSamples:
            continue
        snps[snp2id[snp]] += passed
    return snps, names
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="file(s) to process")
    parser.add_argument("-g", "--gtf", default="", help="intervals in GTF to intersect [%(default)s]")
    parser.add_argument("-e", "--eid", default=0, type=int,
                        help="element of bam name (after . splitting) [%(default)s]")
    parser.add_argument("--ext", default="png", choices=['png', 'svg', 'pdf', 'jpg'], 
                        help="figure extension [%(default)s]")
    parser.add_argument("-w", "--window", default=5000, type=int, help="window size [%(default)s]")
    
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

    gtf = load_gtf(o.gtf, o.window)
            
    for fn in o.fnames:
        snps, names = load_snps(fn, snp2id, o.eid, gtf)
        for i, snp in enumerate(id2snp):
            print "%s\t"%snp + "\t".join(str(x) for x in snps[i])

if __name__=="__main__":
    main()
