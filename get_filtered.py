#!/usr/bin/env python
desc="""Return enrichment of various types of editing
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw/Bratislava/Fribourg, 21/07/2015
"""

import os, sys, gzip
import numpy as np
from remove_dbSNP import load_dbSNP

def get_filtered(fn, snps, minDepth=20, minFreq=0.01, minSamples=3):
    """Return dictionary of snps and their occurencies"""
    i = k = 0
    out = open(fn+".parsed.txt", "w")
    if fn.endswith('.gz'):
        handle = gzip.open(fn)
    else:
        handle = open(fn)
    for l in handle:
        ldata = l.replace('\t\t','\t')[:-1].split('\t')
        if l.startswith('#') or not l.endswith('\n') or len(ldata)<3:
            out.write(l)
            continue
        i += 1
        chrom, pos, snp = ldata[:3]
        sampledata = np.array(map(float, ldata[3:])).reshape(len(ldata[3:])/2, 2)
        passed = sum((sampledata[:, 0]>=minDepth) & (sampledata[:, 1]>=minFreq))# & (sampledata[:, 1]<1.0))
        if passed<minSamples:
            continue
        out.write(l)
        k += 1
    return i, k
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="files to preocess")
    parser.add_argument("-s", "--snps", default="", help="dbSNP file")
    parser.add_argument("-d", "--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("-f", "--minAltfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("-n", "--minsamples", default=3, type=int, help="number of samples [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    snps = load_dbSNP(o.snps)

    print "#fname\tediting\tpassed"
    for fn in o.fnames:
        i, k = get_filtered(fn, snps, o.minDepth, o.minAltfreq, o.minsamples)
        print "%s\t%s\t%s"%(fn, i, k)
    
if __name__=="__main__":
    main()
    
