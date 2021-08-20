#!/usr/bin/env python
desc="""Calculate mutual information for SNP sites
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Fribourg, 10/08/2017
"""

import os, sys, pysam, resource, zlib
from datetime import datetime
from multiprocessing import Pool
import numpy as np

from REDiscover import logger, base2index, alphabet, is_antisense, is_qcfail, is_duplicate
import bam2strandness


def get_blocks(a, start, end, baseq, i):
    """Return tuple of aligned position of query and reference"""
    def _match(refi, readi, bases): return refi+bases, readi+bases, True
    def _insertion(refi, readi, bases): return refi, readi+bases, []
    def _deletion(refi, readi, bases): return refi+bases, readi, []
    def _skip(refi, readi, bases): return refi, readi, []
    code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                     2: _deletion, 3: _deletion, 4: _insertion, 5: _skip}
    
    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # skip if not alignment or 
        if not data or refi<start-1:
            continue
        # typical alignment part
        if prefi<start:
            bases -= start-prefi
            preadi += start-prefi
            prefi = start
        # break if block outside of region
        if refi>end:
            break
        yield prefi, preadi, bases


def get_mi(bams, ref, start, end, snps):
    """Return mutual information for snps"""
    sam = pysam.AlignmentFile(bam)
    
    # process alignments
    for a in sam.fetch(ref, start, end):
        if max(snps)<a.pos or min(snps)>a.aend \
          or is_qcfail(a, mapq) or is_duplicate(a, pa):
            continue        
        # process alignment blocks
        for refi, readi, bases in get_blocks(a, start, end, baseq, i):
            # intersect with snps
            for snpi in set(range(refi, bases+1)).intersection(snps)):

def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-o", "--out", required=1, help="output file")
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), help="output stream [stdin]")
    

       
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
