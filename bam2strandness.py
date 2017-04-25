#!/usr/bin/env python
desc="""Calculate strandness for each bam

It require tab-delimited BED-like file (chr, start, end, strand) exons
that are non-overlapping with antisense genes. This can be generated easily with bedtools:
 bedtools intersect -v -S -a exons.gtf -b genes.gtf | bedtools sort | bedtools merge -i - -s > exons.uniq.bed
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 25/04/2017
"""

import os, sys, pysam, resource, zlib
from datetime import datetime
from multiprocessing import Pool
import numpy as np
from REDiscover import is_antisense, is_qcfail
import random

def load_bed2(fname):
    """Return regions from BED file"""
    chr2intervals = {}
    for l in open(fname):
        if l.startswith('#'):
            continue
        chrom, start, end, strand = l[:-1].split('\t')[:4]
        # start-end    
        start, end = map(int, (start, end))
        # strand
        if strand=="+":
            reverse = False
        else:
            reverse = True
        # store
        if chrom not in chr2intervals:
            chr2intervals[chrom] = []
        chr2intervals[chrom].append((start, end, reverse))
            
    dtype = np.dtype({'names':   ['start',  'end',    'strand'], \
                      'formats': ['uint32', 'uint32', 'bool_']})
    
    for chrom, data in chr2intervals.iteritems():
        data.sort()
        chr2intervals[chrom] = np.array(data, dtype=dtype)        
    return regions

def load_bed(fname):
    """Return regions from BED file"""
    regions = []
    for l in open(fname):
        if l.startswith('#'):
            continue
        ref, start, end, strand = l[:-1].split('\t')[:4]
        start, end = map(int, (start, end))
        regions.append((ref, start, end, strand))
    return regions#[:10000]
'''    
def load_strandness(outfn):
    """Return strandness info from file"""
    ref2strandness = {}
    for l in open(outfn):
        ref, sense, antisense = l[:-1].split('\t')[:3]
        ref2strandness[ref] = map(int, (sense, antisense))
    return ref2strandness
'''    
def bam2strandness(bam, regions, mapq, verbose):
    """Calculate strandness"""
    '''outfn = bam+".strandness.tsv"
    # load
    if os.path.isfile(outfn):
        return load_strandness(outfn)
    
    # calculate from scratch
    ref2strandness = {}'''
    strandness = [0, 0]
    sam = pysam.AlignmentFile(bam)
    for ref, start, end, strand in regions:    
        # stop if ref not in sam file
        if ref not in sam.references:
            return
        # ref: [sense, antisense]
        #if ref not in ref2strandness:
        #    ref2strandness[ref] = [0, 0]
        for a in sam.fetch(ref, start, end):
            if is_qcfail(a, mapq): 
                continue
            # get transcript strand
            if is_antisense(a):
                if strand == "-":
                    i = 0
                else:
                    i = 1
            # sense and +
            elif strand=="+":
                i = 0
            # sense and -
            else:
                i = 1
            #ref2strandness[ref][i] += 1
            strandness[i] += 1
    '''        
    # store
    with open(outfn, "w") as out:
        for chrom, (sense, antisense) in sorted(ref2strandness.iteritems()):
            out.write("%s\t%s\t%s\n"%(chrom, sense, antisense))
    
    return ref2strandness'''
    return strandness

def init_args(*args):
    global regions, mapq, verbose 
    regions, mapq, verbose = args
    
def worker(bam):       
    global regions, mapq, verbose
    return bam2strandness(bam, regions, mapq, verbose)

def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-b", "--bam", nargs="+",  help="input RNA-Seq BAM file(s)")
    parser.add_argument("-e", "--exons", required=1, help="non-overlapping exon file (chr, start, end, strand)")
    parser.add_argument("-s", "--subset", default=1000, type=int, help="select random exons [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=1, type=int, help="number of cores to use [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if o.verbose:
        sys.stderr.write("Loading regions\n")
    regions = load_bed(o.exons)
    if o.verbose:
        sys.stderr.write(" %s exons loaded!\n"%len(regions))
    if o.subset:
        if o.verbose:
            sys.stderr.write("Selecting %s random regions...\n"%o.subset)
        regions = random.sample(regions, o.subset)
        
    if o.verbose:
        sys.stderr.write("Parsing bam file(s)...\n")
    i = 0
    if o.threads<2: # this is useful for debugging
        for bam in o.bam:
            strandness = bam2strandness(bam, regions, o.mapq, o.verbose)#; print strandness
            print bam, sum(strandness), 1.*strandness[0] / sum(strandness)
    else:
        initargs = (regions, o.mapq, o.verbose)
        p = Pool(o.threads, initializer=init_args, initargs=initargs)
        parser = p.imap_unordered(worker, o.bams)
        for i, strandness in enumerate(parser, 0):
            bam = o.bams[i]

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)

