#!/usr/bin/env python
desc="""Calculate strandness for each bam

It require tab-delimited BED-like file (chr, start, end, strand) exons
that are non-overlapping with antisense genes. This can be generated easily with bedtools:
 bedtools intersect -v -S -a <(awk '$3=="exon"' gtf) -b <(awk '$3=="gene"' gtf) | bedtools merge -i - -s > exons.uniq.bed
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 25/04/2017
"""

import os, sys, pysam, random, resource, zlib
from datetime import datetime
from multiprocessing import Pool
import numpy as np
from REDiscover import is_antisense, is_qcfail

def load_exons(fname, verbose=0):
    """Return regions from gtf/gff file"""
    uniqexons = fname+".exons.uniq.bed"
    if not os.path.isfile(uniqexons):
        cmd1 = """awk '$3=="exon"' %s > %s.exon.gtf"""%(fname, fname)
        cmd2 = """awk '$3=="gene"' %s > %s.gene.gtf"""%(fname, fname)
        cmd3 = "bedtools intersect -v -S -a %s.exon.gtf -b %s.gene.gtf | bedtools sort | bedtools merge -i - -s > %s"%(fname, fname, uniqexons)
        if verbose:
            sys.stderr.write(" Generating merged exons that are non-overlapping with antisense genes as: %s\n  %s\n"%(uniqexons, "\n  ".join((cmd1, cmd2, cmd3))))
        os.system(cmd1); os.system(cmd2); os.system(cmd3)
        
    regions = []
    for l in open(uniqexons):
        if l.startswith('#'):
            continue
        ref, start, end, strand = l[:-1].split('\t')[:4]
        start, end = map(int, (start, end))
        regions.append((ref, start, end, strand))
    return regions

def bam2strandness(bam, regions, mapq, verbose):
    """Calculate strandness"""
    strandness = [0, 0]
    sam = pysam.AlignmentFile(bam)
    for ref, start, end, strand in regions:    
        # stop if ref not in sam file
        if ref not in sam.references:
            return

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
            strandness[i] += 1

    return strandness

def init_args(*args):
    global regions, mapq, verbose 
    regions, mapq, verbose = args
    
def worker(bam):       
    global regions, mapq, verbose
    return bam2strandness(bam, regions, mapq, verbose)

def process_bams(bams, gtf, mapq, subset, threads, verbose=0):
    """Process all bam files and yield number of reads and strandness"""
    if verbose:
        sys.stderr.write("Loading regions\n")
    regions = load_exons(gtf, verbose)
    if verbose:
        sys.stderr.write(" %s exons loaded!\n"%len(regions))
    if subset:
        n = int(round(subset*len(regions)))+1
        if verbose:
            sys.stderr.write(" Selecting %s random regions...\n"%n)
        regions = random.sample(regions, n)

    if not regions:
        sys.exit("No regions!\n")
        
    if verbose:
        sys.stderr.write("Parsing bam file(s)...\n")
    if threads<2: # this is useful for debugging
        for bam in bams:
            strandness = bam2strandness(bam, regions, mapq, verbose)
            yield bam, sum(strandness), 1.*strandness[0] / sum(strandness)
    else:
        initargs = (regions, mapq, verbose)
        p = Pool(threads, initializer=init_args, initargs=initargs)
        parser = p.imap_unordered(worker, bams)
        for i, strandness in enumerate(parser, 0):
            bam = bams[i]
            yield bam, sum(strandness), 1.*strandness[0] / sum(strandness)
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-b", "--bam", nargs="+",  help="input RNA-Seq BAM file(s)")
    parser.add_argument("-g", "--gtf", required=1, help="gtf/gff file with genes and exons")
    parser.add_argument("-s", "--subset", default=0.01, type=float, help="select fraction of exons [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="number of cores to use [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    print "#fname\tno. of reads sampled\tfraction of reads from sense strand"
    for data in process_bams(o.bam, o.gtf, o.mapq, o.subset, o.threads, o.verbose):
        print "\t".join(map(str, data))

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)

