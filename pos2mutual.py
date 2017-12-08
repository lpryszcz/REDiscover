#!/usr/bin/env python
desc="""Calculate mutual information between likely edited sites and ignore those
being SNPs or alignment issues

TBA:
- utilise PE info
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 6/12/2017
"""

import os, sys, pysam, resource, zlib, gzip
from datetime import datetime
from itertools import izip
from collections import Counter
from multiprocessing import Pool
import numpy as np
from REDiscover import logger, alphabet, base2index, code2function, get_consecutive, is_qcfail, is_duplicate

def header2bams(handle):
    """Retrieve bam fnames and strandness"""
    header = []
    for l in handle:
        if not l.startswith('##'):
            break
        header.append(l[3:-1])
    return header[1].split(), header[2].split()  

def chr2pos(handle, out=sys.stdout, window=10000): 
    """Load editing candidates"""
    pchrom, pos = '', []
    for l in handle:
        if l.startswith('#'):
            out.write(l)
            continue
        # get chrom and pos
        chrom, p = l[:-1].split()[:2]
        p = int(p)-1 # 1-off to 0-off
        # report    
        if chrom != pchrom:
            if pos:
                for pos in get_consecutive(pos, window):
                    yield pchrom, list(pos)
            pchrom, pos = chrom, []
        # store if different than previous
        if not pos or p!=pos[-1]:
            pos.append(p)
    if pos:
        for pos in get_consecutive(pos, window):
            yield pchrom, list(pos)

def get_pos2base(a, start, end, baseq, insint=5, delint=6):
    """Return dictionary of position and corresponding base. Include DEL but ignore INS."""
    pos2base = {}
    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # typical alignment part
        if data and refi>start:
            if prefi<start:
                bases -= start-prefi
                preadi += start-prefi
                prefi = start
            if refi>end:
                bases -= refi-end
            if bases<1:
                break
            for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases]), prefi):
                if q<baseq or b not in base2index:
                    continue
                pos2base[ii] = base2index[b]+1
        # deletion
        elif code==2:
            pos2base[prefi] = delint
    return pos2base

def _match_bases(arr):
    """Return alphabet fitted between both positions"""
    d = {}
    c = Counter(tuple(arr[:, i]) for i in range(arr.shape[1]))
    for (a, b), _c in c.most_common():
        if b not in d and a not in d.values():
            if not d:
                arr[1][arr[1]==b] = -1
                a0 = a
            else:
                arr[1][arr[1]==b] = a
            d[b] = a
    arr[1][arr[1]==-1] = a0
    return arr
        
def update_mutual_info(mutual_info, calls, i, pos, minCommonReads=5):
    """Calculate mutual information between given position and all other positions"""
    # process only positions having at least minCommonReads
    for j, in np.argwhere(np.sum(calls[i+1:,calls[i]>0]>0, axis=1)>=minCommonReads)+i+1:
        arr = calls[[i,j]][:, np.all(calls[[i,j]], axis=0)]
        # subsample for low freq sites
        altc = ii = 0
        c = Counter(arr[0])#; print c
        if len(c)<2 or c.most_common(2)[1][1]<minCommonReads:
            c = Counter(arr[1])
            if len(c)<2 or c.most_common(2)[1][1]<minCommonReads:
                continue
            ii = 1
        altc = c.most_common(2)[1][1]
        cols = [col[0] for b, _c in c.most_common() for col in np.argwhere(arr[ii]==b)[:altc]]#; print c, cols
        arr = arr[:, cols]
        # calculate hamming distance between calls & get corresponding bases - hth        
        hammingdist = (_match_bases(arr)[:, None] != arr).sum()/2. 
        # store mutual information
        #print arr.shape[1], arr
        mi = 1 - 1.* hammingdist / arr.shape[1]#; print pos[i], pos[j], hammingdist, mi, c.most_common()
        if mi>mutual_info[i]: mutual_info[i] = mi 
        if mi>mutual_info[j]: mutual_info[j] = mi 
    return mutual_info
    
def bam2mutual_info(bam, stranded, ref, pos, mapq=15, baseq=20, maxdepth=100000):
    """Get mutual info from positions"""
    start, end = pos[0], pos[-1]+1
    sam = pysam.AlignmentFile(bam)
    calls = np.zeros((len(pos), maxdepth), dtype="int8", order='F')
    posset = set(pos)
    pos2idx = {p: idx for idx, p in enumerate(pos)}
    mutual_info = np.zeros(len(pos)) 
    iread = pa = 0
    for a in sam.fetch(ref, start, end):
        if is_qcfail(a, mapq) or is_duplicate(a, pa): 
            continue
        pos2base = get_pos2base(a, start, end, baseq)
        if not posset.intersection(pos2base):
            continue
        # make sure not too many reads
        iread += 1
        if iread>=maxdepth:
            # count how many empty lines
            empty = calls.sum(axis=0)==0#; print empty.sum()
            iread -= empty.sum()
            # strip info about past reads and add new
            logger(" Resizing array: %s rows left"%iread)
            calls = calls[:, empty]
            calls = np.hstack((calls, np.zeros((len(pos), empty.sum()), dtype="int8", order='F')))
        # report
        while a.pos>pos[0]:
            mutual_info = update_mutual_info(mutual_info, calls, 0, pos)
            yield mutual_info[0] 
            # update calls & pos2idx
            pos, calls, mutual_info = pos[1:], calls[1:], mutual_info[1:]
            pos2idx = {p: idx for idx, p in enumerate(pos)}
        # store calls
        for p in posset.intersection(pos2base):
            calls[pos2idx[p]][iread] = pos2base[p]
    # report the rest
    for i in range(calls.shape[0]):
        mutual_info = update_mutual_info(mutual_info, calls, i, pos)
        yield mutual_info[i]
    
def line_generator(handle):
    """Report lines"""
    pchrom = pp = ''
    lines = []
    for l in handle:
        if l.startswith('#'):
            continue
        # get chrom and pos
        chrom, p = l[:-1].split()[:2]
        p = int(p)-1 # 1-off to 0-off
        # report    
        if chrom != pchrom or p!=pp:
            if lines:
                # need to combine
                yield lines[0]
            pchrom, pp, lines = chrom, p, []
        lines.append(l)
    if lines:
        yield lines[0]
    
def worker(data):
    # ignore all warnings
    np.seterr(all='ignore')
    bams, strands, ref, pos = data
    logger(" %s:%s-%s"%(ref, pos[0], pos[-1]))    
    #parsers = [bam2mutual_info(bam, strand, ref, pos) for bam, strand in zip(bams, strands)]
    data = [] #d for d in zip(parsers)]
    for bam, strand in zip(bams, strands):
        data.append([d for d in bam2mutual_info(bam, strand, ref, pos)])
    return data

def pos2mutual(fname, out=sys.stdout, threads=4, verbose=1, log=sys.stderr):
    """Filter out positions with high mutual info"""\
    # parse header 
    handle = gzip.open(fname)
    handle2 = gzip.open(fname)
    bams, strands = header2bams(handle)
    # process
    logger("Processing chromosomes...")
    """
    for ref, pos in chr2pos(handle, out):
        logger(" %s:%s-%s"%(ref, pos[0], pos[-1]))
        #if ref!="tRNA.Gly.GCC.Bacillus_subtilis.prokaryotic_cytosol": continue
        parsers = [bam2mutual_info(bam, strand, ref, pos) for bam, strand in zip(bams, strands)]        
        for i, data in enumerate(izip(line_generator(handle2), *parsers)):
            l = data[0]
            mi = np.array(data[1:])
            #print "%s:%s\t%.3f\t%s"%(ref, pos[i]+1, mi[mi>0].mean() if mi.sum() else 0, str(mi))#) #, mi[mi>0]
            out.write("%s\t%.3f\t%s"%("\t".join(l.split("\t")[:3]), mi[mi>0].mean() if mi.sum() else 0, "\t".join(l.split("\t")[3:])))
    """        
    # process
    p = Pool(threads)
    for parsers in p.imap(worker, [(bams, strands, ref, pos) for ref, pos in chr2pos(handle, out)]):
        for i, data in enumerate(izip(line_generator(handle2), *parsers)):
            l = data[0]
            mi = np.array(data[1:])
            out.write("%s\t%.3f\t%s"%("\t".join(l.split("\t")[:3]), mi[mi>0].mean() if mi.sum() else 0, "\t".join(l.split("\t")[3:])))
    #"""    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")
    parser.add_argument('--version', action='version', version='1.2a')
    parser.add_argument("-i", "--input", default='', help="REDiscover variation file")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"), help="output file")
    parser.add_argument("-t", "--threads", default=4, type=int, help="number of cores to use [%(default)s]")

    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    pos2mutual(o.input, o.out, o.threads, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
