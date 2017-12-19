#!/usr/bin/env python
desc="""Calculate mutual information between likely edited sites (REDiscover output). 

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

def header2bams(handle, out=None):
    """Retrieve bam fnames and strandness"""
    header = []
    for l in handle:
        if not l.startswith('##'):
            if out:
                out.write("%s\tmutual information\t%s"%('\t'.join(l.split('\t')[:3]), '\t'.join(l.split('\t')[3:])))
            break
        if out:
            out.write(l)
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

def get_pos2base(a, start, end, baseq, pos2idx, insint=5, delint=6):
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
                if ii not in pos2idx or q<baseq or b not in base2index:
                    continue
                pos2base[ii] = base2index[b]+1
        # deletion
        elif code==2 and prefi in pos2idx:
            pos2base[prefi] = delint
    return pos2base

def _match_bases(arr):
    """Return alphabet fitted between both positions"""
    d, dvals = {}, set()
    c = Counter(tuple(arr[:, i]) for i in xrange(arr.shape[1]))
    for (a, b), _c in c.most_common():
        if b not in d and a not in dvals:
            if not d:
                arr[1][arr[1]==b], a0 = -1, a
            else:
                arr[1][arr[1]==b] = a
            d[b] = a
            dvals.add(a)
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
        # subsample so alt base is as freq as major allele
        altc = c.most_common(2)[1][1]
        cols = [col[0] for b, _c in c.most_common() for col in np.argwhere(arr[ii]==b)[:altc]]
        arr = arr[:, cols]
        # calculate hamming distance between calls & get corresponding bases - hth        
        hammingdist = (_match_bases(arr)[:, None] != arr).sum()/2. 
        # store mutual information
        #print arr.shape[1], arr
        mi = 1 - 1.* hammingdist / arr.shape[1]#; print pos[i], pos[j], hammingdist, mi, c.most_common()
        if mi>mutual_info[i]: mutual_info[i] = mi 
        if mi>mutual_info[j]: mutual_info[j] = mi 
    return mutual_info

def bams2mutual_info(bams, ref, pos, mapq=15, baseq=20, maxdepth=1000000, minfree=100000, maxcov=600):
    """Get mutual info from positions"""
    start, end = pos[0], pos[-1]+1
    sams = [pysam.AlignmentFile(bam).fetch(ref, start, end) for bam in bams]
    calls = np.zeros((len(pos), maxdepth), dtype="int8", order='F')
    posset = set(pos)
    pos2idx = {p: idx for idx, p in enumerate(pos)}
    mutual_info = np.zeros(len(pos))
    iread = p = stops = 0
    while True:
        pa = 0
        for sam in sams:
            ppos = cov = 0
            for a in sam:
                if is_qcfail(a, mapq): # or is_duplicate(a, pa): 
                    continue
                # check cov
                if ppos == a.pos:
                    if cov>maxcov:
                        continue
                    cov += 1
                else:
                    ppos, cov = a.pos, 0
                pa = a
                pos2base = get_pos2base(a, start, end, baseq, pos2idx)
                if not pos2base: 
                    continue
                # make sure not too many reads
                iread += 1
                if iread>=maxdepth:
                    # count how many empty lines
                    notempty = calls[p:].sum(axis=0)>0
                    iread = notempty.sum()
                    # just in case no empty; remove 10% of first reads
                    if iread+minfree > maxdepth:
                        iread -= minfree
                        calls = calls[:, minfree:]                
                    else:
                        calls = calls[:, notempty]
                    # strip info about past reads and add new
                    logger("  %s:%s-%s: resizing array: %s rows left"%(ref, start, end, iread))
                    calls = np.hstack((calls, np.zeros((len(pos), maxdepth-iread), dtype="int8", order='F')))#; print calls.shape
                # store calls
                for _p, _b in pos2base.iteritems():
                    calls[pos2idx[_p]][iread] = _b 
                if a.pos>pos[p]:
                    break
        # calculate mutual info
        while a.pos>pos[p]:
            mutual_info = update_mutual_info(mutual_info, calls, p, pos)#; print bam, ref, pos[p], mutual_info[p]
            p += 1
        if not pa:
            break
    # calculate mutual info
    for p in range(p, len(mutual_info)):
        mutual_info = update_mutual_info(mutual_info, calls, p, pos)
    return mutual_info

def combine_lines(lines):
    """Combine genotypes for the same position. Handle strands separately."""
    if len(lines)<2:
        return lines
    strand2snps = {}
    for l in lines:
        ldata = l[:-1].split("\t")#; print ldata
        snp = ldata[2]
        strand = snp[-1]
        # store new snp
        if strand not in strand2snps:
            strand2snps[strand] = ldata[:2] + [snp[:-1]] + ldata[3:]
        else:
            # add alt allel
            strand2snps[strand][2] += snp[-2]
            # add freq
            for i in range(4, len(ldata), 2):
                strand2snps[strand][i] += ";%s"%ldata[i]
    # combine snp info
    lines = []        
    for strand, l in strand2snps.iteritems():
        l[2] += strand
        lines.append("\t".join(l)+"\n")
    return lines
        
def line_generator(handle):
    """Report lines"""
    pchrom = pp = ''
    lines = []
    for l in handle:
        if l.startswith('#'):
            continue
        # get chrom and pos
        chrom, p = l[:-1].split("\t")[:2]
        p = int(p)-1 # 1-off to 0-off
        # report    
        if chrom != pchrom or p!=pp:
            if lines:
                yield combine_lines(lines)
            pchrom, pp, lines = chrom, p, []
        lines.append(l)
    if lines:
        yield combine_lines(lines)
    
def worker(data):
    # ignore all warnings
    np.seterr(all='ignore')
    bams, ref, pos = data
    logger(" %s:%s-%s"%(ref, pos[0], pos[-1]))    
    return bams2mutual_info(bams, ref, pos)

def pos2mutual(fname, out=sys.stdout, threads=4, verbose=1, log=sys.stderr):
    """Filter out positions with high mutual info"""
    if type(out) is str:
        out = open(out, "w")
    if out.name.endswith(".gz"):
        out = gzip.open(out.name, "w")

    # parse header 
    handle = gzip.open(fname)
    handle2 = gzip.open(fname)
    bams, strands = header2bams(handle, out)
    # process
    if threads<2:
        import itertools
        p = itertools
    else:
        p = Pool(threads)
    logger("Processing chromosomes...")
    for parsers in p.imap(worker, [(bams, ref, pos) for ref, pos in chr2pos(handle)]):
        for lines, mi in izip(line_generator(handle2), parsers):
            for l in lines:
                out.write("%s\t%.3f\t%s"%("\t".join(l.split("\t")[:3]), mi, "\t".join(l.split("\t")[3:])))
    out.close()
            
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
