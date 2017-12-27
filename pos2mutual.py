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

def store_pos2base(a, start, end, baseq, pos2idx, calls, iread, insint=5, delint=6):
    """Store calls from individual reads. Include DEL but ignore INS."""
    stored = False
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
                if ii in pos2idx and q>=baseq and b in base2index:
                    calls[pos2idx[ii], iread] = base2index[b]+1
                    stored = True
        # deletion
        elif code==2 and prefi in pos2idx:
            calls[pos2idx[prefi], iread] = delint
    if stored:
        iread += 1
    return iread, calls

def _match_bases0(arr):
    """Return alphabet fitted between both positions"""
    d, dvals = {}, set()
    c = Counter(tuple(arr[:, i]) for i in xrange(arr.shape[1]))
    for (a, b), _c in c.most_common():
        if not d:
            arr[1][arr[1]==b], a0 = -1, a
        # skip stored match
        elif b in d or a in dvals:
            arr[1][arr[1]==b] = a
        d[b] = a
        dvals.add(a)
    arr[1][arr[1]==-1] = a0
    return arr
    
def _match_bases(arr, verbose=0):
    """Return alphabet fitted between both positions"""
    b2a, a2b, bset = {}, {}, set()
    c = Counter(tuple(arr[:, i]) for i in xrange(arr.shape[1]))
    if verbose: print c.most_common()
    for (a, b), _c in c.most_common():
        if b in b2a or a in a2b:
            continue
        b2a[b] = a
        a2b[a] = b
        arr[1][arr[1]==b] = -a
        bset.add(b)
    # change missing b
    for b in filter(lambda x: x not in b2a, bset):
        arr[1][arr[1]==b] = 4+b
    # * -1
    arr[1] *= -1
    return arr
        
def update_mutual_info(mutual_info, calls, i, pos, minCommonReads=5):
    """Calculate mutual information between given position and all other positions"""
    #return mutual_info
    # process only positions having at least minCommonReads
    for j, in np.argwhere(np.sum(calls[i+1:,calls[i]>0]>0, axis=1)>=minCommonReads)+i+1:
        arr = calls[[i,j]][:, np.all(calls[[i,j]], axis=0)]
        # subsample for low freq sites
        altc = ii = 0
        c, c1 = Counter(arr[0]), Counter(arr[1])#; print c
        if len(c)<2 or len(c1)>1 and c.most_common(2)[1][1]<c1.most_common(2)[1][1]:
            c = c1
            ii = 1
        if len(c)<2 or c.most_common(2)[1][1]<minCommonReads:
            continue
        # subsample so alt base is as freq as major allele
        altc = c.most_common(2)[1][1]
        arr = np.delete(arr, np.argwhere(arr[ii]==c.most_common(1)[0][0])[altc:].T[0], 1)
        arr0 = np.copy(arr)
        # calculate hamming distance between calls & get corresponding bases - hth        
        hammingdist = (_match_bases(arr)[:, None] != arr).sum()/2. 
        # store mutual information
        #print arr.shape[1], arr
        mi = 1 - 1.* hammingdist / arr.shape[1]#; print pos[i], pos[j], hammingdist, mi, c.most_common()
        if mi>mutual_info[i]: mutual_info[i] = mi 
        if mi>mutual_info[j]: mutual_info[j] = mi
        if mi>.9:
            if 263<pos[i]<266 or 263<pos[j]<266:
                print pos[i], pos[j], mi, c.most_common(4)
                print arr.shape, arr, Counter(arr[0]).most_common(4), Counter(arr[1]).most_common(4)
                print Counter(arr0[0]).most_common(4), Counter(arr0[1]).most_common(4), _match_bases(arr0, 1)
    return mutual_info

def bams2mutual_info(bams, ref, pos, mapq=15, baseq=20, maxcov=600, dtype="int8", order="C"):
    """Get mutual info from positions"""
    maxdepth = maxcov*(len(bams)+1)
    minfree = maxdepth/10
    calls = np.zeros((len(pos), maxdepth), dtype=dtype, order=order)
    mutual_info = np.zeros(len(pos))
    if len(pos)<2: return mutual_info
    start, end = pos[0], pos[-1]+1
    sams = [pysam.AlignmentFile(bam).fetch(ref, start, end) for bam in bams]
    posset = set(pos)
    pos2idx = {p: idx for idx, p in enumerate(pos)}
    iread = p = stops = 0
    while True:
        pa = 0
        for sam in sams:
            ppos = cov = 0
            for a in sam:
                if a.pos>pos[p]:
                    break
                # check cov
                if cov>maxcov or is_qcfail(a, mapq) or is_duplicate(a, pa): 
                    continue
                cov += 1
                pa = a
                # store calls and count stored reads
                iread, calls = store_pos2base(a, start, end, baseq, pos2idx, calls, iread)
                # make sure not too many reads
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
                    calls = np.hstack((calls, np.zeros((len(pos), maxdepth-iread), dtype=dtype, order=order)))
                
        # calculate mutual info
        while a.pos>pos[p]:
            mutual_info = update_mutual_info(mutual_info, calls, p, pos)
            #if 263<pos[p]<266: print ref, pos[p], mutual_info[p]
            p += 1
        if not pa:
            break
    # calculate mutual info
    for p in range(p, len(mutual_info)):
        mutual_info = update_mutual_info(mutual_info, calls, p, pos)
        #if 263<pos[p]<266: print ref, pos[p], mutual_info[p]
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
            # add alt allele
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
