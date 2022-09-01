#!/usr/bin/env python3
desc="""Calculate mutual information between likely edited sites (REDiscover output). 

TBD:
- since we already captured all positions with variants, now we can focus only on those
this should greatly speed up this step
- still having problem if there is low freq mis-alignment and true motif
ie. rRNA.LSU.23S.Escherichia_coli.prokaryotic_cytosol:272 A>CG. C is mis-alignment, but G is true editing.
- modifications with RT signature up-/down-stream
ie. K0 P-1 rRNA.LSU.23S.Escherichia_coli.prokaryotic_cytosol:745     G>ACTd. 0.607   743     0.017;0.020;0.074;0.061
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw/Mizerow, 6/12/2017
"""

import os, sys, pysam, resource, zlib, gzip
from datetime import datetime
from collections import Counter
from multiprocessing import Pool
import numpy as np
from REDiscover import logger, alphabet, base2index, code2function, get_consecutive, is_qcfail, is_duplicate

def header2bams(handle, out=None):
    """Retrieve bam fnames and strandness"""
    header = []
    for l in handle:
        if not l.startswith('#'):
            if out:
                out.write("%s\tmutual information\t%s"%('\t'.join(l.split('\t')[:3]), '\t'.join(l.split('\t')[3:])))
            break
        if out:
            out.write(l)
        header.append(l[:-1].split()[3:])
    return header[-3], header[-2]

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

def store_pos2base(a, start, end, baseq, pos2idx, calls, iread, delint=6):
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
        return 1
    return 0

def _match_bases(arr, verbose=0):
    """Return alphabet fitted between both positions"""
    b2a, a2b = {}, {}
    arr[1] += -len(alphabet)
    c = Counter(tuple(arr[:, i]) for i in range(arr.shape[1]))
    for (a, b), _c in c.most_common():
        if b not in b2a and a not in a2b:
            b2a[b] = a
            a2b[a] = b
            arr[1][arr[1]==b] = a
    if verbose: return c.most_common()
    return arr
    
def update_mutual_info(mutual_info, calls, i, pos, minCommonReads=5, posi=2203):
    """Calculate mutual information between given position and all other positions
    
    Note, This can be phased all at once by change in read order and computing some similarity.
    """
    # process only positions having at least minCommonReads
    for j, in np.argwhere(np.sum(calls[i+1:,calls[i]>0]>0, axis=1)>=minCommonReads)+i+1:
        arr = calls[[i,j]][:, np.all(calls[[i,j]], axis=0)]
        # subsample for low freq sites
        altc = ii = 0
        c, c1 = Counter(arr[0]), Counter(arr[1])
        if len(c)<2 or len(c1)<2 or c.most_common(2)[1][1]<minCommonReads or \
           c1.most_common(2)[1][1]<minCommonReads:
            continue
        # match bases at both positions
        _match_bases(arr)
        # subsample so alt base is as freq as major allele
        ## get all mutual bases at two positions
        c = Counter(zip(arr[0], arr[1]))
        # count how many alternative bases
        altc = arr.shape[1] - c.most_common(1)[0][1]
        # get subset of major allele
        arr = np.delete(arr, np.argwhere(np.all(arr==c.most_common(1)[0][0][0], 0))[altc:].T[0], 1)
        # store mutual information - 0.5-1.0 instead of 0.5-1.0
        #mi = 1 - 1.* np.sum(arr[0] != arr[1]) / arr.shape[1]
        mi = 2*(0.5 - 1.* np.sum(arr[0] != arr[1]) / arr.shape[1])
        #mutual_info[i].append(mi); mutual_info[j].append(mi)
        if mi>mutual_info[i]: mutual_info[i] = mi 
        if mi>mutual_info[j]: mutual_info[j] = mi
    return mutual_info

def clean_calls(bams, calls, pairs, ref, pos, p, iread, minfree, maxdepth, dtype, order):
    """Remove empty rows from calls to store new reads."""
    if iread>=maxdepth:
        logger("  cleaning cache (%s:%s)..."%(ref, pos[p]))
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
        pairs = [{} for bam in bams]
    return calls, pairs, iread
    
def bams2mutual_info(bams, ref, pos, mapq=15, baseq=20, maxcov=600, dtype="int8", order="C"):
    """Get mutual info from positions"""
    logger(" %s:%s-%s"%(ref, pos[0], pos[-1]))
    maxdepth = maxcov*(len(bams)+1)
    minfree = int(maxdepth/10)
    calls = np.zeros((len(pos), maxdepth), dtype=dtype, order=order)
    mutual_info = np.zeros(len(pos))
    if len(pos)<2: return mutual_info
    start, end = pos[0], pos[-1]+1
    #print(bams)
    sams = [pysam.AlignmentFile(bam).fetch(ref, start, end) for bam in bams]
    pairs = [{} for bam in bams]
    posset = set(pos)
    pos2idx = {p: idx for idx, p in enumerate(pos)}
    iread = p = 0
    while True:
        pa = 0
        for sami, sam in enumerate(sams):
            cov = 0
            for a in sam:
                pa = a
                if a.pos>pos[p]:
                    break
                # check cov
                if cov>maxcov or is_qcfail(a, mapq): # or is_duplicate(a, pa): 
                    continue
                cov += 1
                # store calls and count stored reads
                ## pair of previously stored read - store in row that pair occupies
                if a.qname in pairs[sami]:
                    store_pos2base(a, start, end, baseq, pos2idx, calls, pairs[sami][a.qname])
                ## new row for a pair
                elif store_pos2base(a, start, end, baseq, pos2idx, calls, iread):
                    pairs[sami][a.qname] = iread
                    iread += 1
                # make sure not too many reads
                calls, pairs, iread = clean_calls(bams, calls, pairs, ref, pos, p, iread, minfree, maxdepth, dtype, order)    
                
        # calculate mutual info
        while a.pos>pos[p]:
            mutual_info = update_mutual_info(mutual_info, calls, p, pos)
            p += 1
        # make sure not too many reads
        calls, pairs, iread = clean_calls(bams, calls, pairs, ref, pos, p, iread, minfree, maxdepth, dtype, order)    
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
        elif snp[-2] not in strand2snps[strand][2]: #else:
            # add alt allele
            strand2snps[strand][2] += snp[-2]
            # add freq
            for i in range(4, len(ldata), 2):
                strand2snps[strand][i] += ";%s"%ldata[i]
    # combine snp info
    lines = []
    for strand, l in strand2snps.items():
        l[2] += strand
        lines.append("\t".join(l)+"\n")
    return lines
        
def line_generator(handle):
    """Report lines"""
    pchrom = pp = ''
    lines = []
    for l in handle:
        if l.startswith(('#', 'chromosome\tposition')):
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
    
def worker(args):
    # ignore all warnings
    #np.seterr(all='ignore')
    return bams2mutual_info(*args) 

def pos2mutual(fname, out=sys.stdout, threads=4, mapq=15, bcq=20, maxcov=600, verbose=1, log=sys.stderr):
    """Filter out positions with high mutual info"""
    if type(out) is str: out = open(out, "wt")
    if out.name.endswith(".gz"): out = gzip.open(out.name, "wt")

    # parse header 
    handle = gzip.open(fname, "rt")
    bams, strands = header2bams(handle, out)#; print(bams, strands)
    # process
    if threads<2:
        import itertools
        p = itertools
    else:
        p = Pool(threads)
        
    logger("Processing chromosomes...")
    mutual_info = []
    for outdata in p.imap(worker, [(bams, ref, pos, mapq, bcq, maxcov) for ref, pos in chr2pos(handle)]):
        mutual_info += list(outdata)
    #logger("Saving...")
    for mi, lines in zip(mutual_info, line_generator(gzip.open(fname, "rt"))):
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
    parser.add_argument("-q", "--mapq", default=3, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("-Q", "--bcq", default=20, type=int, help="basecall quality [%(default)s]")
    parser.add_argument("--maxcov", default=600, type=int, help="max coverage per sample [%(default)s]")

    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    pos2mutual(o.input, o.out, o.threads, o.mapq, o.bcq, o.maxcov, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
