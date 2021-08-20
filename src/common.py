#!/usr/bin/env python
desc="""Identify RNA editing sites from RNAseq and DNAseq alignements (.bam).
Alternatively, reference genome can be used instead of DNAseq,
but at the cost of higher false positive. 

TBD:
- editing from heterozygous sites?
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw/Bratislava/Fribourg, 21/07/2015
"""

import os, sys, pysam, resource
from datetime import datetime
from multiprocessing import Pool
import numpy as np

alphabet = "ACGT"
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()):
    base2index[b] = i

def is_duplicate(a, pa):
    """Return True if read is duplicate"""
    if pa and a.pos==pa.pos and a.flag==pa.flag and a.cigarstring==pa.cigarstring and a.isize==pa.isize and a.seq==pa.seq:
        return True

def is_heavily_clippped(a, maxfrac=0.05):
    """Return True if fraction of read being soft-clipped is larger than maxfrac"""
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    clipped_bases = sum(bases for code, bases in a.cigar if code==4)
    clipped_frac = 1.*clipped_bases/a.rlen
    if clipped_frac>maxfrac:
        return True
        
def is_qcfail(a, mapq=15):
    """Return True if alignment record fails quality checks"""
    if a.mapq<mapq or a.flag&3840 or is_heavily_clippped(a): 
        return True
        
def get_blocks_simple(a, start, end, baseq, i, basesize):
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
        if refi>end:
            bases -= refi-end
        if bases<1:
            break
        block = [0]*basesize*bases 
        for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases])):
            if q<baseq or b not in base2index:
                continue
            block[ii*basesize+base2index[b]+i] += 1
        yield prefi, block

def get_blocks(a, start, end, baseq, i, basesize):
    """Return tuple of aligned position of query and reference"""
    def _match(refi, readi, bases): return refi+bases, readi+bases, True
    def _insertion(refi, readi, bases): return refi, readi+bases, []
    def _softclip(refi, readi, bases): return refi, readi+bases, True
    def _deletion(refi, readi, bases): return refi+bases, readi, []
    def _skip(refi, readi, bases): return refi, readi, []
    code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                     2: _deletion, 3: _deletion, 4: _softclip, 5: _skip}

    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # mark soft-clipped & indels regions
        if code in (1, 2, 4):
            if code in (1, 2) and bases<25:
                bases = 25
            ## mark on both sides
            prefi-=bases
            bases = 2*bases
            if code==2:
                bases += bases
            _refi = prefi + bases
            if prefi<start:
                bases -= start-prefi
                prefi = start
            if _refi>end:
                bases -= _refi-end
            yield prefi, code, bases
            continue
        # skip if not alignment or 
        if not data or refi<start-1:
            continue
        # typical alignment part
        if prefi<start:
            bases -= start-prefi
            preadi += start-prefi
            prefi = start
        if refi>end:
            bases -= refi-end
        if bases<1:
            break
        block = [0]*basesize*bases 
        for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases])):
            if q<baseq or b not in base2index:
                continue
            block[ii*basesize+base2index[b]+i] += 1
            # add fwd strand signals
            if not a.is_reverse: block[ii*basesize+base2index[b]+2*len(alphabet)+i] += 1
        yield prefi, block, 0

def is_antisense(a):
    """Return 1 if read pair is from antisense strand,
    meaning first-in-pair or unpaired aligned to reverse strand, 
    or second-in-pair aligned to forward strand.
    """
    # antisense if reverse
    if a.is_reverse:
        # and first in pair or unpaired
        if a.is_read1 or not a.is_paired:
            return 1
    # or forward and second in pair
    elif a.is_read2:
        return 1
    # otherwise sense
    return 0
        
def bam2calls(bam, ref, start, end, mapq=15, baseq=20, offset=33):
    """Return 2D array of basecalls from BAM file, as follows
    - 1D positions from start to end
    - 2D base counts for ACGT from sense and antisense strand at given position
    """
    sam = pysam.AlignmentFile(bam)
    # ACGT x2 for each strand
    basesize = 2*len(alphabet)
    n =  basesize * (end-start+1)
    calls = np.zeros(n, dtype="int64") 
    # stop if ref not in sam file
    if ref not in sam.references:
        return calls.reshape((end-start+1, basesize))
    pa = None  
    for a in sam.fetch(ref, start, end):
        if is_qcfail(a, mapq) or is_duplicate(a, pa):
            continue
        pa = a
        # get transcript strand
        i = 0 # for +/for i == 0; for -/rev i==len(alphabet)+1
        if is_antisense(a):
            i = len(alphabet)
        for refi, block in get_blocks(a, start, end, baseq, i, basesize):
            s, e = basesize*(refi-start), basesize*(refi-start)+len(block)
            calls[s:e] += block
    return calls.reshape((end-start+1, basesize))

def get_combined_calls(bams, ref, start, end, mapq, baseq, stranded=0):
    """Combine basecalls from several files"""
    parsers = (bam2calls(bam, ref, start, end, mapq, baseq) for bam in bams)
    for call in np.sum(parsers, axis=0): 
        if stranded:
            yield (call[:len(alphabet)], call[len(alphabet):])
        else:
            yield (call[:len(alphabet)] + call[len(alphabet):])
        
def fasta2calls(fastafn, ref, start, end, cov=100):
    """Return list of basecalls from FastA file."""
    fasta = pysam.FastaFile(fastafn)
    if ref not in fasta.references:
        raise StopIteration
    for b in fasta.fetch(ref, start, end):
        call = [0]*len(alphabet)
        if b in base2index:
            call[base2index[b]] += cov
        yield call

def get_calls(dna, rna, fasta, stranded, ref, start, end, mapq, baseq, minDepth):
    """Return basecalls from multiple BAM (& FastA) file(s)"""
    # define strands
    if not stranded:
        strands = "."  # unstranded
    elif stranded=="firststrand":
        strands = "-+" # dUTP, NSR, NNSR
    else:
        strands = "+-" # Illumina or Standard Solid
    # get parsers    
    rnaparser = get_combined_calls(rna, ref, start, end, mapq, baseq, stranded)
    if dna:
        dnaparser = get_combined_calls(dna, ref, start, end, mapq, baseq, stranded=0)
    else:
        dnaparser = fasta2calls(fasta, ref, start, end)
    # process
    for pos, (dnacall, rnacalls) in enumerate(zip(dnaparser, rnaparser), start+1):
        if sum(dnacall)<minDepth:
            continue
        '''for strand, rnacall in zip(strands, rnacalls):
            if sum(rnacall)<minDepth:
                continue
            yield ref, pos, strand, dnacall, rnacall'''
        strand_info = sorted(zip(rnacalls, strands), key=lambda x: sum(x[0]), reverse=1)
        yield ref, pos, dnacall, strand_info

def get_major_alleles(bases, freqs, maxfrac=0.33):
    """Return two lists: major alleles and their frequeinces"""
    if len(bases)<2:
        return bases, freqs
    tokeep = set(i for i, f in enumerate(freqs) if f >= maxfrac * max(freqs))
    return [b for i, b in enumerate(bases) if i in tokeep], [f for i, f in enumerate(freqs) if i in tokeep]
                    
def get_allele_freqs(counts, minFreq=0.03, minCount=3):
    """Return two lists: alleles that passed filtering and their frequencies"""
    bases, freqs = [], []
    for c, b in zip(counts, alphabet):
        freq = 1.*c/sum(counts)
        # skip if alt base calling if less than 3 reads or low freq
        if c >= minCount and freq >= minFreq:
            bases.append(b)
            freqs.append(freq)
    return bases, freqs

def region2editing(dna, rna, fasta, stranded, minDepth, minDNAfreq, minRNAfreq, mapq, baseq, 
                   report_alternatives, verbose, ref="", start="", end=""):
    """Return RNA editing positions"""
    parser = get_calls(dna, rna, fasta, stranded, ref, start, end, mapq, baseq, minDepth)
    #for contig, pos, strand, refbases, sbases in parser:
    for contig, pos, refbases, strand_info in parser:
        baseRef, refFreq = get_allele_freqs(refbases, minDNAfreq)
        if not baseRef: continue
        # keep only major allele(s) - not needed currently
        #baseRef, refFreq = get_major_alleles(baseRef, refFreq)

        pbases = set()
        for sbases, strand in strand_info:
            if sum(sbases)<minDepth: continue
            
            bases, freqs = get_allele_freqs(sbases, minRNAfreq)
            if not bases or bases==baseRef: continue
            
            # remove ref base from alternative bases
            refBase = baseRef[0]
            if refBase in bases:
                idx = bases.index(refBase)
                bases.pop(idx)
                freqs.pop(idx)

            if not report_alternatives:
                # keep only major allele(s) from sample, but after removing ref allele
                bases, freqs = get_major_alleles(bases, freqs)

                if len(bases) > 1:
                    info = "[WARNING] More than 1 alternative allele: %s:%s %s %s %s\n"
                    sys.stderr.write(info%(contig, pos, ",".join(baseRef), ",".join(bases), str(freqs)))
                    continue

            # report 
            for altbase, altfreq in zip(bases, freqs):
                # skip signal from antisense if altbase already added and 5x lower coverage
                ## this is not perfect, as it'll ignore also higher freq editing from sense
                ## if antisense has much higher coverage but low freq
                if altbase in pbases: # and 2.*sum(sbases)<pcov:
                    continue
                yield contig, pos, strand, refBase, altbase, sum(refbases), refFreq[0], sum(sbases), altfreq
            # store variables
            pbases = set(bases)

def init_args(*args):
    global dna, rna, fasta, stranded, minDepth, minDNAfreq, minRNAfreq, mapq, bcq, report_alternatives, verbose
    dna, rna, fasta, stranded, minDepth, minDNAfreq, minRNAfreq, mapq, bcq, report_alternatives, verbose = args
    
def worker(args):
    global dna, rna, fasta, stranded, minDepth, minDNAfreq, minRNAfreq, mapq, bcq, report_alternatives, verbose
    ref, start, end = args
    totdata = []
    for data in region2editing(dna, rna, fasta, stranded, minDepth, minDNAfreq, minRNAfreq,
                               mapq, bcq, report_alternatives, verbose, ref, start, end):
        totdata.append(data)
    return totdata

def get_consecutive(data, stepsize=1):
    """Return consecutive windows allowing given max. step size"""
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)
            
def get_covered_regions(bams, mincov=3, mapq=10, stranded=0, maxdist=16000, step=1000000,
                        window_size=100, min_strandness=0.99):
    """Return chromosome regions covered by at least mincov"""
    sam = pysam.Samfile(bams[0])
    references, lengths = sam.references, sam.lengths
    for ref, length in zip(references, lengths):
        #if ref.startswith('chr'): continue 
        #if ref!="chr21": continue
        coverage = np.zeros(length, dtype='uint16')
        strands = np.zeros(((length/window_size)+1, 2), dtype='float32')
        for i, bam in enumerate(bams):
            sam = pysam.Samfile(bam)
            for a in sam.fetch(reference=ref):
                if is_qcfail(a, mapq):
                    continue
                # store strand info
                strands[a.pos/window_size, is_antisense(a)] += 1
                # add alg blocks
                for s, e in a.blocks:
                    coverage[s:e] += 1
            # get strand enrichment
            sums = strands.sum(axis=1)
            frac = (strands[sums>1].max(axis=1)/sums[sums>1]).mean()#; print ref, bam, frac
            if stranded and frac<min_strandness:
                sys.stderr.write("[WARNING] Low strandness (%.3f) in %s from chromosome %s\n"%(frac, bam, ref))
        # get regions with coverage
        covered = np.where(coverage>=mincov)[0]
        for positions in get_consecutive(covered, maxdist):
            if len(positions)<1:
                continue
            s, e = positions[0]+1, positions[-1]+1
            # further split regions for max 1M windows
            while s < e-step:
                yield ref, s, s+step
                s += step
            yield ref, s, e

def logger(info, add_timestamp=1, add_memory=1, out=sys.stderr):
    """Report nicely formatted stream to stderr"""
    memory = timestamp = ""
    if add_timestamp:
        timestamp = "[%s] "%datetime.ctime(datetime.now())
    if add_memory:
        selfmem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.
        childrenmem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.
        memory = " [memory: %7.1f Mb]"%(childrenmem + selfmem, ) #; %7.1f Mb self
    out.write("%s%s%s\n"%(timestamp, info, memory))

def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-o", "--output", required=True,  help="output file")
    parser.add_argument("-r", "--rna", nargs="+",  help="input RNA-Seq BAM file(s)")
    refpar = parser.add_mutually_exclusive_group(required=True)
    refpar.add_argument("-d", "--dna", nargs="*", default = [],  help="input DNA-Seq BAM file(s)")
    refpar.add_argument("-f", "--fasta", default = None,  help="reference FASTA file")
    parser.add_argument("-s", "--stranded", "-fr-secondstrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. Illumina or Standard Solid")
    parser.add_argument("-fr-firststrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. dUTP, NSR, NNSR")
    parser.add_argument("-a", "--report_alternatives", default=False, action="store_true", 
                        help="report sites with more than 1 alternative alleles")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    #parser.add_argument("--minAltReads", default=3,  type=int,
    #                    help="minimum no. of reads with alternative base to call RNA editing [%(default)s]")
    parser.add_argument("--minRNAfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("--minDNAfreq",  default=0.99, type=float,
                        help="min frequency for genomic base [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("--bcq", default=20, type=int, help="basecall quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="number of cores to use [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
     
    # mark stranded protocol
    if o.fr_firststrand:
        o.stranded = "firststrand"
    
    # check if all input files exists
    for fn in o.dna+o.rna:
        if not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)

    # check if outfile exists and not empty
    if o.output=="-":
        output = sys.stdout
    elif os.path.exists(o.output) and open(o.output).readline():
        sys.stderr.write("The output file %s exists!\n"%o.output)
        sys.exit(1)
    else:
        output = open(o.output, "w")

    runinfo = " ".join(sys.argv)
    header = "## %s\n# chr\tpos\tstrand\tref\talt\tref cov\tref freq\talt cov\talt freq\n"%runinfo
    output.write(header)
    output.flush()
    
    logger("Indexing bam file(s)...")
    for fn in o.dna + o.rna:
        if not os.path.isfile(fn+".bai"):
            cmd = "samtools index %s"%fn
            if o.verbose:
                sys.stderr.write(" %s\n"%cmd)
            os.system(cmd)    

    logger("Genotyping...")
    info = "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%.3f\n"
    regions = get_covered_regions(o.rna, o.minDepth, o.mapq, o.stranded)
    if o.threads<2: # this is useful for debugging
        for ref, start, end in regions:
            parser = region2editing(o.dna, o.rna, o.fasta, o.stranded, o.minDepth, o.minDNAfreq, o.minRNAfreq, \
                                    o.mapq, o.bcq, o.report_alternatives, o.verbose, ref, start, end)
            for data in parser:
                output.write(info%data)
    else:
        initargs = (o.dna, o.rna, o.fasta, o.stranded, o.minDepth, o.minDNAfreq, o.minRNAfreq, \
                    o.mapq, o.bcq, o.report_alternatives, o.verbose)
        p = Pool(o.threads, initializer=init_args, initargs=initargs)
        parser = p.imap(worker, regions)
        for data in parser:
            output.write("".join(info%d for d in data))

    output.write("#Finished!\n")
    logger("Done!")
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
