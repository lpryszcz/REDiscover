#!/usr/bin/env python3
desc="""Identify differential RNA editing sites from multiple RNAseq experiments (.bam).

Dependencies: pysam numpy matplotlib pandas

TBD:
- cap at some max coverage ie. 800X
- process ref bam in the same filtering run
- optimise pos2mutual
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw/Bratislava/Fribourg/Barcelona, 21/07/2015
"""

import gzip, os, sys, pysam, resource, zlib
from datetime import datetime
from multiprocessing import Pool
import numpy as np
from bam2strandness import bam2strandness, is_antisense, is_qcfail, is_duplicate

VERSION = '1.0a'
alphabet = "ACGT" #idse" 
base2name = {"A": "A", "C": "C", "G": "G", "T": "T",
             "i": "insertions", "d": "deletions",
             "s": "read_starts", "e": "read_ends"}
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()): base2index[b] = i

# CIGAR operations
"""Op BAM Description +1Q +1R
M 0 alignment match (can be a sequence match or mismatch) yes yes
I 1 insertion to the reference yes no
D 2 deletion from the reference no yes
N 3 skipped region from the reference no yes
S 4 soft clipping (clipped sequences present in SEQ) yes no
H 5 hard clipping (clipped sequences NOT present in SEQ) no no
P 6 padding (silent deletion from padded reference) no no
= 7 sequence match yes yes
X 8 sequence mismatch yes yes
    """
def _match(refi, readi, bases): return refi+bases, readi+bases, True
def _insertion(refi, readi, bases): return refi, readi+bases, False
def _deletion(refi, readi, bases): return refi+bases, readi, False
def _skip(refi, readi, bases): return refi, readi, False
code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                 2: _deletion, 3: _deletion, 4: _insertion, 5: _skip}
'''
def store_blocks(a, start, end, baseq, i, calls):
    """Store base calls from aligned blocks. INDEL aware."""
    # store start & end counts - only if not soft-clipped and if withing calls (could be another block)
    if a.cigar[0][0]!=4 and 0<=a.pos-start<calls.shape[0]:
        calls[a.pos-start, i, 6] += 1
    if a.cigar[-1][0]!=4 and 0<=a.aend-start<calls.shape[0]:
        calls[a.aend-start-1, i, 7] += 1
    # process read blocks and store bases counts and indels as well
    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # skip if current before start
        if refi<=start:
            continue
        # typical alignment part
        if data:
            if prefi<start:
                bases -= start-prefi
                preadi += start-prefi
                prefi = start
            if refi>end:
                bases -= refi-end
            if bases<1:
                break
            for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases])):
                if q>=baseq and b in base2index:
                    calls[prefi-start+ii, i, base2index[b]] += 1
        elif start<prefi<end:
            # insertion
            if code==1:
                calls[prefi-start, i, 4] += 1
            # deletion
            elif code==2:
                calls[prefi-start, i, 5] += 1
    return calls

def name2barcode(name):
    """Return barcode from read name"""
    return name.split(".")[0]

def bam2calls_old(sam, stranded, ref, start, end, mapq=15, baseq=20, barcoded=False):
    """Return 3D array of basecalls from BAM file, as follows:
    - 1D positions from start to end of the ref
    - 2D sense and antisense strand
    - 3D base counts for ACGTidse
    """
    # position, strand, ACGTid
    calls = np.zeros((end-start+1, 2, len(alphabet)), dtype="int64")#; print(calls.shape)
    # stop if ref not in sam file
    if ref not in sam.references:
        if ref.startswith('chr') and ref[3:] in sam.references:
            ref = ref[3:]
        elif 'chr%s'%ref in sam.references:
            ref = 'chr%s'%refgg
        else:
            return calls
    pa = None
    barcodes = {}        
    for a in sam.fetch(ref, start, end):
        if is_qcfail(a, mapq):
            continue
        # skip PCR duplicates from barcoded experiments
        if barcoded:
            if name2barcode(a.qname) in barcodes and barcodes[name2barcode(a.qname)]!=a.qname:
                continue
            barcodes[name2barcode(a.qname)] = a.qname
        #elif is_duplicate(a, pa):
        #    continue
        pa = a
        # get transcript strand
        i = 0 # for +/for i == 0; for -/rev i==1
        if stranded=="firststrand": 
            if not is_antisense(a):
                i = 1
        # unstranded or secondstrand
        elif is_antisense(a):
            i = 1
        # store alignment blocks
        calls = store_blocks(a, start, end, baseq, i, calls)
    return calls
'''

def bam2calls(sam, stranded, ref, start, end, mapq, baseq, barcoded=False):
    """Return 3D array of basecalls from BAM file, as follows:
    - 1D positions from start to end of the ref
    - 2D sense and antisense strand
    - 3D base counts for alphabet

    Note, this doesn't report insertions, deletions, reads starts and ends, 
    and doesn't support UMI, but it's more than 10x faster than bam2calls_old. 
    """
    def is_positive(a): return not is_qcfail(a) and not is_antisense(a)
    def is_negative(a): return not is_qcfail(a) and is_antisense(a)

    if stranded=="firststrand": is_positive, is_negative = is_negative, is_positive
    else: is_positive, is_negative = is_positive, is_negative
    
    calls = np.zeros((end-start, 2, len(alphabet)), dtype="uint16")
    calls[:, 0, :4] = np.array(sam.count_coverage(ref, start, end, quality_threshold=baseq,
                                                  read_callback=is_positive)).T
    calls[:, 1, :4] = np.array(sam.count_coverage(ref, start, end, quality_threshold=baseq,
                                                  read_callback=is_negative)).T
    return calls

def get_combined_calls(sams, ref, start, end, mapq, baseq, stranded, bam2calls):
    """Combine basecalls from several files""" 
    parsers = (bam2calls(sam, stranded, ref, start, end, mapq, baseq) for sam in sams) 
    for call in np.sum(parsers, axis=0):
        if stranded: yield (call[0], call[1])
        else: yield call[0] + call[1]

def fasta2calls(faidx, ref, start, end, cov=100):
    """Return list of basecalls from FastA file."""
    if ref not in faidx: raise StopIteration
    calls = np.zeros((end-start, len(alphabet)))
    seq = faidx.fetch(ref, start, end).upper()
    for p, b in enumerate(seq):
        if b in base2index: calls[p, base2index[b]] += cov
        yield calls[p]
    
def diff_editing(position, fasta, minDepth, minDNAfreq, minAltfreq, minAltReads,
                 stranded, mapq, baseq, barcoded, allpositions, verbose):
    """Return RNA editing positions

    Note, dna_sams [pysam.AlignmentFile(bam) for bam in dna] 
    and sams [pysam.AlignmentFile(bam) for bam in rna] 
    are provided as global vars to save time (~0.5s per region). 
    """
    # define strands - strands from individual files are handled by bam2calls
    if not stranded[0]: strands = "." # unstranded
    else: strands = "+-"
    info = []
    ref, start, end = position
    start -= 1 # here pos should be 0-based, but 1-based in output
    if dna_sams:
        refparser = get_combined_calls(dna_sams, ref, start, end, mapq, baseq, 0, bam2calls)
    else:
        #faidx = pysam.FastaFile(fasta)
        refparser = fasta2calls(faidx, ref, start, end)
    #sams = [pysam.AlignmentFile(bam) for bam in bams]
    parsers = [bam2calls(sam, _stranded, ref, start, end, mapq, baseq, barcoded)
               for sam, _stranded in zip(sams, stranded)]
    for pos, calls in enumerate(zip(refparser, *parsers), start+1):
        refcall, bamcalls = calls[0], np.array(calls[1:])#; print pos, calls
        # low dna coverage ie. N
        refcov = sum(refcall)
        if not allpositions and refcov < minDepth:
            continue
        refi = np.argmax(refcall)
        reffreq = 1.*refcall[refi]/refcov
        if not allpositions and reffreq < minDNAfreq:
            continue
        refbase = alphabet[refi]
        # collapse both strands for unstranded data, without changing dimensions
        if not stranded[0]:
            bamcalls[:, 0] += bamcalls[:, 1]
        #print(pos, refbase, bamcalls[:, 0], bamcalls.shape)
        # process strands
        for si, strand in enumerate(strands):
            # skip if low coverage in all samples
            coverage = bamcalls[:, si].sum(axis=1)
            if not allpositions and coverage.max() < minDepth: continue
            # process alt bases
            altbases = []
            for bi, base in enumerate(alphabet[:4]): #:6
                # skip reference and if less than 3 reads per allele
                if bi==refi or bamcalls[:, si, bi].max()<minAltReads: continue
                # store if freq of at least one larger than minfreq
                freqs = 1.* bamcalls[:, si, bi] / coverage
                if np.nanmax(freqs) >= minAltfreq: altbases.append(base)
            if altbases or allpositions:
                text = "%s\t%s\t%s>%s%s\t%s\n"%(ref, pos, refbase, "".join(altbases), strand, 
                                                "\t".join("\t".join(map(str, bamcalls[i, si]))
                                                          for i in range(bamcalls.shape[0])))
                if text:
                     info.append(text)
    logger(" %s:%s-%s"%position)
    return "".join(info)
    
def worker(args):
    # ignore all warnings
    np.seterr(all='ignore')  
    return diff_editing(*args)
    
def worker_coverage(args):
    """return coverage"""
    bam, ref, mapq = args
    sam = pysam.Samfile(bam)
    ref2len = {r: l for r, l in zip(sam.references, sam.lengths)}    
    return get_coverage(sam, ref, mapq, ref2len[ref])

def get_coverage(sam, ref, mapq, rlen):
    coverage = np.zeros(rlen, dtype='uint16')
    for a in sam.fetch(reference=ref):
        if is_qcfail(a, mapq): continue
        if a.blocks:
            for s, e in a.blocks: coverage[s:e] += 1
        else: coverage[a.pos:a.aend] += 1
    return coverage

def get_consecutive(data, stepsize=1):
    """Return consecutive windows allowing given max. step size"""
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)

def get_covered_regions_per_bam(bams, mincov=3, mapq=15, threads=4, chrs=[], verbose=0,
                                maxdist=100, step=10000):
    """Return chromosome regions covered by at least mincov."""
    sams = [pysam.AlignmentFile(bam, threads=threads) for bam in bams]
    sam = sams[0]
    references, lengths = sam.references, sam.lengths
    ref2len = {r: l for r, l in zip(sam.references, sam.lengths)}    
    for ref, length in zip(references, lengths):
        if chrs and ref not in chrs:
            if verbose:
                sys.stderr.write(" skipped %s\n"%ref)
            continue
        coverage = np.zeros(length, dtype='uint16')
        for sam in sams:
            coverage = np.max([coverage, get_coverage(sam, ref, mapq, length)], axis=0)
        # get regions with coverage
        covered = np.where(coverage>=mincov)[0]
        for positions in get_consecutive(covered, maxdist):
            if len(positions)<1:
                continue
            s, e = positions[0]+1, positions[-1]+1
            # further split regions for max step size separated by maxdist
            while s < e-step:
                yield ref, s, s+step
                s += step
            yield ref, s, e
    
def load_bed(fname):
    """Return regions from BED file"""
    if os.path.isfile(fname):
        for l in open(fname, "rt"):
            if l.startswith('#') or not l[:-1]:
                continue
            ldata = l[:-1].replace(',','').split('\t')#; print ldata
            if len(ldata) >= 3:
                ref, start, end = ldata[:3]
            else:
                ref, se = ldata[0].split(':')
                start, end = se.split('-')
            start, end = map(int, (start, end))
            yield ref, start, end
    else:
        ref, se = fname.replace(',','').split('\t')[0].split(':')
        start, end = se.split('-')
        start, end = map(int, (start, end))
        yield ref, start, end

def init_args(*args):
    """Share globals with pool of workers"""
    global sams, faidx, dna_sams
    bams, fasta, dna = args
    sams = [pysam.AlignmentFile(bam) for bam in bams]
    dna_sams = [pysam.AlignmentFile(bam) for bam in dna]
    faidx = pysam.FastaFile(fasta)
        
def get_differential_editing(outfn, regionsfn, fasta, dna, rna, minDepth, minDNAfreq, minAltfreq, minAltReads, 
                             stranded, mapq, bcq, threads, barcoded, allpositions, verbose, chrs):
    """Get alternative base coverage and frequency for every bam file"""
    out = gzip.open(outfn, "wt")
    header  = "##\n# REDiscover (ver. %s)\n"%VERSION
    header += "#\n# If you have any suggestions or questions, please use https://github.com/lpryszcz/REDiscover/issues \n##\n"
    header += "# command: %s\n"%" ".join(sys.argv)
    header += "# DNA/RNA alphabet (T=T or U, i=insertion, d=deletion, s=read_start, e=read_end): %s\n"%alphabet
    header += "# BAM files: %s\n" % "\t".join(rna)
    header += "# strands: %s\n##\n" % "\t".join(stranded)
    header += "chromosome\tposition\tvariation\t%s\n"%"\t".join("%s %s"%(bam, base2name[b])
                                                                for bam in rna for b in alphabet)
    out.write(header)

    logger("Genotyping...")
    if regionsfn:
        regions = load_bed(regionsfn)
    else:
        regions = get_covered_regions_per_bam(rna, minDepth, mapq, threads, chrs, verbose) 
      
    p = Pool(threads, initializer=init_args, initargs=(rna, fasta, dna), maxtasksperchild=1000)
    i = 0
    parser = p.imap(worker, ((pos, fasta, minDepth, minDNAfreq, minAltfreq, minAltReads, stranded, mapq, bcq, barcoded, allpositions, verbose) for pos in regions)) #_unordered , chunksize=10
    for i, data in enumerate(parser, 1):
        out.write(data)
            
    logger("%s regions processed"%i)
    out.close()
    p.close()

def compute_strandness(rna, gtf, mapq, threads, verbose, subset=0.01, limit=0.66):
    """Compute strandness"""
    strands = []
    skipUnstranded = False
    toSkip = set()
    for bam, reads, strandness in bam2strandness(rna, gtf, mapq, subset, threads, verbose):
        strand = ""
        if strandness>=limit:
            strand = "secondstrand"
        elif 1.-strandness>=limit:
            strand = "firststrand"
        strands.append(strand)
        if strand and .05<strandness<0.95:
            sys.stderr.write(" [WARNING] %s poorly stranded: %.3f\n"%(bam, strandness))
        if not reads:
            sys.stderr.write(" [WARNING] %s without alignments!\n"%(bam, ))
            toSkip.add(bam)
            
    # notify about mixed strandness
    strandtypes = set(strands)
    if len(strandtypes)>1 and "" in strandtypes:
        sys.stderr.write("[WARNING] Mixing stranded and unstranded BAM files is not allowed: %s!\n"%str(strandtypes))
        skipUnstranded = True
        
    if skipUnstranded or toSkip:
        _rna, _strands = [], []
        for fn, strand in zip(rna, strands):
            if skipUnstranded and not strand or fn in toSkip:
                sys.stderr.write(" %s skipped!\n"%fn)
                continue
            _rna.append(fn)
            _strands.append(strand)
        rna, strands = _rna, _strands
        
    return strands, rna

def is_any_lib_single(bams):
    """Return True if any lib is single-end"""
    for bam in bams:
        r = pysam.Samfile(bam).next()
        if not r.is_paired:
            return True
    
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
    
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument("-o", "--out", required=1, help="output file")
    parser.add_argument("-q", "--mapq", default=3, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("-Q", "--bcq", default=20, type=int, help="basecall quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="number of cores to use [%(default)s]")
    
    refpar = parser.add_mutually_exclusive_group(required=True)
    refpar.add_argument("-d", "--dna", nargs="*", default = [], help="input DNA-Seq BAM file(s)")
    refpar.add_argument("-f", "--fasta", default='', help="reference FASTA file")
    parser.add_argument("-r", "--rna", nargs="+", help="input RNA-Seq BAM file(s)")
    
    strand = parser.add_mutually_exclusive_group()
    strand.add_argument("-g", "--gtf", default="", help="GTF/GFF for auto-detection of strandness")
    strand.add_argument("-u", "--unstranded", action="store_true", help="unstranded RNAseq libraries")
    strand.add_argument("-s", "--stranded", "-fr-secondstrand", action="store_true", 
                        help="stranded RNAseq libraries ie. Illumina or Standard Solid")
    strand.add_argument("-fr-firststrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. dUTP, NSR, NNSR")
    
    parser.add_argument("-b", "--regions", "--bed", help="BED file with regions to genotype")
    parser.add_argument("-c", "--chrs", nargs="*", default=[], help="analyse selected chromosomes [all]")
    parser.add_argument("--allpositions", action='store_true', help="report information for all positions (this will ignore all below filters!) [%(default)s]")
    parser.add_argument("--minDepth", default=5, type=int,  help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minDNAfreq", default=0.99, type=float, help="min frequency for DNA base [%(default)s]")
    parser.add_argument("--minAltfreq", default=0.01, type=float, help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("-a", "--minAltReads", default=3, type=int,  help="min number of reads with alternative allele [%(default)s]")
    parser.add_argument("--barcoded", action="store_true", help="barcoded (UMI) samples")
    
    # get enrichment parameters
    parser.add_argument("--dbSNP", default=[], nargs="+", help="dbSNP file")    
    parser.add_argument("--dist", default=300, type=int, help="distance between SNPs in cluster [%(default)s]")

    parser.add_argument("--nomutual", action="store_true", help="disable mutual info calculation")    
    parser.add_argument("-m", "--mimax", default=0.9, type=float, help="max allowed mutual information [%(default)s]")
    parser.add_argument("--maxcov", default=600, type=int, help="max coverage per sample [%(default)s]")
    parser.add_argument("-n", "--minsamples", nargs="+", default=[1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300],
                        type=int, help="number of samples [%(default)s]")
    parser.add_argument("--nofiltering", action="store_true", help="disable final filtering")    

    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # add extensions
    if not o.out.endswith('.gz'):
        o.out += ".gz"
    # and create outdir if not exists
    if os.path.dirname(o.out) and not os.path.isdir(os.path.dirname(o.out)):
        os.makedirs(os.path.dirname(o.out))
        
    # calculate differential editing if file doesn't exist
    if os.path.isfile(o.out) and open(o.out, "rb").readline():
        sys.stderr.write("Outfile exists or not empty: %s\n"%o.out)
    else:
        # check if all input files exists
        for fn in o.dna+o.rna+[o.fasta]:
            if fn and not os.path.isfile(fn):
                sys.stderr.write("No such file: %s\n"%fn)
                sys.exit(1)

        logger("Indexing bam file(s)...")
        for fn in o.dna+o.rna:
            if not os.path.isfile(fn+".bai"):
                cmd = "samtools index %s"%fn
                if o.verbose:
                    sys.stderr.write(" %s\n"%cmd)
                os.system(cmd)

        # mark stranded protocol
        if o.gtf and not o.unstranded:
            logger("Detecting strandness in BAMs...")
            o.stranded, o.rna = compute_strandness(o.rna, o.gtf, o.mapq, o.threads, o.verbose)
        elif o.stranded:
            o.stranded = ["secondstrand"]*len(o.rna)
        elif o.fr_firststrand:
            o.stranded = ["firststrand"]*len(o.rna)
        else:
            o.stranded = [""]*len(o.rna)
            
        # get differential editing
        get_differential_editing(o.out, o.regions, o.fasta, o.dna, o.rna, o.minDepth, o.minDNAfreq, o.minAltfreq, o.minAltReads, 
                                 o.stranded, o.mapq, o.bcq, o.threads, o.barcoded, o.allpositions, o.verbose, o.chrs)

    # mutual info
    outfn = o.out
    if not o.nomutual:
        logger("Calculating mutual info...")
        outfn = o.out+".pos2mi.gz"
        if os.path.isfile(outfn):
            logger(" %s exists!"%outfn)
        else:
            import pos2mutual
            pos2mutual.pos2mutual(o.out, outfn, o.threads, o.mapq, o.bcq, o.maxcov, o.verbose)

    if not o.nofiltering:
        logger("Filtering...")
        import get_enrichment
        get_enrichment.get_enrichment([outfn,], o.dbSNP, o.minDepth, o.minAltfreq, o.minAltReads, o.minsamples,
                                      o.mimax, o.dist, out=sys.stdout)
        
    logger("Finished!")
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
