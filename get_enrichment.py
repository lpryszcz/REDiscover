#!/usr/bin/env python
desc="""Return enrichment of various types of editing

TBD:
- skip mixed regions or low freq SNPs in clusters ie <30bp apart
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw/Bratislava/Fribourg, 21/07/2015
"""

import os, sys, gzip
import numpy as np
from remove_dbSNP import load_dbSNP
from collections import Counter

bases = "ACGT"
strands = "+-"
bins = 20
base2rc= {"A": "T", "T": "A", "C": "G", "G": "C", ">": ">", "+": "-", "-": "+"}

def parser_diff(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples):
    """SNP generator"""
    ppos = 0
    snps = []
    for l in handle: 
        ldata = l.replace('\t\t','\t')[:-1].split('\t')
        if l.startswith('#') or not l.endswith('\n') or len(ldata)<3:
            for out in outs.itervalues():
                out.write(l)
            continue
        chrom, pos, snp = ldata[:3]
        if ppos!=pos:
            if snps:
                yield snps
            ppos = pos
            snps = []    
        # unstranded
        if snp[-1] == ".":
            snp = snp.replace(".", "+")
        if not chrom.startswith('chr'):
            chrom = "chr%s"%chrom
        # skip if present in dbSNP
        if chrom in dbSNP and int(pos) in dbSNP[chrom]:
            continue
        sampledata = np.array(map(float, ldata[3:])).reshape(len(ldata[3:])/2, 2)
        # enough depth, frequency and more than 3 reads in alt allele    
        passed = sum((sampledata[:, 0]>=minDepth) & (sampledata[:, 1]>=minFreq) \
                     & (sampledata[:, 0]*sampledata[:, 1]>=minAltReads))
        # store
        snps.append((l, passed, chrom, pos, snp))
    if snps:
        yield snps

def _is_antisense(snps):
    """Return true if snps truly antisense"""
    change1, strand1 = snps[0][-1][:-1], snps[0][-1][-1]
    change2, strand2 = snps[1][-1][:-1], snps[1][-1][-1]
    if change1==change2 and strand1!=strand2:
        #print "antisense"
        return True

def get_counts(l):
    """Return snp counts in all samples"""
    c = 0.0
    ldata = l.split('\t')
    for i in range(3, len(ldata), 2):
        cov, freq = int(ldata[i]), float(ldata[i+1])
        if cov:
            c += cov * freq
    return c

def clean_cluster(cluster, minFreq=0.67):
    # make sure if clean cluster - remove low freq SNPs
    c = Counter(c[-1] for c in cluster)
    if len(c)>1:
        snp, snpi = c.most_common(1)[0]
        if 1.*snpi/len(cluster) < minFreq:
            #print "lowfreq cluster", c, "%s:%s"%cluster[0][2:4]#, cluster
            return []
        #print "filtering only %s"%snp, c, "%s:%s"%cluster[0][2:4]#cluster
        cluster = filter(lambda x: x[-1]==snp, cluster)
    return cluster
    
def get_clusters(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, validSNPs, step=30):
    """SNP cluster generator. Combines cluster"""
    cluster = []
    pchrom = ppos = 0
    for snps in parser_diff(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples):
        # skip mulit-SNP
        if len(snps)>2:
            continue
        # try to rescue correct strand
        elif len(snps)==2:
            # skip if not antisense snps meaning only A>G+ and A>G- allowed, but not A>G+ and A>C-
            if not _is_antisense(snps):
                continue
            counts = [get_counts(l) for l, passed, chrom, pos, snp in snps]
            maxi = np.argmax(counts)#; print counts, maxi
            l, passed, chrom, pos, snp = snps[maxi]
        else:
            l, passed, chrom, pos, snp = snps[0]
            
        # skip if unexpected editing type
        if snp not in validSNPs:
            continue
        # report cluster
        pos = int(pos)    
        if pchrom!=chrom or pos>ppos+step:
            if cluster and clean_cluster(cluster):
                yield clean_cluster(cluster)
            cluster = []
            pchrom, ppos = chrom, pos
        # populate cluster
        cluster.append((l, passed, chrom, pos, snp))

    if cluster and clean_cluster(cluster):
        yield clean_cluster(cluster)
    
def txt2changes_clusters(editing, handle, dbSNP, minDepth=20, minFreq=0.01, minAltReads=3, minSamples=[3], template="%s>%s%s"):
    """Return dictionary of snps and their occurencies"""
    outs = {n: gzip.open(handle.name+".n%s.gz"%n, "w") for n in minSamples}
    snp2c = {n: {template%(a, b, s): 0 for a in bases for b in bases for s in strands if a!=b} for n in minSamples}
    validSNPs = set(snp2c[minSamples[0]].keys())
    for cluster in get_clusters(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, validSNPs):
        for l, passed, chrom, pos, snp in cluster:
            # store only if enough passed samples
            for n in filter(lambda x: x<=passed, snp2c.keys()): 
                snp2c[n][snp] += 1
                outs[n].write(l)
    # close outs
    for out in outs.itervalues():
        out.close()
    return snp2c
    
def txt2changes_noclusters(editing, handle, dbSNP, minDepth=20, minFreq=0.01, minAltReads=3, minSamples=[3], template="%s>%s%s"):
    """Return dictionary of snps and their occurencies"""
    outs = {n: gzip.open(handle.name+".n%s.gz"%n, "w") for n in minSamples}
    snp2c = {n: {template%(a, b, s): 0 for a in bases for b in bases for s in strands if a!=b} for n in minSamples}
    for snps in parser_diff(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples):
        # skip mulit-SNP
        if len(snps)>2:
            continue
        # try to rescue correct strand
        elif len(snps)==2:
            # skip if not antisense snps meaning only A>G+ and A>G- allowed, but not A>G+ and A>C-
            if not _is_antisense(snps):
                continue
            counts = [get_counts(l) for l, passed, chrom, pos, snp in snps]
            maxi = np.argmax(counts)#; print counts, maxi
            l, passed, chrom, pos, snp = snps[maxi]
        else:
            l, passed, chrom, pos, snp = snps[0]
            
        # skip if unexpected editing type
        if snp not in snp2c[minSamples[0]]:
            continue

        # store only if enough passed samples
        for n in filter(lambda x: x<=passed, snp2c.keys()): 
            snp2c[n][snp] += 1
            outs[n].write(l)
    # close outs
    for out in outs.itervalues():
        out.close()
    return snp2c        

def get_enrichment(fnames, snps, minDepth, minAltfreq, minAltReads, minsamples, snptypes, noclusters=0, out=sys.stdout):
    """Filter RNA editing and compute enrichment"""
    if noclusters:
        txt2changes = txt2changes_noclusters
    else:
        txt2changes = txt2changes_clusters
    # process all files
    editing = {}
    for fn in fnames:
        if fn.endswith('.gz'):
            handle = gzip.open(fn)            
        else:
            handle = open(fn)
        
        minSamplesSNP2c = txt2changes(editing, handle, snps, minDepth, minAltfreq, minAltReads, minsamples)

        for n in sorted(minsamples):
            snp2c = minSamplesSNP2c[n]
            total = sum(snp2c.itervalues())
            if not total:
                sys.stderr.write("[WARNING] No editing in %s\n"%(fn, ))
                continue
            # get freqs and strand enrichment
            strands = []
            freqs = []
            for snp in sorted(filter(lambda x: x[-1]=="+", snp2c)):
                snprc = "".join(base2rc[b] for b in snp)
                freq = 1.*(snp2c[snp]+snp2c[snprc])/total
                freqs.append(freq)
                strands.append((snp2c[snp], snp2c[snp[:-1]+"-"]))
            strands = np.array(strands, dtype="float32")
            frac = sum(strands.max(axis=1))/strands.sum()
            print "%s %s\t%s\t%s\t%s"%(fn, n, total, frac, "\t".join(map(str, freqs)))

def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="files to process")
    parser.add_argument("-s", "--snps", default=[], nargs="+", help="dbSNP file")
    parser.add_argument("-d", "--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("-f", "--minAltfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("-a", "--minAltReads",  default=3, type=int,
                        help="min number of reads with alternative base [%(default)s]")
    parser.add_argument("-n", "--minsamples", nargs="+", default=[1], type=int, help="number of samples [%(default)s]")
    parser.add_argument("-c", "--noclusters", action="store_true", help="don't use clusters filtering")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    snps = {}
    for snpfn in o.snps:
        _snps = load_dbSNP(snpfn)
        for chrom in _snps:
            if chrom not in snps:
                snps[chrom] = _snps[chrom]
            else:
                snps[chrom] = snps[chrom].union(_snps[chrom])
    print "%s SNPs loaded in total."%sum(map(len, snps.itervalues()))
        
    snptypes = ["%s>%s"%(a, b) for a in bases for b in bases if a!=b]
    print "#fname\tsites\tstrand enrichment\t"+"\t".join(snptypes)
    get_enrichment(o.fnames, snps, o.minDepth, o.minAltfreq, o.minAltReads, o.minsamples, snptypes, o.noclusters)
    
if __name__=="__main__":
    main()
    
