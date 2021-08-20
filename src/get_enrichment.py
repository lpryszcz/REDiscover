#!/usr/bin/env python3
desc="""Return enrichment of various types of editing
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
base2idx = {b: i for i, b in enumerate(bases)}

def _is_antisense(snps):
    """Return true if snps truly antisense"""
    change1, strand1 = snps[0][-1][:-1], snps[0][-1][-1]
    change2, strand2 = snps[1][-1][:-1], snps[1][-1][-1]
    if change1==change2 and strand1!=strand2:
        return True

def get_counts(l):
    """Return snp counts in all samples"""
    c = 0.0
    ldata = l.split('\t')
    # no mutual info    
    if len(ldata)%2:
        gi = 3
    else:
        gi = 4
    for i in range(gi, len(ldata), 2):
        cov, freq = int(ldata[i]), max(map(float, ldata[i+1].split(";")))
        if cov:
            c += cov * freq
    return c

def clean_cluster(cluster, minFreq=0.501):
    # make sure if clean cluster - remove low freq SNPs
    c = Counter(c[-1] for c in cluster)
    if len(c)>1:
        snp, snpi = c.most_common(1)[0]
        #this causes problems - but why?
        if 1.*snpi/len(cluster) < minFreq:
            #print("lowfreq cluster", c, "%s:%s"%cluster[0][2:4])#, cluster
            return []
        #print("filtering only %s"%snp, c, "%s:%s"%cluster[0][2:4])#cluster
        cluster = filter(lambda x: x[-1]==snp, cluster)
    return cluster
    
def get_clusters(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, mimax, validSNPs, dist=30):
    """SNP cluster generator. Combines cluster"""
    cluster = []
    pchrom = ppos = 0
    for snps in parser_diff(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, mimax):
        # skip mulit-SNP
        if len(snps)>2:
            continue
        # try to rescue correct strand
        elif len(snps)==2:
            # skip if not antisense snps meaning only A>G+ and A>G- allowed, but not A>G+ and A>C-
            if not _is_antisense(snps):
                continue
            #print(snps)
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
        if pchrom!=chrom or pos>ppos+dist:
            if cluster and clean_cluster(cluster):
                yield clean_cluster(cluster)
            cluster = []
            pchrom, ppos = chrom, pos
        # populate cluster
        cluster.append((l, passed, chrom, pos, snp))

    if cluster and clean_cluster(cluster):
        yield clean_cluster(cluster)
    
def txt2changes(editing, handle, dbSNP, minDepth=20, minFreq=0.01, minAltReads=3, \
                minSamples=[3], mimax=1.0, dist=30, template="%s>%s%s"):
    """Return dictionary of snps and their occurencies"""
    outs = {n: gzip.open(handle.name+".n%s.gz"%n, "wt") for n in minSamples}
    snp2c = {n: {template%(a, b, s): 0 for a in bases for b in bases for s in strands if a!=b} for n in minSamples}
    validSNPs = set(snp2c[minSamples[0]].keys())
    for cluster in get_clusters(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, mimax, validSNPs, dist):
        for l, passed, chrom, pos, snp in cluster:
            # store only if enough passed samples
            for n in filter(lambda x: x<=passed, snp2c.keys()): 
                snp2c[n][snp] += 1
                outs[n].write(l)
    # close outs
    for out in outs.values():
        out.close()
    return snp2c, outs

def parser_diff(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, mimax):
    """SNP generator"""
    ppos, snps = '', []
    for l in handle: 
        ldata = l[:-1].split('\t')
        if l.startswith("#") or not l.endswith('\n') or len(ldata)<3: continue
        elif l.startswith('chromosome\tposition'):
            samples = [s.split()[0] for s in ldata[4::8]]#; print(samples); sys.exit()
            for out in outs.values():
                out.write("chrom\tpos\tvar\tmi\t%s\n"%"\t".join("%s cov\t%s freq"%(s, s) for s in samples))
            continue
        chrom, pos, snp = ldata[:3]
        # no mutual info    
        if len(ldata)%2:
            gi = 3
            mi = "-"
        else:
            gi = 4
            mi = float(ldata[3])
            if mimax and mi>mimax: continue
        if ppos!=pos:
            if snps: yield snps
            ppos, snps = pos, []
        # unstranded
        if snp[-1] == ".":
            snp = snp.replace(".", "+")
        if not chrom.startswith('chr'):
            chrom = "chr%s"%chrom
        # skip if present in dbSNP
        if chrom in dbSNP and int(pos) in dbSNP[chrom]: continue
        # get A, C, G & T counts
        ref = snp[0]
        strand = snp[-1]
        refi = base2idx[ref] # this column is reference, thus will be skipped
        sampledata = np.array(list(map(float, ldata[gi:]))).reshape(-1, 8)[:, :4]#; print(sampledata)
        cov = sampledata.sum(axis=1)#; print(cov)
        for alt in snp[2:-1]:
            # skip indels
            if alt not in base2idx: continue 
            alti = base2idx[alt]
            freq = sampledata[:, alti] / cov#; print(chrom, pos, snp, alt, freq)
            # enough depth, frequency and more than 3 reads in alt allele
            passed = sum((cov>=minDepth) & (freq>=minFreq) & (sampledata[:, alti]>=minAltReads))#; print(passed)
            # "\t".join(map(str, ["%.3f"%f for f in freq]))
            sinfo = "\t".join(map(str, ["%i\t%.3f"%(c, f) for c, f in zip(cov, freq)]))
            newl = "%s\t%s\t%s>%s%s\t%s\t%s\n"%(chrom, pos, ref, alt, strand, mi, sinfo)#; print(newl)
            # store
            snps.append((newl, passed, chrom, pos, snp))

    if snps:
        yield snps
    
def txt2changes_noclusters(editing, handle, dbSNP, minDepth=20, minFreq=0.01, minAltReads=3,
                           minSamples=[3], mimax=1.0, dist=0, template="%s>%s%s"):
    """Return dictionary of snps and their occurencies"""
    outs = {n: gzip.open(handle.name+".n%s.gz"%n, "wt") for n in minSamples}
    snp2c = {n: {template%(a, b, s): 0 for a in bases for b in bases for s in strands if a!=b} for n in minSamples}
    for snps in parser_diff(handle, dbSNP, outs, minDepth, minFreq, minAltReads, minSamples, mimax):
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

        # store only if enough passed samples
        for n in filter(lambda x: x<=passed, snp2c.keys()):
            outs[n].write(l)
            # count only standard SNP type
            if snp in snp2c[minSamples[0]]:
                snp2c[n][snp] += 1
    # close outs
    for out in outs.values():
        out.close()
    return snp2c, outs     

def get_enrichment(fnames, dbSNP, minDepth, minAltfreq=.01, minAltReads=3, minsamples=[3], mimax=0.75, dist=0, out=sys.stdout):
    """Filter RNA editing and compute enrichment"""
    snps = {}
    for snpfn in dbSNP:
        _snps = load_dbSNP(snpfn)
        for chrom in _snps:
            if chrom not in snps:
                snps[chrom] = _snps[chrom]
            else:
                snps[chrom] = snps[chrom].union(_snps[chrom])
    print("%s SNPs loaded in total."%sum(map(len, snps.values())))
        
    snptypes = ["%s>%s"%(a, b) for a in bases for b in bases if a!=b]
    print("#fname\tsites\tstrand enrichment\t"+"\t".join(snptypes))
    # process all files
    editing = {}
    for fn in fnames:
        if fn.endswith('.gz'): handle = gzip.open(fn, "rt")
        else: handle = open(fn, "rt")
        
        #minSamplesSNP2c, outs = txt2changes_noclusters(editing, handle, snps, minDepth, minAltfreq, minAltReads, minsamples, mimax, dist)
        minSamplesSNP2c, outs = txt2changes(editing, handle, snps, minDepth, minAltfreq, minAltReads, minsamples, mimax, dist)

        for n in sorted(minsamples):
            snp2c = minSamplesSNP2c[n]
            total = sum(snp2c.values())
            # unlink outfile if not snps
            if not total:
                os.unlink(outs[n].name)
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
            print("%s %s\t%s\t%s\t%s"%(fn, n, total, frac, "\t".join(map(str, freqs))))

def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="files to process")
    parser.add_argument("-s", "--snps", "--dbSNP", default=[], nargs="+", help="dbSNP file")
    parser.add_argument("-d", "--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("-f", "--minAltfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("-a", "--minAltReads",  default=3, type=int,
                        help="min number of reads with alternative base [%(default)s]")
    parser.add_argument("-m", "--mimax", default=1.0, type=float, help="max allowed mutual information [%(default)s]")
    parser.add_argument("-n", "--minsamples", nargs="+", default=[1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300], type=int, help="number of samples [%(default)s]")
    parser.add_argument("--dist", default=300, type=int, help="distance between SNPs in cluster [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    get_enrichment(o.fnames, o.snps, o.minDepth, o.minAltfreq, o.minAltReads, o.minsamples, o.mimax, o.dist)
    
if __name__=="__main__":
    main()
    
