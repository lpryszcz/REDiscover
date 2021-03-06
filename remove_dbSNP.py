#!/usr/bin/env python
# Strip snps present in dbSNP
# USAGE: remove_dbSNP.py 00-common_all.vcf.gz file1 [file2 ... fileN]

import os, sys, gzip, pickle

def load_dbSNP(vcf, add_chr=1):
    """Return dictionary of SNP chromosomes and positions"""
    # load pickle
    if os.path.isfile(vcf+".pickle"):
        snps=pickle.load(open(vcf+".pickle"))
        print "Loaded %s SNPs from %s"%(sum(map(len, snps.values())), vcf)
        return snps
    
    snps = {}
    j = k = 0
    for i, l in enumerate(gzip.open(vcf), 1):
        #if i>1000: break
        if l.startswith('#'):
            continue
        j += 1
        chrom, pos, rs, ref, alt = l[:-1].split('\t')[:5]
        #if len(ref)!=len(alt): continue
        # 1 -> chr1
        if add_chr and not chrom.startswith('chr'): 
            chrom = "chr%s"%chrom
        if chrom not in snps:
            snps[chrom] = set()
        snps[chrom].add(int(pos))
        k += 1
    print "Loaded %s out of %s SNPs"%(k, j)
    pickle.dump(snps, open(vcf+".pickle", "w"), 2)
    return snps
    
def parse_editing(fn, snps, out):
    """Ignore positions with known SNPs"""
    j = k = 1
    for l in open(fn):
        if l.startswith("#"):
            out.write(l)
            continue
        if len(l[:-1].split('\t'))<5:
            continue
        j += 1
        chrom, pos, strand, ref, alt = l[:-1].split('\t')[:5]
        if chrom in snps and int(pos) in snps[chrom]:
            continue
        out.write(l)
        k += 1
    return j, k

def main():
    vcf = sys.argv[1]
    fnames = sys.argv[2:]
    
    snps = load_dbSNP(vcf)

    print "#fname\tnot in dbSNP\tall"
    for fn in filter(lambda fn: not fn.endswith('.parsed.txt'), fnames):
        outfn = fn+".parsed.txt"
        # skip if outfn exists and newer than fn
        if os.path.isfile(outfn) and os.stat(fn).st_mtime < os.stat(outfn).st_mtime:
            continue
        out = open(outfn, "w")
        j, k = parse_editing(fn, snps, out)
        print "\t".join(map(str, (fn, k, j)))
    
if __name__=="__main__":
    main()
