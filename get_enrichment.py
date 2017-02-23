#!/usr/bin/env python
# Return enrichment of cannonical editing

# USAGE: ./get_enrichment.py [-s snps.vcf] file1.txt [... fileN.txt]

import os, sys
import numpy as np
from remove_dbSNP import load_dbSNP

bases = "ACGT"
strands = "+-"
bins = 20
base2rc= {"A": "T", "T": "A", "C": "G", "G": "C", ">": ">", "+": "-", "-": "+"}

def txt2changes(editing, handle, snps, minDepth=20, template="%s>%s%s"):
    """Return dictionary of snps and their occurencies"""
    snp2c = {template%(a, b, s): 0 for a in bases for b in bases for s in strands if a!=b}
    for l in handle:
        ldata = l[:-1].split('\t')
        if l.startswith('#') or not l.endswith('\n') or len(ldata)<3:
            continue
        # REDiscover output
        if len(ldata)>4:
            chrom, pos, strand, ref, alt = ldata[:5]
            altcov, altfreq = ldata[-2:]
            altcov, altfreq = int(altcov), float(altfreq)
            # skip if known snp or low cov
            if altcov < minDepth or chrom in snps and int(pos) in snps[chrom]:
                continue
            snp = template%(ref, alt, strand)
        # common.txt
        else:
            chrom, pos, snp = ldata[:3]

        if snp not in snp2c:
            continue
        snp2c[snp] += 1
        # store editing
        k = "%s:%s"%(pos, snp)
        if chrom not in editing:
            editing[chrom] = {k: 1}
        elif k not in editing[chrom]:
            editing[chrom][k] = 1
        else:
            editing[chrom][k] += 1
    return editing, snp2c

def main():
    if sys.argv[1]=="-s":
        snps = load_dbSNP(sys.argv[2])
        fnames = sys.argv[3:]
    else:
        snps = {}
        fnames = sys.argv[1:]

    # process all files
    editing = {}
    snptypes = ["%s>%s"%(a, b) for a in bases for b in bases if a!=b]
    print "#fname\tsites\tstrand enrichment\t"+"\t".join(snptypes)
    for fn in fnames:
        try:
            editing, snp2c = txt2changes(editing, open(fn), snps, minDepth=0)
        except Exception, e:
            sys.stderr.write("[ERROR] Couldn't parse %s with error: %s\n"%(fn, str(e)))
            continue
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
        print "%s\t%s\t%s\t%s"%(fn, total, frac, "\t".join(map(str, freqs)))

    # don't output common txt for single file
    if len(fnames)>1:
        with open("common.txt", "w") as out:
            for chrom in editing:
                for k, v in editing[chrom].iteritems():
                    out.write("%s\t%s\t%s\n"%(chrom, k.replace(':', '\t'), v))
    
if __name__=="__main__":
    main()
    
