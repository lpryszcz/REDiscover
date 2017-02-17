#!/usr/bin/env python
# Strip snps present in dbSNP
# USAGE: ./remove_dbSNP.py 00-common_all.vcf.gz file1 [file2 ... fileN]

import os, sys, gzip, pickle

vcf = sys.argv[1]
fnames = sys.argv[2:]

def load_dbSNP(vcf, add_chr=1):
  """Return dictionary of SNP chromosomes and positions"""
  # load pickle
  if os.path.isfile(vcf+".pickle"):
    snps=pickle.load(open(vcf+".pickle"))
    print "Loaded %s snps"%sum(map(len, snps.values()))
    return snps
  
  snps = {}
  j = k = 0
  for i, l in enumerate(gzip.open(vcf), 1):
    #if i>1000: break
    if l.startswith('#'):
      continue
    j += 1
    chrom, pos, rs, ref, alt = l[:-1].split('\t')[:5]
    if len(ref)!=len(alt):
      continue
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
  
snps = load_dbSNP(vcf)

for i, fn in enumerate(fnames, 1):
  outfn = fn+".parsed.txt"
  #if os.path.isfile(outfn):
  #  continue
  out = open(outfn, "w")
  j, k = parse_editing(fn, snps, out)
  print fn, k, j
  
