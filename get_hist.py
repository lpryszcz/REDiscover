#!/usr/bin/env python
# Return histogram of each type of editing

# USAGE: ./get_hist.py [-s snps.vcf] file1.txt [... fileN.txt]

import os, sys
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt

def load_snps(handle):
  """Return dictionary of chromosomes and set of snps positions"""
  snps = {}
  for l in handle:
    if l.startswith('#'):
      continue  
    chrom, pos = l.split('\t')[:5]
    if chrom not in snps:
      snps[chrom] = {}
    snps[chrom].add(int(pos))
  return snps
  
def txt2changes(handle, snps, template="%s>%s%s"):
  """Return dictionary of snps and their occurencies"""
  snp2freq = {template%(a, b, s): [] for a in bases for b in bases 
                                     for s in strands if a!=b}
  for l in handle:
    if l.startswith('#'):
      continue
    try:
      chrom, pos, strand, ref, alt, refcov, refQ, reffreq, altcov, altQ, altfreq = l[:-1].split('\t')
      altfreq = float(altfreq)
    except:
      sys.stderr.write("[WARNING] Wrong line %s\n"%str(l[:-1].split('\t')))
      break
    # skip if known snp
    if chrom in snps and int(pos) in snps[chrom]:
      continue
    snp = template%(ref, alt, strand)
    if snp not in snp2freq:
      continue
    snp2freq[snp].append(altfreq)
  return snp2freq


bases = "ACGT"
strands = "+-"
bins = 20
base2rc= {"A": "T", "T": "A", "C": "G", "G": "C", ">": ">", "+": "-", "-": "+"}

if sys.argv[1]=="-s ":
  snps = load_snps(open(sys.argv[2]))
  fnames = sys.argv[3:]
else:
  snps = {}
  fnames = sys.argv[1:]
  
# process all files
for i, fn in enumerate(fnames, 1):
  fig = plt.figure(figsize=(15, 15)) # figsize=(24,16)) # figsize=(11.69,8.27)) #
  fig.suptitle(fn, size=20)
  print i, fn
  outfn = "%s.hist.png"%fn
  if os.path.isfile(outfn):
    continue
  # get snps  
  snp2freq = txt2changes(open(fn), snps)
  if not sum(len(x) for x in snp2freq.values()):
    sys.stderr.write(" No editing found!\n")
    continue
  # plot hist
  for ii, snp in enumerate(sorted(filter(lambda x: x[-1]=="+", snp2freq)), 1):
    ax = fig.add_subplot(4, 3, ii)
    snprc = "".join(base2rc[b] for b in snp)#; print snp, snprc
    data = [[], []]
    if snp in snp2freq:
      data[0] = snp2freq[snprc]
    if snprc in snp2freq:
      data[1] = snp2freq[snprc]
    n, bins, patches = ax.hist(data, bins, normed=0, stacked=True, color=['blue', 'red',], label=["+", "-"])
    ax.set_title(snp[:-1])
    ax.grid(True)
  plt.savefig(outfn, dpi=300, transparent=False) #orientation='landscape', 
  
  
