#!/usr/bin/env python
# Return histogram of each type of editing at min depth of coverage 20X

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
  
def txt2changes(handle, snps, minDepth=20, template="%s>%s%s"):
  """Return dictionary of snps and their occurencies"""
  snp2freq = {template%(a, b, s): [] for a in bases for b in bases 
                                     for s in strands if a!=b}
  for l in handle:
    if l.startswith('#'):
      continue
    ldata = l[:-1].split('\t')
    try:
      chrom, pos, strand, ref, alt = ldata[:5]
      altcov, altfreq = ldata[-2:]
      altcov, altfreq = int(altcov), float(altfreq)
    except:
      sys.stderr.write("[WARNING] Wrong line %s\n"%str(ldata))
      break
    # skip if known snp or low cov
    if altcov < minDepth or chrom in snps and int(pos) in snps[chrom]:
      continue
    snp = template%(ref, alt, strand)
    if snp not in snp2freq:
      continue
    snp2freq[snp].append(altfreq)
  return snp2freq

def plot_hist(bins, fn, outfn, snp2freq, colors, log=0):
  """Plot hist"""
  fig, subplots = plt.subplots(4, 3, sharey='all', sharex='all', figsize=(15, 15))
  fig.suptitle(fn, size=20)
  for ii, snp in enumerate(sorted(filter(lambda x: x[-1]=="+", snp2freq))):
    ax = subplots[ii/3][ii%3] 
    snprc = "".join(base2rc[b] for b in snp)
    data = [[], []]
    if snp in snp2freq:
      data[0] = snp2freq[snprc]
    if snprc in snp2freq:
      data[1] = snp2freq[snprc]
    n, bins, patches = ax.hist(data, bins, normed=0, stacked=True, color=colors, label=["+", "-"])
    ax.set_title(snp[:-1])
    ax.grid(True)
  # make y axis log
  if log:  
    plt.yscale('log')
  plt.savefig(outfn, dpi=300, transparent=False) #orientation='landscape', 
  
  
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
colors = ['red', 'red'] #["grey", "grey"] #['blue', 'red',]
for i, fn in enumerate(fnames, 1):
  # skip if outfn exists and newer than fn
  outfn = "%s.hist.png"%fn
  if os.path.isfile(outfn) and os.stat(fn).st_mtime < os.stat(outfn).st_mtime:
    continue
  print i, fn
  # get snps  
  snp2freq = txt2changes(open(fn), snps)
  if not sum(len(x) for x in snp2freq.values()):
    sys.stderr.write(" No editing found!\n")
    continue
  # plot hist
  plot_hist(bins, fn, outfn, snp2freq, ['blue', 'grey'], log=0) 
  plot_hist(bins, fn, outfn.replace('.png', '.log.png'), snp2freq, ['blue', 'blue'], log=1)
  
