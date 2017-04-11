#!/usr/bin/env python
desc="""Return enrichment of various types of editing
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw/Bratislava/Fribourg, 21/07/2015
"""

import os, sys, gzip
import numpy as np
from remove_dbSNP import load_dbSNP

bases = "ACGT"
strands = "+-"
bins = 20
base2rc= {"A": "T", "T": "A", "C": "G", "G": "C", ">": ">", "+": "-", "-": "+"}

def txt2changes(editing, handle, snps, minDepth=20, minFreq=0.01, minSamples=3, template="%s>%s%s"):
    """Return dictionary of snps and their occurencies"""
    out = gzip.open(handle.name+".n%s.gz"%minSamples, "w")
    snp2c = {template%(a, b, s): 0 for a in bases for b in bases for s in strands if a!=b}
    for l in handle:
        ldata = l.replace('\t\t','\t')[:-1].split('\t')
        if l.startswith('#') or not l.endswith('\n') or len(ldata)<3:
            continue
        # REDiscover output
        if ldata[2] in '+-.': #len(ldata)>4
            chrom, pos, strand, ref, alt = ldata[:5]
            if not chrom.startswith('chr'):
                chrom = "chr%s"%chrom
            altcov, altfreq = ldata[-2:]
            altcov, altfreq = int(altcov), float(altfreq)
            # skip if known snp or low cov
            if altcov < minDepth or chrom in snps and int(pos) in snps[chrom]:
                continue
            snp = template%(ref, alt, strand)
        # REDiscover.diff2 output # get
        elif len(ldata)>4:
            chrom, pos, snp = ldata[:3]
            if not chrom.startswith('chr'):
                chrom = "chr%s"%chrom
            # skip if present in dbSNP
            if chrom in snps and int(pos) in snps[chrom]:
                continue
            sampledata = np.array(map(float, ldata[3:])).reshape(len(ldata[3:])/2, 2)
            passed = sum((sampledata[:, 0]>=minDepth) & (sampledata[:, 1]>=minFreq))# & (sampledata[:, 1]<1.0))
            #print passed, ldata[3:]
            if passed<minSamples: continue
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
        out.write(l)
    out.close()
    return editing, snp2c

def get_enrichment(fnames, snps, minDepth, minAltfreq, minsamples, snptypes, out=sys.stdout):
    # process all files
    editing = {}
    for fn in fnames:
        if fn.endswith('.gz'):
            handle = gzip.open(fn)            
        else:
            handle = open(fn)
        try:
            editing, snp2c = txt2changes(editing, handle, snps, minDepth, minAltfreq, minsamples)
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
        print "%s %s\t%s\t%s\t%s"%(fn, minsamples, total, frac, "\t".join(map(str, freqs)))

    # don't output common txt for single file
    if len(fnames)>1:
        with open("common.txt", "w") as out:
            for chrom in editing:
                for k, v in editing[chrom].iteritems():
                    out.write("%s\t%s\t%s\n"%(chrom, k.replace(':', '\t'), v))
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="files to preocess")
    parser.add_argument("-s", "--snps", default=[], nargs="+", help="dbSNP file")
    parser.add_argument("-d", "--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("-f", "--minAltfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("-n", "--minsamples", nargs="+", default=[1], type=int, help="number of samples [%(default)s]")
    
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
                snps[chrom].union(_snps[chrom])
        
    snptypes = ["%s>%s"%(a, b) for a in bases for b in bases if a!=b]
    print "#fname\tsites\tstrand enrichment\t"+"\t".join(snptypes)
    for n in o.minsamples:
        get_enrichment(o.fnames, snps, o.minDepth, o.minAltfreq, n, snptypes)
    
if __name__=="__main__":
    main()
    
