#!/usr/bin/env python
# Remove antisense SNPs

import os, sys

def get_coverage(l):
    """Return SNP coverage"""
    return int(l[:-1].split('\t')[7])

def get_position(l):
    """Return chrom, pos and altbase"""
    ldata = l[:-1].split('\t')
    return ldata[:2]+ldata[4:5]

def compare_lines(out, l, pl):
    """Compare current and previous lines and write previous if not antisense SNP"""
    if pl:
        # don't write anything if antisense SNP
        if get_position(l)==get_position(pl):
            # higher coverage will be written in the next round
            if get_coverage(l)<get_coverage(pl):
                return pl
        else:
            out.write(pl)
    return l

for fn in sys.argv[1:]:
    if os.path.isfile(fn+'.bck'):
        print "File exists: %s"%fn
        continue
    # rename file
    os.rename(fn, fn+'.bck')
    with open(fn, "w") as out:
        pl = ''
        for l in open(fn+'.bck'):
            if l.startswith('#'):
                out.write(l)
                continue
            pl = compare_lines(out, l, pl)
        compare_lines(out, l, pl)
            