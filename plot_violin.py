#!/usr/bin/env python
desc="""Generate violin plot from REDiscover.diff output.

TBA:
- half-violin plot for sense/antisense? https://stackoverflow.com/questions/29776114/half-violin-plot
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Fribourg, 14/08/2017
"""

import os, sys, gzip
import numpy as np
from itertools import izip
from datetime import datetime
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt
from plot_hist import *

def violin_plot(outfn, snps, fnames, id2snp):
    """Generate violin plot"""
    fig, axes = plt.subplots(nrows=1, ncols=12, figsize=(16, .5*len(fnames)+2))

    # process all files
    pos = range(1, len(fnames)+1)
    for ax, _snps, name in izip(axes, snps, id2snp):
        # get ax and set title
        name = "%s\n%s events"%(name[:-1], max([len(x) for x in _snps]))
        ax.set_title(name, fontsize=10)
        
        # violin plot after http://matplotlib.org/examples/statistics/violinplot_demo.html
        ax.violinplot(_snps, pos, points=1000, widths=.5, vert=False, showmeans=True, showextrema=True, showmedians=True, bw_method='silverman')

    # show labels only of first plot
    for i, ax in enumerate(axes.flatten()):
        if i:
            ax.set_yticklabels([])  
        else:
            ax.set_yticks(pos)
            ax.set_yticklabels(fnames)
            
    #fig.suptitle("RNA editing violin plot", fontsize=16)
    fig.subplots_adjust(hspace=.2, left=.15, right=.99, bottom=.05)
    plt.savefig(outfn, dpi=300, transparent=False)
    
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--fnames", nargs="+", help="file(s) to process")
    parser.add_argument("-b", "--bins", default=50, type=int,
                        help="number of bins in histogram [%(default)s]")
    parser.add_argument("-e", "--eid", default=0, type=int,
                        help="element of bam name (after . splitting) [%(default)s]")
    parser.add_argument("-s", "--selected", nargs="+", default=[], help="select samples [all]")
    parser.add_argument("-s2", "--startswith", default='', help="select samples that start with string [all]")
    parser.add_argument("--ext", default="png", choices=['png', 'svg', 'pdf', 'jpg'], 
                        help="figure extension [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))


    # get snp2id
    snp2id, id2snp = {}, []
    s = "+"
    for a in bases:
        for b in bases:
            if a==b: continue
            snp = "%s>%s%s"%(a, b, s)
            snprc = "".join(base2rc[b] for b in snp)
            snp2id[snp] = len(id2snp)
            snp2id[snprc] = len(id2snp)
            id2snp.append(snp)

    # process input files
    for fn in o.fnames:
        snps, names = load_snps(fn, snp2id, o.eid)
        # report stats
        print fn
        print "snp\t" + "\t".join(names)
        oline = "%s\t"*6
        for i, snp in enumerate(id2snp):
            lens = map(len, [filter(lambda x: 0.02<x<0.98, _snps) for _snps in snps[i]])
            print "%s\t"%snp[:-1] + "\t".join(map(str, lens))
            
        outfn = "%s.violin_plot.%s"%(fn, o.ext)
        violin_plot(outfn, snps, names, id2snp)
        sys.stderr.write("[INFO] Violin plot saved as: %s\n"%outfn)

if __name__=="__main__":
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)