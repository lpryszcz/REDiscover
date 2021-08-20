#!/usr/bin/env python
desc="""Identify differential RNA editing sites from multiple RNAseq experiments (.bam). 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw/Bratislava/Fribourg, 21/07/2015
"""

import os, sys
from datetime import datetime
import numpy as np

'''
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition
from collections import Counter
#'''

def fn2tissue(fn): return os.path.basename(fn).split('.')[1]
def fn2donor(fn): return os.path.basename(fn).split('.')[0]
def fn2replica(fn): return '.'.join(os.path.basename(fn).split('.')[:2])
    
def plot_PCA(infile, outbase, minDepth, minAltfreq, verbose, n=3, frac=0.33):
    """Plot PCA"""
    f = open(infile)
    cmd, fields = f.readline(), f.readline()[:-1].split('\t')
    bams = [fields[i].split()[0] for i in range(4, len(fields), 2)]; print bams[:3]
    # load array
    d = np.loadtxt(infile, usecols=range(3,3+len(bams)*2))
    # reshape bam x snps x (cov, freq)
    d2 = d.reshape(len(bams), d.shape[0], 2)
    selected_samples = np.sum(d2[:,:,0]>=minDepth, axis=1)>=frac*d2.shape[1]
    d2 = d2[selected_samples]
    selected_snps = np.sum(d2[:,:,0]>=minDepth, axis=0)>=frac*d2.shape[0]
    X = d2[:, selected_snps, 1]>=minAltfreq
    
    classes, fn2class = [], []
    tissues, t2c = {}, {}
    for i, fn in enumerate(bams):
        # skip if not in selected
        if not selected_samples[i]:
            sys.stderr.write(" %s skipped.\n"%fn)
            continue
        tissue = fn2tissue(fn) # fn2donor(fn) # fn2replica(fn) #
        if tissue not in tissues:
            tissues[tissue] = len(tissues)
            t2c[tissue] = 1
        else:
            t2c[tissue] += 1
        classes.append(tissues[tissue])
            
    y = np.array(classes)#, dtype=float)
    print d2.shape, X.shape, len(y), len(tissues), tissues

    # print stats
    print len(t2c), sorted(t2c.iteritems(), key=lambda x: x[1], reverse=1); return
        
    # http://scikit-learn.org/stable/auto_examples/decomposition/plot_pca_iris.html    
    fig = plt.figure(1, figsize=(4, 3))
    ax = plt.subplot(111) # Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134) #

    plt.cla()
    pca = decomposition.PCA(n_components=n)
    pca.fit(X)
    X = pca.transform(X)
    print X.shape, pca.explained_variance_ratio_
    
    colors = plt.cm.Paired(np.linspace(0, 1, len(tissues)))
    markers = 'xopsv<^>*h'
    for tissue, label in tissues.iteritems():
        ax.scatter(X[y==label, 0], X[y==label, 1], #X[y==label, 2],
                   label=tissue, color=colors[label], marker=markers[label%len(markers)])

    ax.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8)#, bbox_to_anchor=(0, 0))
    plt.show()
  
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--infile", required=1, help="input file")
    parser.add_argument("-o", "--outbase", default=sys.stdout, help="output file [stdout]")
    parser.add_argument("-b", "--regions", "--bed", help="BED file with regions to genotype")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minAltfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
     
    plot_PCA(o.infile, o.outbase, o.minDepth, o.minAltfreq, o.verbose)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
