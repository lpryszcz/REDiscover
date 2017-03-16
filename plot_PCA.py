


def fn2tissue(fn): return os.path.basename(fn).split('.')[1]
def fn2donor(fn): return os.path.basename(fn).split('.')[0]
def fn2replica(fn): return '.'.join(os.path.basename(fn).split('.')[:2])
    
def plot_PCA(tmpfn, outbase, bams, minDepth, minAltfreq, verbose, n=3, frac=0.33):
    """Plot PCA"""
    # load array
    d = np.loadtxt(tmpfn, usecols=range(3,3+len(bams)*2))
    # reshape bam x snps x (cov, freq)
    d2 = d.reshape(len(bams), d.shape[0], 2)
    selected_samples = np.sum(d2[:,:,0]>=minDepth, axis=1)>=frac*d2.shape[1]
    d2 = d2[selected_samples]
    selected_snps = np.sum(d2[:,:,0]>=minDepth, axis=0)>=frac*d2.shape[0]
    X = d2[:, selected_snps, 1]>=minAltfreq
    
    classes = []
    tissues = {}
    fn2class = []
    for i, fn in enumerate(bams):
        # skip if not in selected
        if not selected_samples[i]:
            sys.stderr.write(" %s skipped.\n"%fn)
            continue
        tissue = fn2tissue(fn) # fn2donor(fn) # fn2replica(fn) #
        if tissue not in tissues:
            tissues[tissue] = len(tissues)
        classes.append(tissues[tissue])
    y = np.array(classes)#, dtype=float)
    print d2.shape, X.shape, len(y), len(tissues), tissues
        
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
    parser.add_argument("-o", "--outbase", required=1, help="output file")
    parser.add_argument("-b", "--regions", "--bed", help="BED file with regions to genotype")
    refpar = parser.add_mutually_exclusive_group(required=True)
    refpar.add_argument("-d", "--dna", nargs="*", default = [],  help="input DNA-Seq BAM file(s)")
    refpar.add_argument("-f", "--fasta", default = None,  help="reference FASTA file")
    parser.add_argument("-r", "--rna", nargs="+",  help="input RNA-Seq BAM file(s)")
    parser.add_argument("-s", "--stranded", "-fr-secondstrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. Illumina or Standard Solid")
    parser.add_argument("-fr-firststrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. dUTP, NSR, NNSR")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minAltfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("--bcq", default=20, type=int, help="basecall quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="number of cores to use [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
     
    # mark stranded protocol
    if o.fr_firststrand:
        o.stranded = "firststrand"
    
    # check if all input files exists
    for fn in o.rna:
        if not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)

    # calculate differential editing if file doesn't exist
    tmpfn = o.outbase+".diff.tsv"
    if not os.path.isfile(tmpfn) or not open(tmpfn).readline():
        logger("Calculating differential editing...")
        get_differential_editing(tmpfn, o.regions, o.fasta, o.rna, o.minDepth, o.minAltfreq, o.stranded, o.mapq, o.bcq,
                                 o.threads, o.verbose)
    logger("Done!")
    return
    logger("Processing...")
    plot_PCA(tmpfn, o.outbase, o.rna, o.minDepth, o.minAltfreq, o.verbose)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
