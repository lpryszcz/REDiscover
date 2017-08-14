REDiscover
==========
Tool for RNA editing discovery from NGS data.

.. contents:: Table of Contents

===========
Methodology
===========
REDiscover reports differences between transcriptome and underlying genome, these are putative RNA editing sites.
To achieve that, genome and transcriptome are genotyped simultanously and basecalls are compared.

REDiscover is:

- **easy-to-use** - the programme will auto-detect and estimate all necessary parameters ie. strandness of your library.
- **fast** & **lightweight**, multi-core support and memory-optimised, so it can be run even on the laptop
- **flexible** toward many sequencing technologies and experimental designs ie. stranded and unstranded RNA-Seq, multiple genomes and/or transcriptomes are accepted as input
- **reliable** - the tools was tested extensively on vertebrates *D. rerio*

.. image:: /docs/flowchart.png
           :width: 600 px 
           :align: right
           :scale: 75

                   
By default, REDiscover filters:

- QC failed reads
- duplicates
- reads with mapping quality (mapQ) below 15 

REDiscover reports only regions fulfilling several stringency criteria:

- depth-of-coverage
- mean basecall quality

Finally, reads with basecall quality below 20 (0.01 probability of error) for given positiong are ignored. 

.. [//]: # "For more information have a look at the [poster](/docs/poster.pdf) or [manuscript](/docs/manuscript.pdf)."

=============
Prerequisites
=============
- Python 2.7+
- `samtools <http://www.htslib.org/>`_
- `pysam <https://github.com/pysam-developers/pysam>`_
- `FastaIndex <https://github.com/lpryszcz/FastaIndex>`_
- `numpy <http://www.numpy.org/>`_
- `matplotlib <http://matplotlib.org/>`_ for plotting

All above can be easily installed with `bioconda <https://bioconda.github.io/>`_:

.. code-block:: bash

conda install samtools pysam FastaIndex numpy matplotlib


=====
Usage
=====
REDiscover input consists of **aligned NGS reads** (BAM) from genome(s) and transcriptomes(s).
REDiscover will return a list of putative **RNA editign sites** and their depth of coverage
and frequency across samples. 
Note, you can run REDiscover with RNA-seq reads alone, then you need to provide reference FastA. 
REDiscover will detect strandness of your library if you provide it with exon annotation (GTF or GFF).
Note, mixing of stranded and unstranded libraries is not allowed! 

Parameters
~~~~~~~~~~
Most of REDiscover parameters can be adjusted manually (default values are given in square brackets []):  

- General options
  
  -h, --help            show this help message and exit
  -v, --verbose         verbose
  --version             show program's version number and exit
  -o OUT, --out OUT     output file
  -q MAPQ, --mapq MAPQ  mapping quality [3]
  -Q BCQ, --bcq BCQ     basecall quality [20]
  -t THREADS, --threads THREADS
                        number of cores to use [4]

- Reference genome (BAM or FastA)
  
  -d DNA, --dna DNA     input DNA-Seq BAM file(s)
  -f FASTA, --fasta FASTA
                        reference FASTA file
                        
- Aligned RNA-seq reads & strandness information
  
  -r RNA, --rna RNA     input RNA-Seq BAM file(s)  
  -g GTF, --gtf GTF     GTF/GFF for auto-detection of strandness
  -u, --unstranded      unstranded RNAseq libraries
  -s, --stranded, -fr-secondstrand
                        stranded RNAseq libraries ie. Illumina or Standard Solid
  -fr-firststrand       stranded RNAseq libraries ie. dUTP, NSR, NNSR

- Analyse only subset of regions
  
  -b REGIONS, --regions REGIONS, --bed REGIONS
                        BED file with regions to genotype
  -c CHRS, --chrs CHRS  analyse only sublset of chromosomes [all]

- Filtering
  
  --minDepth MINDEPTH   minimal depth of coverage [5]
  --minDNAfreq MINDNAFREQ
                        min frequency for DNA base [0.99]
  --minAltfreq MINALTFREQ
                        min frequency for RNA editing base [0.01]
  -m MAXSTRANDBIAS, --maxStrandBias MAXSTRANDBIAS
                        max allowed strand bias [0.1]
  -a, --advancedFiltering
                        enable advanced filtering (slightly more accurate, but much slower)
  --dbSNP               dbSNP file
  --dist DIST           distance between SNPs in cluster [300]



                        
Test run
~~~~~~~~
To run the test example, first download & unpack `the test dataset <http://zdglab.iimcb.gov.pl/lpryszcz/REDiscover/test.tgz>`_:

.. code-block:: bash

   wget http://zdglab.iimcb.gov.pl/lpryszcz/REDiscover/test.tgz
   tar xpfvz test.tgz


Then execute `REDiscover.diff`:

.. code-block:: bash

   # discover editing in RNA-seq samples (*.bam) without reference sequencing (ref.fa needed)
   ~/src/REDiscover/REDiscover.diff -f test/ref.fa -r test/star/*.bam -o test/editing.gz

   # discover editing in RNA-seq samples (*.bam) with reference sequencing (ref*.bam needed)
   ~/src/REDiscover/REDiscover.diff -d test/ref*.bam -r test/star/*.bam -o test/editing.ref.gz

   # if you want to ignore dbSNP sites, just add `--dbSNP snps.vcf.gz` to above commands
   # or recompute only last step using `./get_enrichment.py`
   ## you can alter also `--minDepth`, `--minAltfreq` and many more...
   ~/src/REDiscover/get_enrichment.py -i test/editing.gz --dbSNP snps.vcf.gz
   
   # violin plots for editing sites present in at least 2 samples
   ~/src/REDiscover/plot_violin.py -i test/editing.gz.n2.gz

   # histograms for editing sites present in at least 5 samples
   ~/src/REDiscover/plot_hist.py -i test/editing.gz.n5.gz

For more details have a look in `test directory </test>`_. 

=====
Tools
=====
Along with REDiscover, we provide a bunch of usefull tools for characterisation of RNA editing.
More details about these can be find in `tools directory </tools>`_. 

===
FAQ
===


========
Citation
========
Pryszcz LP, Bochtler M, Winata CL. (In preparation) REDiscover: Robust & efficient detection of RNA editing from large NGS datasets. 
