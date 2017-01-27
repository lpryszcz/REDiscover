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

- **fast** & **lightweight**, multi-core support and memory-optimised, so it can be run even on then laptop
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
- Python 2.7+ & numpy `sudo easy_install -U numpy`
- `samtools <http://www.htslib.org/>`_

=====
Usage
=====
REDiscover input consists of **aligned NGS reads** (BAM) from genome(s) and transcriptomes(s).
REDiscover will return a list of putative **RNA editign sites**. 

Parameters
~~~~~~~~~~
Most of REDiscover parameters can be adjusted manually (default values are given in square brackets []):  

  -h, --help            show this help message and exit
  -v, --verbose         verbose
  --version             show program's version number and exit
  -o OUTPUT, --output OUTPUT
                        output stream   [stdout]
  -r RNA, --rna RNA
                        input RNA-Seq BAM file(s)
  -d DNA, --dna DNA
                        input DNA-Seq BAM file(s)
  -f FASTA, --fasta FASTA
                        reference FASTA file
  --minDepth MINDEPTH   minimal depth of coverage [5]
  --minAltReads MINALTREADS
                        minimum no. of reads with alternative base to call RNA editing [3]
  --minRNAfreq MINRNAFREQ
                        min frequency for RNA editing base [0.01]
  --minDNAfreq MINDNAFREQ
                        min frequency for genomic base [0.99]
  --mpileup_opts MPILEUP_OPTS
                        options passed to mpileup         [-I -q 15 -Q 20]
  -t THREADS, --threads THREADS
                        number of cores to use [1] NOT IMPLEMENTED YET!


Test run
~~~~~~~~
To run the test example, execute:

.. code-block:: bash

    ###
    # RNAseq + DNAseq
    ./REDiscover -r test/RNA.bam -d test/DNA.bam 
    
    # filter by min. frequency and cluster optionally
    
    
    ###
    # RNAseq alone (high false positive expected!)
    ./REDiscover -r test/RNA.bam -f test/ref.fa
    
    # discard known SNPs ie. using dbSNP


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
Pryszcz LP, Bochtler M, Winata CL. (In preparation) Detection of RNA editing from NGS. 
