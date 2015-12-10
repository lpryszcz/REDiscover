### Table of Contents
- **[REDiscover](#rediscover)**  
  - **[Methodology](#methodology)**  
  - **[Prerequisites](#prerequisites)**  
  - **[Usage](#usage)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[FAQ](#faq)**  
  - **[Citation](#citation)**  

# REDiscover
Tool for **R**NA **e**diting **discover**y from NGS data. 

## Methodology

REDiscover reports differences between transcriptome and underlying genome, these are putative RNA editing sites.
To achieve that, genome and transcriptome are genotyped simultanously and basecalls are compared.

REDiscover is:
- **fast** & **lightweight**, multi-core support and memory-optimised, 
so it can be run even on then laptop
- **flexible** toward many sequencing technologies and experimental designs ie. stranded and unstranded RNA-Seq, multiple genomes and/or transcriptomes are accepted as input
- **reliable** - the tools was tested extensively on vertebrates (D. rerio) 


By default, REDiscover filters 
- QC failed reads
- duplicates
- reads with mapping quality (mapQ) below 15 

REDiscover reports only regions fulfilling several stringency criteria:
- depth-of-coverage
- mean basecall quality

Finally, reads with basecall quality below 20 (0.01 probability of error) for given positiong are ignored. 

For more information have a look at the [poster](/docs/poster.pdf) or [manuscript](/docs/manuscript.pdf).

![Flowchart](/docs/flowchart.png)

## Prerequisites
- Python 2.7+ & numpy `sudo easy_install -U numpy`
- [samtools](http://www.htslib.org/)

## Usage
REDiscover input consists of **aligned NGS reads** (BAM) from genome(s) and transcriptomes(s).
REDiscover will return a list of putative **RNA editign sites**. In addition, ... 

### Parameters
Most of REDiscover parameteers can be adjusted manually (default values are given in square brackets []):  
```
  -h, --help            show this help message and exit
  -v, --verbose         verbose
  --version             show program's version number and exit
  -o OUTPUT, --output OUTPUT
                        output stream   [stdout]
  -r RNA [RNA ...], --rna RNA [RNA ...]
                        input RNA-Seq BAM file(s)
  -d [DNA [DNA ...]], --dna [DNA [DNA ...]]
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
```

### Test run
To run the test example, execute: 
```bash
./REDiscover -d test/DNA.bam -r test/RNA.bam
```

For more details have a look in [test directory](/test). 

## FAQ

## Citation
Pryszcz LP, Botchler M, Winata CL (In preparation) Detection of RNA editing from NGS. 
