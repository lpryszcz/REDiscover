CHANGELOG
=========

v1.20
~~~~~
- mutual SNP info

v1.15
~~~~~
- auto-detection of stranded protocol (`-g/--gtf` needed)
  - supporting two types of stranded protocols `-fr-secondstrand` and `-fr-firststrand` (for more info: http://salmon.readthedocs.io/en/latest/library_type.html)
- both reads (read1 & read2) are processed for stranded libraries
  - this wasn't possible with samtools mpileup, so using pysam instead of `samtools mpileup` subprocess
- using `pysam` required further optimisation for performance (5-6x faster, as fast as `samtools mpileup`)
  - using `numpy.array`
  - adding read blocks instead of individual bases 
- stripped mean basecall quality from output
- major allele is rescued for low frequency alternative alleles
  - optionally report sites with multiple alternative alleles `-a / --report_alternatives`
- scripts
  - removing known SNP sites ie. dbSNP
  - plotting frequency histogram
  - annotating from GTF/GFF (TBD)
- more accurate
  - recognise and ignore duplicates (even those not annotated)
  - warn about samples with strandness
  - skip SNPs from antisense having coverage lower than sense strand (useful for poorly stranded samples)

v1.14
~~~~~
- stranded libraries support (only read1 is processed)
- optimised for performance

v1.13
~~~~~
- multithreading support

v1.12
~~~~~
- skip alt base calling if less than 3 reads (--minAltReads)
- -o / --output option added

v1.11
~~~~~
- use fasta file as reference
