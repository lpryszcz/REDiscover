CHANGELOG
=========

v1.15
~~~~~
- both reads (read1 & read2) are processed for stranded libraries
- major allele is rescued for low frequency alternative alleles
- further optimisation for performance
  - pysam instead of `samtools mpileup` subprocess

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

TBD
~~~
- rescue major alt haplotype if the other is only very low freq
- add tools ie cosslinking with GTF and removing dbSNPs
