# Test directory
This directory hosts test datasets.
To run the test example, execute: 
```bash
./REDiscover -d test/DNA.bam -r test/RNA.bam


./REDiscover -v -r test/PRJNA345214/*.bam -f test/PRJNA345214/dna.fa -o test/PRJNA345214 -t 4 -fr-secondstrand


s=PRJNA345214; #PRJNA266390
./REDiscover -v -r test/$s/*.bam -f test/$s/dna.fa -o test/$s.gz -t 4 -fr-secondstrand
./fasta2intersection.py -r test/$s/rna.fa -i test/$s.gz.pos2mi.gz -o test/$s.gz.pos2mi.gz.mods

```


**Table of contents**  
- [Generation of test set](#Generation-of-test-set)
  - [REDiscover pipeline](#rediscover-pipeline)
- [Accuracy estimation](#accuracy-estimation)


## Generation of test set


### REDiscover pipeline


## Accuracy estimation
