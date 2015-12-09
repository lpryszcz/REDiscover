# REDiscover
Tool for **R**NA **e**diting **discover**y from NGS data. 

## Methodology

REDiscover uses internally (samtools mpileup)[http://www.htslib.org/] to genotype genome and transcriptome. 

By default, REDiscover filters out alignments failing quality criteria:
- secondary alignments
- QC failed reads
- duplicates
- reads with mapping quality (mapQ) below 15 

REDiscover reports only regions fulfilling several stringency criteria:
- depth-of-coverage
- mean basecall quality

Finally, reads with basecall quality below 20 (0.01 probability of error) for given positiong are ignored. 

## Usage

### Test dataset

## FAQ

## Citation
