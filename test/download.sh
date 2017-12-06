#!/bin/bash
# Download REDiscover data sets

# B. subtilis
www=zdglab.iimcb.gov.pl/cluster/rna_modifications/subtilis/last/PRJNA345214
wget -nc -q -r $www; ln -s $www

# E. coli
www=zdglab.iimcb.gov.pl/cluster/rna_modifications/ecoli/last/PRJNA266390
#wget -nc -q -r $www; ln -s $www
