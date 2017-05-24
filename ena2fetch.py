#!/usr/bin/env python
# Fetches data from ENA

import os, sys, re

# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/001/ERR1442851/ERR1442851_1.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/002/ERR1442562/ERR1442562_1.fastq.gz
www = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/00%s/%s/%s_%s.fastq.gz"

fn = sys.argv[1]
for l in open(fn):
    ldata = l[:-1].split('\t')
    ldata = [re.sub('[: -]+', "_", x) for x in ldata]#; print ldata
    strain, tissue, stage, taxid, sample, experiment, project, run = ldata[1:9]
    stagedesc = stageid = ""
    if "_ZFS_" in stage:
        stagedesc, stageid = stage.split("_ZFS_")
        strain = tissue
    for i in [1, 2]:
        _www = www%(run[:6], run[-1], run, run, i)
        outfn = "%s.%s.%s.%s.%s_%s.fq.gz"%(strain, stageid, stagedesc, sample, run, i)
        cmd = "wget -nc -O %s %s"%(outfn, _www)
        print cmd

    