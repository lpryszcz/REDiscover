#!/usr/bin/env python
desc="""Return intersection between modified bases and SNPs
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 8/12/2017
"""

import gzip, os, sys
import numpy as np
from collections import Counter
from datetime import datetime

from modifications2rna import fasta_parser, table2modifications
from modifications2signatures import load_modifications

def fasta2intersection(out, infname, rna, table, verbose):
    """Get intersection between modified bases and SNPs"""
    mod2base, mod2name = table2modifications(table)

    # load modifications
    mods, unmods = load_modifications(rna)

    pos2mod = {p: m for m in mods for p in mods[m]}
    
    for l in gzip.open(infname): 
        ldata = l[:-1].split('\t')
        if l.startswith('##') or len(ldata)<3:
            continue
        ref, pos, snp = ldata[:3]
        coord = "%s:%s"%(ref, pos)
        name = "-"
        if coord in pos2mod:
            mod = pos2mod[coord]
            name = mod2name[mod]
        out.write("%s\t%s"%(name, l))
    

def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.1') 
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), help="output stream [stdout]")
    parser.add_argument("-i", "--input", required=1, help="REDiscover output file to process")
    #parser.add_argument("-d", "--dna", required=1, help="DNA FastA")
    parser.add_argument("-r", "--rna", required=1, help="RNA FastA")
    parser.add_argument("-t", "--table", default="test/modifications.txt", help="modification table [%(default)s]" )
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    fasta2intersection(o.output, o.input, o.rna, o.table, o.verbose)

if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    