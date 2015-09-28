#! /usr/bin/env python
"""
usage:
  tb.py primer [--ref=<reference>] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from clint.textui import colored, puts, indent
from utils.vcf import *
from utils.fasta import *
import sys
import os
from glob import glob

class af_vcf(vcf):
    """
        Subclass of vcf that calculates frequency of alleles by strain
    """
    def __init__(self, filename, reference):
        vcf.__init__(self,filename)


debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]


if __name__ == '__main__':
    #print debug
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    
    # Locate Reference
    v = af_vcf(args["<vcf>"], args["--ref"])
    