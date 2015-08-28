#! /usr/bin/env python
"""
usage:
  tb.py tajima --option

Example

command:
  tajima        Calculate Tajima's D

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from subprocess import call
from utils.vcf import *
import sys



debug = None
if len(sys.argv) == 1:
    debug = ['tajima', '--option']



if __name__ == '__main__':
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=True)
    print(args)
    print "TAJIMA"

    x = vcf(args["<vcf>"])
