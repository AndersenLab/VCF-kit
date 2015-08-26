#!/usr/bin/env python
"""VCF-Toolbox 0.1

Usage:
  tb.py tajima <vcf>
  tb.py -h | --help
  tb.py --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
"""
from docopt import docopt
from subprocess import call
from utils.vcf import *
import sys



debug = None
if len(sys.argv) == 1:
    debug = ['tajima', '20150731_WI_PASS.vcf.gz']



if __name__ == '__main__':
    print(sys.argv)
    args = docopt(__doc__, version='VCF-Toolbox v0.1', argv = debug)
    print(args)
    x = vcf(args["<vcf>"])
    print(x.next())
    for i in x.window(windowsize=20, shift_type = "interval", shift_kind = "SNP"):
      print i 



