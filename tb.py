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
from utils.vcf import *
import sys



debug = None
if len(sys.argv) == 1:
    debug = ['tajima', 'test.vcf.gz']



if __name__ == '__main__':
    print(sys.argv)
    args = docopt(__doc__, version='VCF-Toolbox v0.1', argv = debug)
    print(args)
    x = vcf(args["<vcf>"])

    print ["SNP", "sliding"]
    for i in x.window(windowsize=20, shift_method = ["SNP", "sliding"]):
      print i 

    x = vcf(args["<vcf>"])
    print ["POS", "sliding"]
    for i in x.window(windowsize=200, shift_method = ["POS", "sliding"]):
      print i 

    x = vcf(args["<vcf>"])
    print ["SNP", "interval"]
    for i in x.window(windowsize=20, shift_method = ["SNP", "interval"]):
      print i 



