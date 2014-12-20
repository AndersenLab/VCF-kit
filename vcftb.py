#!/usr/bin/env python
"""VCF-Toolbox.

Usage:
  vcftb.py listvars <vcf>          
  vcftb.py plot <vcf> <x> [<y>]      [options]
  vcftb.py QC <vcf> <x> [<y>]        [options]
  vcftb.py compare <vcflist>...          [options]
  vcftb.py report <vcf>              [options]
  vcftb.py -h | --help
  vcftb.py --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
  --title=<title>             Set Custom plot titles.
  --region=<region>           Restrict analysis to a particular region.
  --include=<filter-expr>     Use a custom filtering string with bcftools.
  --facet=<facet-var>         Facet analysis on a categorical variable.
  --split-format              When plotting genotype FORMAT fields, facet by sample.

"""
from docopt import docopt
from vcf import vcf
from utils import *

if __name__ == '__main__':
    args = docopt(__doc__, version='Naval Fate 2.0')
    #print(args)

    v = vcf(args["<vcf>"])

    if args["listvars"] == True:
      v.list_vars()
    elif args["plot"] == True:
      if args["<y>"] is None:
        # Single Variable Plot
        v.query(args["<x>"])
      else:
        # Two Variable Plot
        print "Two Variable Plot"

    elif args["QC"] == True:
      print("List Variables")

    elif args["compare"] == True:
      print("List Variables")

    elif args["report"] == True:
      print("Generate summary report")
      