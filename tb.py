#!/usr/bin/env python
"""VCF-Toolbox 0.1

Usage:
  tb.py listvars <vcf>          
  tb.py plot <vcf> <x> [<y>]      [options]
  tb.py qc <vcf>                  [options]
  tb.py concordance <vcf> [--vcf2=<vcf2>] [--x=<x>] [--pairs=<pairset>]
  tb.py -h | --help
  tb.py --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
  --title=<title>             Set Custom plot titles.
  --region=<region>           Restrict analysis to a particular region.
  --include=<filter-expr>     Use a custom filtering string with bcftools.
  --facet=<facet-var>         Facet analysis on a categorical variable.
  --split-sample              When plotting genotype FORMAT fields, facet by sample.

"""
from docopt import docopt
from subprocess import call
from vcf import vcf
from utils import *
from plots import *


class opts:
  """ Defines options that can be overridden """
  functions = ""
  binwidth = ""
  add = ""


if __name__ == '__main__':
    args = docopt(__doc__, version='Naval Fate 2.0')
    print(args)

    v = vcf(args["<vcf>"])

    if args["listvars"] == True:
      v.list_vars()
    elif args["plot"] == True:
      if args["<y>"] is None:

        #======================#
        # Single Variable Plot #
        #======================#

        filename, r = v.query(args["<x>"])
        if args["<x>"] == "POS":
          # Facet by Chromosome Automatically
          print("")
          print(bc("Plotting Position; Automatically facetting by Chromosome","BOLD"))
          print("")
          # Setup Plot for chromosome.
          opts.add += " + \n facet_grid(.~CHROM, scales='free_x')"
          opts.add += " + \n scale_x_continuous(labels = genetic_scale) "
          opts.functions += genetic_scale
        if r["number"] == 1 and r["type"] in ["Integer","Float"]:
          print(bc("Creating Histogram of %s" % r["df"],"BOLD"))
          var1 = r["df"]
          Rcode = histogram.format(**locals())
        elif r["number"] == 1 and r["type"] in ["String"]:
          print(bc("Creating Bar Chart of %s" % r["df"],"BOLD"))
          var1 = r["df"]
          Rcode = barchart.format(**locals())
        
        with open(filename + ".R","w") as R:
          R.write(Rcode)
        call(["Rscript",filename + ".R"])

          
      else:

        #======================#
        # Two Variable Plot    #
        #======================#

        r = v.query(args["<x>"], args["<y>"])


    elif args["qc"] == True:
      print("List Variables")

    elif args["concordance"] == True:
      print v.compare_vcf(variable = args["--x"], pairs = args["--pairs"])

    else:
      pass
      
  # Run R script to generate plots
