#!/usr/bin/env python
"""VCF-Toolbox 0.1

Usage:
  tb.py listvars <vcf>          
  tb.py plot <vcf> <x> [<y>]      [options]               
  tb.py concordance <vcf> [--vcf2=<vcf2>] [--x=<x>] [--pairs=<pairset>]
  tb.py qc <vcf>   
  tb.py tstv <vcf> [--x=<x>]
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


class opts:
  """ Defines options that can be overridden """
  title = ""
  functions = ""
  binwidth = ""
  add = ""

  # Log
  log_x_open = ""
  log_x_close = ""


if __name__ == '__main__':
    args = docopt(__doc__, version='Naval Fate 2.0')
    print(args)

    v = vcf(args["<vcf>"])
    #===============#
    # Parse Options #
    #===============#

    if args["--title"] is not None:
      opts.title = args["--title"]

    # Log Transformations
    if args["<x>"] is not None:
      if args["<x>"].startswith("log:"):
        opts.log_x_open = "log10("
        opts.log_x_close = ")"
        args["<x>"] = args["<y>"].replace("log:", "")
    if args["<y>"] is not None:
      if args["<y>"].startswith("log:"):
        opts.log_y_open = "log10("
        opts.log_y_close = ")"
        args["<y>"] = args["<y>"].replace("log:", "")

    # Check that specified regions exist
    if args["--region"] is not None:
      region_check = [x.split(":") for x in args["--region"].split(",")]
      for chrom,bp_range in region_check:
        if chrom not in v.contigs.keys():
          error("The CHROM you specified (%s) does not exist." % chrom)
        within_lower = within_range(int(bp_range.split("-")[0]), 0, v.contigs[chrom]["length"])
        within_upper = within_range(int(bp_range.split("-")[1]), 0, v.contigs[chrom]["length"])
        if not within_lower or not within_upper:
          error("The range (%s) falls outside of CHROM length (%s)" % (bp_range, v.contigs[chrom]["length"]))
    else:
      args["--region"] = ""

    if args["--include"] is None:
      args["--include"] = ""

    if args["listvars"] == True:
      v.list_vars()
    elif args["plot"] == True:
      if args["<y>"] is None:

        #======================#
        # Single Variable Plot #
        #======================#

        query, analysis_dir, filename, r = v.query(args["<x>"], region = args["--region"])
        os.chdir(analysis_dir)

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
          histogram = get_plot("histogram")
          Rcode = histogram.format(**locals())
        elif r["number"] == 1 and r["type"] in ["String"]:
          print(bc("Creating Bar Chart of %s" % r["df"],"BOLD"))
          var1 = r["df"]
          barchart = get_plot("barchart")
          Rcode = barchart.format(**locals())
        
        with open(filename + ".R","w") as R:
          R.write(Rcode)
        call(["Rscript",filename + ".R"])

          
      else:

        #======================#
        # Two Variable Plot    #
        #======================#

        r = v.query(args["<x>"], args["<y>"])


    elif args["tstv"] == True:
      """ Produces a tstv report """
      v.tstv(args["--x"])

    elif args["concordance"] == True:
      print v.compare_vcf(variable = args["--x"], pairs = args["--pairs"])

    else:
      pass
      
  # Run R script to generate plots
