#! /usr/bin/env python
"""
usage:
  tb.py phylo fasta <vcf>
  tb.py phylo tree (nj|upgma) [--plot] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from utils.vcf import *
from subprocess import Popen
import sys
from collections import defaultdict
from clint.textui import colored, puts, indent, puts_err

import os

debug = None
if len(sys.argv) == 1:
    debug = ['genome', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)
    module_path = os.path.split(os.path.realpath(__file__))[0]
    v = vcf(args["<vcf>"])
    samples = v.samples
    if args["fasta"] or args["tree"]:
        """
            Generate an aligned fasta from a VCF file.
        """
        seqs = {}
        for sample in samples:
          seqs[sample] = []
        for line in v:
          if line.is_snp:
            non_missing = [x.replace(".","-") for x in line.gt_bases]
            sample_gt = zip(samples,[x[-1] for x in non_missing])
            for sample, gt in sample_gt:
              seqs[sample].append(gt)
        if not args["tree"]:
          for sample, seq in seqs.items():
            print(">" + sample)
            print(''.join(seq))
        elif args["tree"]:
          """
            Generate a phylogenetic tree using an aligned fasta with muscle.
          """
          fasta = ""
          with indent(4):
              puts_err(colored.blue("\nGenerating Fasta\n"))
          for sample, seq in seqs.items():
              fasta += ">" + sample + "\n" + ''.join(seq) + "\n"
          tree_type = "upgma" # default is upgma
          if args["nj"]:
              tree_type = "neighborjoining"
          with indent(4):
              puts_err(colored.blue("\nGenerating " + tree_type + " Tree\n"))
          tree, err = Popen(["muscle","-maketree","-in","-","-cluster",tree_type], stdin = PIPE, stdout = PIPE).communicate(input = fasta)
          print(tree)
          if args["--plot"]:
              # R code for plotting here!
              pass



