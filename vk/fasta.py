#! /usr/bin/env python
"""
usage:
  vk tajima [--no-header --extra] <window-size> <step-size> <vcf> 
  vk tajima [--no-header --extra] <window-size> --sliding <vcf> 

options:
  -h --help                   Show this screen.
  --version                   Show version.  
  --window-size               blah
  --step-size                 blash 
  --extra                     display extra

command:
  tajima        Calculate Tajima's D


output:
    CHROM
    BIN_START
    BIN_END
    N_Sites
    N_SNPs
    TajimaD


"""
from docopt import docopt
from subprocess import call, Popen, PIPE
from itertools import combinations
from utils.vcf import *
from math import isinf
import sys
import os


debug = None



if __name__ == '__main__':
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug)
    if args["<vcf>"] == "":
      print(__doc__)
    wz = int(args["<window-size>"].replace(",",""))
    sz = None
    if not args["--sliding"]:
      sz = int(args["<step-size>"].replace(",",""))
    if args["--no-header"] == False:
        header_line = ["CHROM",
                       "BIN_START",
                       "BIN_END",
                       "N_Sites",
                       "N_SNPs",
                       "TajimaD"]
        if args["--extra"]:
            header_line += ["filename",
                            "window_size",
                            "step_size"]
        print "\t".join(header_line)
    for i in tajima(args["<vcf>"]).calc_tajima(wz, sz, extra = args["--extra"]):
        print(i)



