#! /usr/bin/env python
"""
usage:
  tb.py <command> [<args>...]
  tb.py -h | --help
  tb.py --version

commands:
  tajima
  primer
  genome
  phylo
  freq
  geno
  vcf2tsv

"""
from docopt import docopt
from subprocess import call
from utils.vcf import *
import sys
import os


debug = None
if len(sys.argv) == 1:
    debug = [""]

def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

if __name__ == '__main__':
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=True)
    argv = [args['<command>']] + args['<args>']
    if args["<command>"] == "":
      print(__doc__)
    elif args['<command>'] in ['tajima', 'primer','genome','phylo','freq',"geno","vcf2tsv"]:
        comm = ['python', getScriptPath() + '/' + args["<command>"] + ".py"] + argv
        exit(call(comm))

