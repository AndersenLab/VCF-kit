#! /usr/bin/env python
"""
usage:
  tb.py <command> [<args>...]
  tb.py -h | --help
  tb.py --version

commands:
  tajima
  primer

"""
from docopt import docopt
from subprocess import call
from utils.vcf import *
import sys



debug = None
if len(sys.argv) == 1:
    debug = [""]


if __name__ == '__main__':
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=True)
    argv = [args['<command>']] + args['<args>']
    if args["<command>"] == "":
      print(__doc__)
    elif args['<command>'] == 'tajima':
        exit(call(['python', 'tajima.py'] + argv))
    elif args['<command>'] == 'primer':
        exit(call(['python', 'tajima.py'] + argv))
