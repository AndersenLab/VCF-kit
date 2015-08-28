#! /usr/bin/env python
"""
usage:
  tb.py <command> [<args>...]
  tb.py -h | --help
  tb.py --version

commands:
  tajima

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
    print(args)
    argv = [args['<command>']] + args['<args>']
    if args["<command>"] == "":
      print(__doc__)
    elif args['<command>'] == 'tajima':
        # In case subcommand is a script in some other programming language:
        exit(call(['python', 'tajima.py'] + argv))
