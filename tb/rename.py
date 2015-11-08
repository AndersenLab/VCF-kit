#! /usr/bin/env python
"""
usage:
  tb rename <vcf> [--prefix=<prefix> --suffix=<suffix> <orig:replace>...]

options:
  -h --help                   Show this screen.
  --version                   Show version.  


"""
from docopt import docopt
from utils.vcf import *
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

debug = None

if __name__ == '__main__':
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug)
    print args
    if args["<vcf>"] == "":
        print(__doc__)
    v = vcf(args["<vcf>"])
    for line in v.output_raw():
        if line.startswith("#CHROM"):
            line = line.split("\t")
            if args["--prefix"]:
                line[9:] = [args["--prefix"] + x for x in line[9:]]
            if args["--suffix"]:
                line[9:] = [x + args["--suffix"] for x in line[9:]]
            if args["<orig:replace>"]:
                find_replace = [x.split(":") for x in args["<orig:replace>"]]
                for orig, replacement in find_replace:
                    for n, sample in enumerate(line[9:]):
                        if sample == orig:
                            line[9+n] = replacement
            print '\t'.join(line)
        else:
            print line,