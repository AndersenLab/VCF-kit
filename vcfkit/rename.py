#! /usr/bin/env python
"""
usage:
  vk rename [--prefix=<prefix> --suffix=<suffix> --subst=<subst>...] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
import re
from signal import SIG_DFL, SIGPIPE, signal

from docopt import docopt
from vcfkit import __version__

from vcfkit.utils.vcf import *

signal(SIGPIPE, SIG_DFL)


def main(debug=None):
    args = docopt(__doc__,
                  argv=debug,
                  version=__version__)
    if args["<vcf>"] == "":
        print(__doc__)
    v = vcf(args["<vcf>"])
    for line in v.output_raw():
        if line.startswith("#CHROM"):
            line = line.split("\t")
            if args["--subst"]:
                find_replace = [re.split("[:=,]", x) for x in args["--subst"]]
                for orig, replacement in find_replace:
                    for n, sample in enumerate(line[9:]):
                        if sample == orig:
                            line[9+n] = replacement
            if args["--prefix"]:
                line[9:] = [args["--prefix"] + x for x in line[9:]]
            if args["--suffix"]:
                line[9:] = [x + args["--suffix"] for x in line[9:]]
            print('\t'.join(line))
        else:
            print((line.strip()))

if __name__ == '__main__':
    main()
