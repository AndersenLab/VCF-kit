#! /usr/bin/env python
"""
usage:
  tb.py vcf2tsv (wide|long) [--print-header] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from utils.vcf import *
from subprocess import Popen, PIPE
import sys
import os
import re
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

# Info
r_info = re.compile(r'''\#\#INFO=<
  ID=(?P<id>[^,]+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>Integer|Float|Flag|Character|String),
  Description="(?P<desc>[^"]*)".*
  >''', re.VERBOSE)

# Format
r_format = re.compile(r'''\#\#FORMAT=<
  ID=(?P<id>.+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>.+),
  Description="(?P<desc>.*)".*
  >''', re.VERBOSE)

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
    field_set = []
    if len(sys.argv) == 1:
        print("Specify a VCF File")
        sys.exit()
    else:
        info = [m.groupdict()["id"] for m in r_info.finditer(v.raw_header)]
        format = [m.groupdict()["id"] for m in r_format.finditer(v.raw_header)]

    # Construct Query String
    print_header = ""
    if args["--print-header"]:
        print_header = "--print-header"
    if args["long"]:
        query_start = repr("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t" + \
                      '\t'.join(['%' + x for x in info]) + "\t[%SAMPLE\t" + \
                      '\t'.join(['%' + x for x in format]) + "\n]").strip("'")
    elif args["wide"]:
        query_start = repr("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t" + \
                      '\t'.join(['%INFO/' + x for x in info]) + \
                      "[\t%SAMPLE\t" + '\t'.join(['%' + x for x in format]) + "]\n").strip("'")
    comm = filter(len,["bcftools", "query", print_header, "-f", query_start, v.filename])
    comm = comm[0:4] + ['"' + comm[4] + '"'] + comm[5:]
    #print ' '.join(comm)
    comm = Popen(comm, stdout = PIPE, stderr = PIPE)
    for n, line in enumerate(comm.stdout):
        if n == 0 and args["--print-header"] and args["wide"]:
            print re.sub("\[[0-9]+\]", "", line).strip("#\n ").replace(":","_")
        elif n == 0 and args["long"] and args["--print-header"]:
          # Fix header for long format.
          line = re.sub("\[[0-9]+\]", "", line).strip("#\n ").split("\t")
          header = []
          for var in line:
              if ":" in var:
                  var = var.split(":")[1]
              header.append(var)
          print("\t".join(header))
        elif n < len(v.samples) and args["long"]:
            pass
        else:
            print(line.strip("\n"))