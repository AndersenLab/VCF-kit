#! /usr/bin/env python
"""
usage:
  vk filter (REF|HET|ALT|MISSING) [--min=<min> --max=<max> --soft-filter=<soft> --mode=(+|x)] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from clint.textui import colored, puts, indent
from utils.vcf import *
from utils.fasta import *
from utils import message
from collections import defaultdict
import sys
import os
from glob import glob
from pprint import pprint as pp
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 


def main(debug = None):
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    if args["--soft-filter"] and not args["--mode"]:
        exit(message("Must Specify --mode with soft-filter"))
    v = vcf(args["<vcf>"])
    n_samples = len(v.samples) * 1.0
    f = {}
    filter_s = [x for x in args.values() if x in ["REF","HET","ALT","MISSING"]][0]
    # Filter by rate or by number?
    if args["--min"]:
        direction = "<"
        if int(float(args["--min"])) != float(args["--min"]):
            filter_key_min = "r_" + filter_s
            filter_val_min = float(args["--min"])
            filter_type = "FREQUENCY"
        else:
            filter_key_min = filter_s
            filter_val_min = int(float(args["--min"]))
            filter_type = "COUNT"
        filter_value = filter_val_min
    if args["--max"]:
        direction = ">"
        if int(float(args["--max"])) != float(args["--max"]):
            filter_key_max = "r_" + filter_s
            filter_val_max = float(args["--max"])
            filter_type = "FREQUENCY"
        else:
            filter_key_max = filter_s
            filter_val_max = int(float(args["--max"]))
            filter_type = "COUNT"
        filter_value = filter_val_max

    # Output header
    header = v.raw_header.splitlines()
    for n, i in enumerate(header):
        if i.startswith("##FILTER") and args["--soft-filter"]:
            filter_name = args["--soft-filter"]
            filter_line = """##FILTER=<ID={filter_name},Description="Apply filter if {filter_type}({filter_s}) {direction} {filter_value}">""".format(**locals())
            header.insert(n+1, filter_line)
            break
    header = '\n'.join(header) + "\n"
    sys.stdout.write(header)
    for line in v:
        filtered = False
        f["ALT"] = line.num_hom_alt
        f["HET"] = line.num_het
        f["REF"] = line.num_hom_ref
        f["MISSING"] = int(n_samples - line.num_called)
        f["r_ALT"] = f["ALT"] / n_samples
        f["r_HET"] = f["HET"] / n_samples
        f["r_REF"] = f["REF"] / n_samples
        f["r_MISSING"] = f["MISSING"] / n_samples
        if args["--min"]:
            if f[filter_key_min] < filter_val_min:
                filtered = True
        if args["--max"]:
            if f[filter_key_max] > filter_val_max:
                filtered = True
        if args["--soft-filter"]:
            line = str(line).split("\t")
            if args["--mode"] == "x":
                line[6] = "PASS"
            if filtered is False:
                sys.stdout.write('\t'.join(line))
            else:
                if args["--mode"] == "+":
                    if line[6] == "PASS":
                        line[6] = ""
                    line[6] = ';'.join([line[6]] + [args["--soft-filter"]]).strip(";")
                elif args["--mode"] == "x":
                    line[6] = args["--soft-filter"]
                sys.stdout.write('\t'.join(line))
        elif filtered is False:
            sys.stdout.write(str(line))

if __name__ == '__main__':
    main()