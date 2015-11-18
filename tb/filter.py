#! /usr/bin/env python
"""
usage:
  tb filter (REF|HET|ALT|MISSING) [--min=<min> --max=<max>] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from clint.textui import colored, puts, indent
from utils.vcf import *
from utils.fasta import *
from collections import defaultdict
import sys
import os
from glob import glob
from pprint import pprint as pp
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 


debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]


if __name__ == '__main__':
    #print debug
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug,
                  options_first=False)
    
    v = vcf(args["<vcf>"])
    n_samples = len(v.samples) * 1.0
    f = {}
    filter_s = [x for x in args.values() if x in ["REF","HET","ALT","MISSING"]][0]
    # Filter by rate or by number?
    print args
    if args["--min"]:
        if int(float(args["--min"])) != float(args["--min"]):
            filter_key_min = "r_" + filter_s
            filter_val_min = float(args["--min"])
        else:
            filter_key_min = filter_s
            filter_val_min = int(float(args["--min"]))
    if args["--max"]:
        if int(float(args["--max"])) != float(args["--max"]):
            filter_key_max = "r_" + filter_s
            filter_val_max = float(args["--max"])
        else:
            filter_key_max = filter_s
            filter_val_max = int(float(args["--max"]))
    for line in v:
        ok_to_print = True
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
                ok_to_print = False
        if args["--max"]:
            if f[filter_key_max] > filter_val_max:
                ok_to_print = False
        if ok_to_print:
            sys.stdout.write(str(line) + str(f) + "\t" + str(filter_val_max))