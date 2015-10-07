#! /usr/bin/env python
"""
usage:
  tb.py phylo fasta <vcf>
  tb.py phylo nexus <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from utils.vcf import *
import sys
from collections import defaultdict
from clint.textui import colored, puts, indent
from nexus import NexusWriter

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
    if args["fasta"]:
        seqs = {}
        for sample in samples:
          seqs[sample] = []
        for line in v:
          if line.is_snp:
            non_missing = [x.replace(".","-") for x in line.gt_bases]
            sample_gt = zip(samples,[x[-1] for x in non_missing])
            for sample, gt in sample_gt:
              seqs[sample].append(gt)
        for sample, seq in seqs.items():
            print(">" + sample)
            print(''.join(seq))
    elif args["nexus"]:
        n = NexusWriter()
        for line in v:
            if line.is_snp:
                non_missing = [x.replace(".","?") for x in line.gt_bases]
                sample_gt = zip(samples,[x[-1] for x in non_missing])
                for sample, gt in sample_gt:
                    n.add(sample, "chrom", gt)
        print n.make_nexus(interleave=True, charblock=True) 
