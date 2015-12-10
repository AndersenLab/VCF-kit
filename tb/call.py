#! /usr/bin/env python
"""
usage:
    tb call <seq.fasta> [<vcf>] [--ref=<reference>]

options:
    -h --help                   Show this screen.
    --version                   Show version.

"""
from docopt import docopt
import tb
from Bio import SeqIO
from utils.blastn import blast_call
from utils.vcf import *
from subprocess import Popen
import sys
from collections import defaultdict
from clint.textui import colored, puts, indent, puts_err
import os

from Bio.Blast.Applications import NcbiblastxCommandline


debug = None
if len(sys.argv) == 1:
        debug = ['genome', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                version='VCF-Toolbox v0.1',
                argv=debug,
                options_first=False) 

    # Reference path - check that it exists
    module_path = os.path.split(os.path.realpath(__file__))[0]
    reference = module_path + "/reference/" + args["--ref"] + ".fa.gz"

    # Open fasta and read
    handle = open(args["<seq.fasta>"], "rU")
    print "chrom\tpos\tref\talt\tsample\tdescription"
    for record in SeqIO.parse(handle, "fasta") :
        bl = blast_call(reference)
        for variant in bl.blast(record.seq):
            print '\t'.join(map(str, variant + [record.name, record.description]))

    if args["<vcf>"]:
        concordance = True
        v = vcf(args["<vcf>"])
        samples = v.samples



