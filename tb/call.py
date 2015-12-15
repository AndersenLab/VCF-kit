#! /usr/bin/env python
"""
usage:
    tb call <seq.fasta> [<vcf>] [--ref=<reference> --all-sites]
    tb call alignments <seq.fasta>  [--ref=<reference>]

options:
    -h --help                   Show this screen.
    --version                   Show version.

"""
from docopt import docopt
import tb
from Bio import SeqIO
from utils.blastn import blast
from utils.vcf import *
from subprocess import Popen
import sys
from collections import defaultdict
from clint.textui import colored, puts, indent, puts_err
import os
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

from Bio.Blast.Applications import NcbiblastxCommandline


debug = None
if len(sys.argv) == 1:
    debug = ['tb','call', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)

    # Reference path - check that it exists
    module_path = os.path.split(os.path.realpath(__file__))[0]
    reference = module_path + "/reference/" + args["--ref"] + ".fa.gz"
    handle = open(args["<seq.fasta>"], "rU")

    # if args["alignments"]:
    #   for record in SeqIO.parse(handle, "fasta") :
    #     bl = blast_call(reference)
    #     for n, variant in enumerate(bl.blast(record.seq)):
    #         pass

    if args["<vcf>"]:
        concordance = True
        v = vcf(args["<vcf>"])
        samples = v.samples
    # Open fasta and read
    b = blast(reference)

    b.call_header(qual=True, sample=True)
    for record in SeqIO.parse(handle, "fastq"):
        sample, description = re.split("[ \|]{1}", record.name, 1)
        sample = sample.strip(">")
        blast_results = b.blast_call(record, all_sites = True)
        for n, variant in enumerate(blast_results):
            if variant is None:
                puts_err(colored.red("No Results for " + sample + " " + description))
                continue
            if n == 0 and args["<vcf>"]:
                vcf_variants = []
                for vcf_variant in v.fetch_variants(variant[0], variant[5], variant[6], samples = [sample]):
                    vcf_variants.append(vcf_variant[sample])
            #sanger_variant = [variant[i] for i in [0, 1, 3]]
            print '\t'.join(map(str, variant + [sample, description]))
