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
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

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


    if args["<vcf>"]:
        concordance = True
        v = vcf(args["<vcf>"])
        samples = v.samples

    # Open fasta and read
    handle = open(args["<seq.fasta>"], "rU")
    print "chrom\tpos\tref\tgt\tvariant_type\tsanger_start\tsanger_end\tresult\tsample\tdescription"
    for record in SeqIO.parse(handle, "fasta") :
        bl = blast_call(reference)
        for n, variant in enumerate(bl.blast(record.seq)):
            if n == 0 and args["<vcf>"]:
                region = variant[0] + ":" +  str(variant[5]) + "-" + str(variant[6]) # Clean this up.
                comm = ["bcftools", "query","-f", """%CHROM\t%POS[\t%TGT\t%GT]\n""", "--regions", region, "--samples", record.name, args["<vcf>"]]
                ins, err = Popen(comm, stdout = PIPE, stderr = PIPE).communicate()
                vcf_variants = [x.split("\t") for x in ins.splitlines()]
                vcf_variants = [[chrom, int(pos), list(set(allele.split("/")))[0]] for chrom, pos, allele, gt in vcf_variants if gt == "1/1"]
                print vcf_variants
            sanger_variant = [variant[i] for i in [0,1,3]]
            if sanger_variant in vcf_variants:
                result = "TP"
            else:
                result = "FN"
            print '\t'.join(map(str, variant + [result, record.name, record.description]))



