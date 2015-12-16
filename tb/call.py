#! /usr/bin/env python
"""
usage:
    tb call <seq> --ref=<reference> [<vcf> --all-sites --vcf-targets]
    tb call alignments <seq>  [--ref=<reference>]

options:
    -h --help                   Show this screen.
    --version                   Show version.

"""
from docopt import docopt
import tb
from Bio import SeqIO
from utils.blastn import blast, blast_variant
from utils.vcf import *
from subprocess import Popen
import sys
from collections import defaultdict
from clint.textui import colored, puts, indent, puts_err
import os
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

from Bio.Blast.Applications import NcbiblastxCommandline


def seq_type(filename):
    filename, ext = os.path.splitext(filename.lower())
    if ext in [".fasta", ".fa"]:
        return 'fasta'
    elif ext in [".fastq",".fq"]:
        return 'fastq'
    elif ext in [".ab1"]:
        return 'ab1'
    else:
        raise Exception("Unknown sequence file type: " + filename)

def format_gt(gt):
    # Return homozygous calls as single bases.
    gt = list(set(re.split("[\|\/]", gt)))
    if len(gt) == 1:
        return gt[0]
    else:
        return '/'.join(gt)


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
    handle = open(args["<seq>"], "rU")

    if args["<vcf>"]:
        concordance = True
        v = vcf(args["<vcf>"])
        samples = v.samples
    # Open fasta and read
    b = blast(reference)

    # Set file type:
    sequence_file_type = seq_type(args["<seq>"])

    # Output header
    print("\t".join(blast_variant.output_order + ["classification", "sample", "description"]))

    for record in SeqIO.parse(handle, sequence_file_type):
        rec_split = re.split("[ \|]{1}", record.name, 1)
        sample = record.name.strip(">")
        description = record.description.strip(">")
        blast_results = b.blast_call(record, all_sites = args["--all-sites"])
        for n, variant in enumerate(blast_results):
            if variant is None:
                puts_err(colored.red("No Results for " + sample + " " + description))
                continue
            if args["<vcf>"]:
                if n == 0:
                    vcf_variants = []
                    for vcf_variant in v(variant.region()):
                        # Deal with phasing
                        gt = format_gt(vcf_variant.gt_bases[v.samples.index(sample)])
                        vcf_variants.append([vcf_variant.CHROM,
                                             vcf_variant.POS,
                                             gt])
                        vcf_variant_positions = [x[0:2] for x in vcf_variants]

                chrom_pos =  variant.chrom_pos_allele()[0:2]
                vcf_variant_match = [x for x in vcf_variants if x[0:2] == chrom_pos]
                if vcf_variant_match:
                    vcf_variant_match = vcf_variant_match[0]
                    variant.vcf_gt = vcf_variant_match[2]
                    if variant.REF == vcf_variant_match[2]:
                        classification = "TN"
                    elif variant.gt == vcf_variant_match[2]:
                        classification = "TP" 
                    elif variant.gt != vcf_variant_match[2]:
                        classification = "TP"
                else:
                    classification = "FN"

                if (args["--vcf-targets"] and variant.chrom_pos_allele()[0:2] in vcf_variant_positions) or variant.is_variant:
                    print '\t'.join([str(variant), classification, sample, description])
                elif args["--vcf-targets"] is False:
                    print '\t'.join([str(variant), classification, sample, description])
            else:
                print '\t'.join([str(variant), classification, sample, description])
