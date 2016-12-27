#! /usr/bin/env python
"""
usage:
    vk call <seq> --ref=<reference> (--all-sites|--vcf-sites) <vcf>

options:
    -h --help                   Show this screen.
    --version                   Show version.

"""
from docopt import docopt
import vk
from utils.reference import resolve_reference_genome
from Bio import SeqIO
from utils.blastn import blast, blast_variant
from utils.vcf import *
from subprocess import Popen
import sys
from collections import defaultdict
from clint.textui import colored, puts, puts_err, indent
import os
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

from Bio.Blast.Applications import NcbiblastxCommandline


def seq_type(filename):
    """
        Resolves sequence filetype using extension.
    """
    filename, ext = os.path.splitext(filename.lower())
    if ext in [".fasta", ".fa"]:
        extension = 'fasta'
    elif ext in [".fastq",".fq"]:
        extension = 'fastq'
    elif ext in [".ab1", '.abi']:
        extension = 'abi'
    else:
        raise Exception("Unknown sequence file type: " + filename)

    with indent(4):
        puts_err(colored.green("\nReading sequences as %s\n" % extension.upper()))
    return extension


def resolve_sample_from_line(samples, line):
    """
        Resolves sample names by splitting fasta line
        on non-word characters.
    """
    line = re.split("\W", line)
    matched_sample = [x for x in samples if x in line]
    if len(matched_sample) == 1:
        return matched_sample[0]
    return ""


def format_gt(gt):
    # Return homozygous calls as single bases.
    gt = list(set(re.split("[\|\/]", gt)))
    if len(gt) == 1:
        return gt[0]
    else:
        return '/'.join(gt)

def format_args(args, add_missing_stdin = False):
    if add_missing_stdin:
        pass # Test for vcf

def main(debug=None):
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)


    module_path = os.path.split(os.path.realpath(__file__))[0]
    handle = open(args["<seq>"], "rb")
    reference = resolve_reference_genome(args["--ref"])

    if args["<vcf>"]:
        concordance = True
        v = vcf(args["<vcf>"])
        samples = v.samples


    if args["--vcf-sites"] and args["<vcf>"] is None:
        with indent(4):
            exit(puts_err(colored.red("\nMust specify <vcf> with --vcf-sites\n")))

    # Setup reference for blast call
    b = blast(reference)

    # Set file type:
    sequence_file_type = seq_type(args["<seq>"])

    # Output header
    print("\t".join(blast_variant.output_order))
    for record in SeqIO.parse(handle, sequence_file_type):
        # Resolve sample within fasta line
        sample = resolve_sample_from_line(samples, handle.name)
        if not sample:
            sample = resolve_sample_from_line(samples, record.name)
        blast_results = b.blast_call(record)
        classification = ""
        for n, variant in enumerate(blast_results):
            output_line = False
            if variant is None:
                puts_err(colored.red("No Results for " + sample + " " + record.description))
                continue
            if args["<vcf>"]:
                if n == 0:
                    vcf_variants = []
                    for vcf_variant in v(variant.region()):
                        if sample:
                            gt = format_gt(vcf_variant.gt_bases[v.samples.index(sample)])
                            vcf_variants.append([vcf_variant.CHROM,
                                                 vcf_variant.POS,
                                                 gt,
                                                 vcf_variant.REF,
                                                 vcf_variant.ALT])
                            vcf_variant_positions = [x[0:2] for x in vcf_variants]

                chrom_pos =  variant.chrom_pos_allele()[0:2]
                vcf_variant_match = [x for x in vcf_variants if x[0:2] == chrom_pos]
                if vcf_variant_match:
                    vcf_variant_match = vcf_variant_match[0]
                    variant.vcf_gt = vcf_variant_match[2]
                    variant.REF = vcf_variant_match[3]
                    variant.ALT = ','.join(vcf_variant_match[4])
                    variant.fetch_variant_type()
                    if variant.REF == variant.seq_gt and variant.seq_gt == variant.vcf_gt:
                        variant.classification = "TN"
                    elif variant.REF != variant.seq_gt and variant.seq_gt == variant.vcf_gt:
                        variant.classification = "TP"
                    elif variant.REF == variant.seq_gt and variant.seq_gt != variant.vcf_gt:
                        variant.classification = "FP" 
                    elif variant.REF != variant.seq_gt and variant.seq_gt != variant.vcf_gt:
                        variant.classification = "FN"
                else:
                    variant.REF = ""
                    variant.ALT = ""
                    variant.fetch_variant_type()
                    variant.classification = ""

                if args["--vcf-sites"] and variant.classification != "":
                    output_line = True
                elif args["--all-sites"] is True:
                    output_line = True
            else:
                if args["--all-sites"]:
                    output_line = True
                elif variant.is_variant:
                    output_line = True
            if output_line:

                variant.sample = sample
                if record.description:
                    variant.description = record.description
                else:
                    variant.description = os.path.split(handle.name)[1]
                print '\t'.join([str(variant)])

if __name__ == '__main__':
    main()
