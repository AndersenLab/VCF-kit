#! /usr/bin/env python
"""
usage:
  vk primer template [options] <vcf>
  vk primer sanger   [options] <vcf>
  vk primer snip     [options] <vcf>
  vk primer indel    [options] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --ref=<reference>           Reference Genome
  --region=<region>           Restrict to region.
  --samples=<samples>         Output genotypes for a sample or set of samples. [default: ALL]
  --template=<template>       The sample template to output [default: ALT]
  --size=<int>                Amplicon size [default: 600-800]
  --box-variants              Add second column for the sequence with the variant boxed.
  --polymorphic               Only output variants that are polymorphic across specified samples.
  --enzymes=<enzymes>         snip-SNP only: Specify groups of restriction enzymes or individual enzymes [default: Common]
  --nprimers=<nprimers>       Maximum number of primers to generate [default: 5]

"""
from docopt import docopt
from utils import message
from utils.primer_vcf import primer_vcf
from utils.reference import *
from utils.fasta import *
import sys
from utils import check_program_exists

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]

def main(debug=None):
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)

    check_program_exists("primer3_core")
    check_program_exists("blastn")

    # Ensure user has specified a reference.
    if args["--ref"] is None:
        exit(message("Must specify a reference with --ref", color="red"))

    v = primer_vcf(args["<vcf>"], reference=args["--ref"], use_template=args['--template'], polymorphic=args["--polymorphic"])
    v.enzymes = args['--enzymes']
    v.nprimers = int(args['--nprimers'])
    # Region
    if args["--region"]:
        v.region = args["--region"]
    else:
        v.region = None

    # Samples
    if args["--samples"]:
        if args["--samples"] == "ALL":
            v.output_samples = v.samples
        else:
            v.output_samples = args["--samples"].split(",")
            for sample in v.output_samples:
                if sample not in v.samples + ["REF", "ALT", ]:
                    exit(message(sample + " not found in VCF", "red"))

    v.box_variants = args["--box-variants"] # Needs to be implemented.
    #v.amplicon_size = args["--size"]

    # Check for std. input
    if args["template"]:
        v.mode = "template"
        if '-' in args['--size']:
            try:
                v.region_size = int(args['--size'].split("-")[0])
                message("Warning: region size set to {s}".format(s=v.region_size))
            except:
                exit(message("Error: Invalid --size"))
        elif str.isdigit(args['--size']):
            v.region_size = int(args['--size'])
        v.amplicon_lower = 0
        v.amplicon_upper = 0
    elif args["indel"]:
        message("--size ignored; size is set dynamically when genotyping indels.")
        v.mode = "indel"
    elif args["snip"]:
        message("--size ignored; Set to 600-800 bp.")
        v.mode = "snip"
        v.amplicon_lower = 600
        v.amplicon_upper = 800
        v.region_size = 500 #x2=1000

    elif args["sanger"]:
        v.mode = "sanger"
        if args['--size'] is None:
            size = "600-800"
            message("Warning: --size set to 600-800 for sanger sequencing.")
        elif str.isdigit(args['--size']):
            exit(message("Error: You must specify a amplicon --size range for Sanger; e.g. 600-800."))
        else:
            size = args['--size']
        
        if "-" in size:
            try:
                v.amplicon_lower = int(size.split("-")[0])
                v.amplicon_upper = int(size.split("-")[1])
            except:
                exit(message("Error: Invalid --size"))
        else:
            v.amplicon_lower = 600
            v.amplicon_upper = 800
        v.region_size = (v.amplicon_upper//2) + 100
        if (v.amplicon_lower < 300 or 800 < v.amplicon_upper):
            message("Warning: --size range should (probably) be between 300-800 for sanger sequencing.")

    for variant in v.variant_iterator():
        variant.out()

if __name__ == '__main__':
    main()

