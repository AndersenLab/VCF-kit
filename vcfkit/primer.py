#! /usr/bin/env python
"""
usage:
  vk primer template [options] <vcf>
  vk primer sanger   [options] <vcf>
  vk primer snip     [options] <vcf>
  vk primer indel    [options] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --ref=<reference>           Reference Genome
  --region=<region>           Restrict to region.
  --samples=<samples>         Output genotypes for a sample or set of samples. [default: ALL]
  --template=<sample>         Sequence to use for template. Use REF, ALT, or sample name. [default: ALT]
  --mark-variant              Use brackets to indicate the location of the variant in output.
  --size=<int>                Amplicon size/2 (Upstream and downstream length) [default: 50]
  --box-variants              Add second column for the sequence with the variant boxed.


"""
from docopt import docopt
from utils import message
from utils.primer_vcf import primer_vcf
from utils.reference import *
from utils.fasta import *
import sys
from utils import check_program_exists

import signal
signal.signal(signal.SIGINT, lambda x,y: sys.exit(0))

debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]


if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)

    check_program_exists("primer3_core")

    # Ensure user has specified a reference.
    if args["--ref"] is None:
        exit(message("Must specify a reference with --ref", color="red"))

    v = primer_vcf(args["<vcf>"], reference=args["--ref"], use_template=args["--template"])

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

    v.box_variants = args["--box-variants"]
    v.region_size = int(args["--size"])

    # Check for std. input
    if args["template"]:
        v.mode = "template"

    elif args["indel"]:
        if args["--size"]:
            message("Warning: --size ignored; size is set dynamically when genotyping indels.")
        v.mode = "indel"

    elif args["snip"]:
        if args["--size"]:
            message("Warning: --size ignored; size is set to ~1000 bp templates.")
        v.mode = "snip"
        v.region_size = 500

    for variant in v.variant_iterator():
        print(variant)


