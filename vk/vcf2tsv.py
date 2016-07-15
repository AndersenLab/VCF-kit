#! /usr/bin/env python
"""
usage:
  vk vcf2tsv (wide|long) [--print-header --ANN] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from utils.vcf import *
from subprocess import Popen, PIPE
import sys
import os
import re
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

# Info
r_info = re.compile(r'''\#\#INFO=<
  ID=(?P<id>[^,]+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>Integer|Float|Flag|Character|String),
  Description="(?P<desc>[^"]*)".*
  >''', re.VERBOSE)

# Format
r_format = re.compile(r'''\#\#FORMAT=<
  ID=(?P<id>.+),
  Number=(?P<number>-?\d+|\.|[AG]),
  Type=(?P<type>.+),
  Description="(?P<desc>.*)".*
  >''', re.VERBOSE)

ANN_header = ["allele",
              "effect",
              "impact",
              "gene_name",
              "gene_id",
              "feature_type",
              "feature_id",
              "transcript_biotype",
              "exon_intron_rank",
              "nt_change",
              "aa_change",
              "cDNA_position/cDNA_len",
              "protein_position",
              "distance_to_feature",
              "error"]

debug = None
if len(sys.argv) == 1:
    debug = ['genome', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)
    module_path = os.path.split(os.path.realpath(__file__))[0]
    v = vcf(args["<vcf>"])
    field_set = []
    if len(sys.argv) == 1:
        print("Specify a VCF File")
        sys.exit()
    else:
        info = [m.groupdict()["id"] for m in r_info.finditer(v.raw_header)]
        format = [m.groupdict()["id"] for m in r_format.finditer(v.raw_header)]

    # Construct Query String
    print_header = ""
    if args["--print-header"]:
        print_header = "--print-header"
    if args["long"]:
        query_start = repr("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t" +
                           '\t'.join(['%' + x for x in info]) + "\t[-->%SAMPLE\t" +
                           '\t'.join(['%' + x for x in format]) + "\n]").strip("'")
    elif args["wide"]:
        query_start = repr("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t" +
                           '\t'.join(['%INFO/' + x for x in info]) +
                           "[\t%SAMPLE\t" + '\t'.join(['%' + x for x in format]) + "]\n").strip("'")
    comm = filter(len, ["bcftools", "query", print_header, "-f", query_start, v.filename])
    comm = Popen(comm, stdout=PIPE, stderr=PIPE)
    if args["--ANN"]:
        ANN_loc = info.index("ANN") + 7
    for n, line in enumerate(comm.stdout):
        if n == 0 and args["--print-header"] and args["wide"]:
            # Split out snpeff annotations
            if "ANN" in info and args["--ANN"]:
                line = line.split("\t")
                line = line[:ANN_loc - 1] + ANN_header + line[ANN_loc + 1:]
                line = '\t'.join(line)
            print re.sub("\[[0-9]+\]", "", line).strip("#\n ").replace(":", "_")
        elif n == 0 and args["long"] and args["--print-header"]:
            # Fix header for long format.
            line = re.sub("\[[0-9]+\]", "", line).strip("#\n ").split("\t")
            header = []
            for var in line:
                if ":" in var:
                    var = var.split(":")[1]
                header.append(var)
            if "ANN" in info and args["--ANN"]:
                header = header[:ANN_loc - 1] + ANN_header + header[ANN_loc + 1:]
            # Add Format to FORMAT
            header = header[:-len(format)] + ["F_" + x for x in header[-len(format):]]
            print("\t".join(header).replace("-->", ""))
        elif n < len(v.samples) and args["long"]:
            pass
        elif n >= len(v.samples) and args["long"]:
            line = line.strip("\n").split("\t")
            if line[0].startswith("-->"):
                line = fill_fields + line[len(line) - len(fill_fields):]
            else:
                fill_fields = line[0:(7 + len(info))]
            if "ANN" in info and args["--ANN"]:
                for var_effect in line[ANN_loc].split(","):
                    out_line = line[:ANN_loc - 1] + var_effect.split("|") + line[ANN_loc + 2:]
                    print('\t'.join(out_line).strip("\n").replace("-->", ""))
            else:
                print "\t".join(line).replace("-->", "")
        else:
            # Print wide format
            if "ANN" in info and args["--ANN"]:
                line = line.split("\t")
                for var_effect in line[ANN_loc].split(","):
                    out_line = line[:ANN_loc - 1] + var_effect.split("|") + line[ANN_loc + 2:]
                    print('\t'.join(out_line).strip("\n"))
            else:
                print(line.strip("\n"))
