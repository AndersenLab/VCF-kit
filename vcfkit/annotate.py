#! /usr/bin/env python
"""
usage:
  tb annotate <vcf> 

"""
from docopt import docopt
from subprocess import Popen, PIPE
from signal import signal, SIGPIPE, SIG_DFL
from utils import autoconvert
from utils.vcf import *
import re
from utils.matrix import *
signal(SIGPIPE, SIG_DFL)

debug = None



ann_fields = ["allele", 
              "annotation",
              "putative_impact",
              "gene_name",
              "gene_id",
              "feature_type",
              "feature_id",
              "transcript_biotype",
              "rank_total",
              "hgvs_c",
              "hgvs_p",
              "cdna_position",
              "cds_position",
              "protein_position",
              "distance_to_feature",
              "errors"]

def parse_aa(hgvs_p):
    aa = re.split("\\.|[0-9]+", hgvs_p)[1:]
    aa = [three_to_one_aa[x] for x in aa]
    return tuple(aa)


def calc_grantham(aa_set):
    """
        Calculate Grantham Score
    """
    ret_set = []
    for i in aa_set:
        if i.count("*") == 0:
            ret_set.append(grantham[i])
        else:
            ret_set.append("")
    return "|".join(ret_set)


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug)
    if args["<vcf>"] == "":
        print(__doc__)
    v = vcf(args["<vcf>"])
    v.add_info_to_header({"ID": "grantham", "Description": "Grantham Score", "Type": "String", "Number": 1})
    for line in v:
        if "ANN" in dict(line.INFO):
            ANN = dict(line.INFO)["ANN"]
            ANN = [dict(zip(ann_fields, x.split("|"))) for x in ANN.split(",")]
            aa_set = [parse_aa(x["hgvs_p"]) for x in ANN]
            line.INFO["grantham"] = calc_grantham(aa_set)
            print(line)


