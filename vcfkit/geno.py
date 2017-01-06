#! /usr/bin/env python
"""
usage:
  vk geno transfer-filter <vcf>
  vk geno het-polarization <vcf>


tb.py geno transfer-filter

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from vk import __version__
from docopt import docopt
from utils.vcf import *
from utils.fasta import *
import sys
import math
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)
np.set_printoptions(threshold=np.nan)


def phred2p(phred):
    return 10**(phred / -10.0)


def GL2PL(gl):
    """ Converts Genotype likelyhoods to phred scaled (PL) genotype likelyhoods. """
    return -int(gl * 10)

debug = None
if len(sys.argv) == 1:
    debug = ['primer', "--ref=WBcel235", "test.vcf.gz"]


if __name__ == '__main__':
    # print debug
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)
    # Locate Reference
    v = vcf(args["<vcf>"])
    format_added = False
    if args["transfer-filter"]:
        for line in v.output_raw():
            line = line.strip()
            if line.startswith("#CHROM"):
                # Get Sample information and count
                samples = line.strip().split("\t")[9:]
            elif line.startswith("#"):
                # Add Info line for het polarization flag
                if line.startswith("##FORMAT") and format_added is False:
                    format_added = True
                    line = line + "\n##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype-level filter\">"
            else:
                line = line.split("\t")
                FILTER = line[6]
                line[8] = line[8] + ":" + "FT"
                line[9:] = [x + ":" + FILTER for x in line[9:]]
                line = '\t'.join(line)
            sys.stdout.write(line + "\n")
    elif args["het-polarization"]:
        for line in v.output_raw():
            line = line.strip()
            if line.startswith("#CHROM"):
                # Get Sample information and count
                line = "##VCF-kit(version {v}) command=geno het-polarization\n".format(v = __version__) + line
                samples = line.split("\t")[9:]
            elif line.startswith("##FORMAT=<ID=GL,"):
                line = """##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">"""
            elif line.startswith("#"):
                # Add Info line for het polarization flag
                if line.startswith("##FORMAT") and format_added is False:
                    format_added = True
                    line += "\n##FORMAT=<ID=HP,Number=1,Type=String,Description=\"Flag used to mark whether a variant was polarized\">"
                # Pass comment lines.
            else:
                if line.split('\t')[8].find("GL") > 0:
                    # Replace GL with PL scores.
                    line = line.split('\t')
                    GL_loc = line[8].split(":").index("GL")
                    line[8] = line[8].replace("GL", "PL")
                    geno_set = []
                    for k, v in enumerate(line[9:]):
                        FT = v.split(":")
                        GL_set = FT[GL_loc].split(",")
                        try:
                            FT[GL_loc] = ','.join([str(GL2PL(float(i))) for i in GL_set])
                        except:
                            FT[GL_loc] = ",".join(GL_set)
                        geno_set.append(':'.join(FT))
                        line = line[0:9] + geno_set
                    line = '\t'.join(l)
                if line.find("0/1") > 0:
                    line = line.strip().split("\t")
                    if line[8].find("PL") > -1 and len(line[4].split(",")) == 1 and line[8].find("HP") == -1:
                        PL = line[8].split(":").index("PL")
                        add_HP_flag = 0
                        for k, v in enumerate(line[9:]):
                            PL_set = v.split(":")[PL].split(",")
                            v = v.split(":")
                            if len(PL_set) == 3:
                                PL_set = [phred2p(int(i)) for i in PL_set]
                                log_score = -math.log10(PL_set[0] / PL_set[2])
                                if add_HP_flag == 0:
                                    if line[8].find("HP") == -1:
                                        line[8] = line[8] + ":HP"
                                    add_HP_flag = 1
                                if (log_score < -2):
                                    v[0] = "0/0"
                                    if v[-1] not in ["AA", "AB", "BB"]:
                                        line[k + 9] = v + ["AA"]
                                    else:
                                        line[k + 9] = v[0:-1] + ["AA"]
                                elif (log_score > 2):
                                    v[0] = "1/1"
                                    if v[-1] not in ["AA", "AB", "BB"]:
                                        line[k + 9] = v + ["BB"]
                                    else:
                                        line[k + 9] = v[0:-1] + ["BB"]
                                else:
                                    v[0] = "0/1"
                                    if v[-1] not in ["AA", "AB", "BB"]:
                                        line[k + 9] = v + ["AB"]
                                    else:
                                        line[k + 9] = v[0:-1] + ["AB"]
                            else:
                                if add_HP_flag == 0:
                                    line[k + 9] = v
                                else:
                                    line[k + 9] = v + ["."]
                            line[k + 9] = ":".join(line[k + 9])
                    line = "\t".join(line)
            print(line)
