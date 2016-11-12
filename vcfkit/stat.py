#! /usr/bin/env python
"""
usage:
  vk stat <vcf>

"""
from docopt import docopt
from subprocess import Popen, PIPE
from signal import signal, SIGPIPE, SIG_DFL
from utils import autoconvert
import re
signal(SIGPIPE, SIG_DFL)


debug = None

class out_line:
    def __init__(self):
        self.vcf = ""
        self.e1 = ""
        self.e2 = ""
        self.e3 = ""
        self.val = ""

    def __str__(self):
        return '\t'.join([self.vcf, self.e1, self.e2, self.e3, self.val])


class stat:
    """
      Class for parsing and manipulating bcf
    """


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug)
    if args["<vcf>"] == "":
        print(__doc__)
    out = out_line()
    out.vcf = args["<vcf>"]
    comm = ["bcftools", "stats", "--depth", "0,1000,1", "-s", "-", "-v", args["<vcf>"]]
    result = Popen(comm, stdout = PIPE, stderr = PIPE)
    for n, line in enumerate(result.stdout):
        line = re.sub("\[[0-9]+\]","", line)
        line = line.strip().split("\t")
        if len(line) > 1:
            if line[0].startswith("# SN"):
                out.e1 = "Summary numbers"
            elif line[0].startswith("SN"):
                out.e2 = line[2].strip(":")
                out.val = line[3]
                print(out)
            elif line[0] in ["# TSTV", "# SiS", "# ST"]:
                variables = line[2:]
                out.e2 = line[0].strip("# ")
            elif line[0] in ["# AF", "# DP",  "# HWE", "# PSC"]:
                variables = line[2:]
                out.e2 = line[0].strip("# ")
            else:
                out.e1 = line[0]
                if line[0] in ["TSTV", "SiS"]:
                    out.e3 = ""
                    for k, v in zip(variables, line[2:]):
                        out.e2 = k
                        out.val = v
                        print(out)
                elif line[0] in ["ST"]:
                    out.e2 = line[2]
                    out.e3 = ""
                    out.val = line[3]
                    print(out)
                elif line[0] in ["AF", "QUAL", "DP", "HWE", "PSC"]:
                    out.e1 = line[0]

                    out.e3 = line[2]
                    for k,v in zip(variables, line[2:]):
                        out.e2 = k
                        if type(autoconvert(v)) != str:
                            out.val = v
                            print(out)

