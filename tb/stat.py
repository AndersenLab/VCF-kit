#! /usr/bin/env python
"""
usage:
  tb stat <vcf> 

"""
from docopt import docopt
from subprocess import Popen, PIPE
from signal import signal, SIGPIPE, SIG_DFL
from utils import autoconvert
import re
signal(SIGPIPE, SIG_DFL)


debug = None

class out_line:

    def __str__(self):
        return '\t'.join([self.vcf, self.e1, self.e2, self.val])


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
    comm = ["bcftools", "stats", "-s", "-", "-v", args["<vcf>"]]
    result = Popen(comm, stdout = PIPE, stderr = PIPE)
    for n, line in enumerate(result.stdout):
        line = re.sub("\[[0-9]+\]","", line)
        line = line.strip().split("\t")
        if n < 100:
            if len(line) > 1:
                if line[0].startswith("# SN"):
                    out.e1 = "Summary numbers"
                elif line[0].startswith("SN"):
                    out.e2 = line[2].strip(":")
                    out.val = line[3]
                    print(out)
                elif line[0] in ["# TSTV", "# SiS"]:
                    variables = line[2:]
                    section = line[0]
                else:
                    if line[0] in ["TSTV", "SiS"]:
                        for k, v in zip(variables, line[2:]):
                            out.e1 = section.strip("# ")
                            out.e2 = k
                            out.val = v
                            print(out)

