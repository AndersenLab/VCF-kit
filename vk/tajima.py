#! /usr/bin/env python
"""
usage:
  tb tajima [--no-header --extra] <window-size> <step-size> <vcf> 
  tb tajima [--no-header --extra] <window-size> --sliding <vcf> 

options:
  -h --help                   Show this screen.
  --version                   Show version.  
  --window-size               blah
  --step-size                 blash 
  --extra                     display extra

command:
  tajima        Calculate Tajima's D


output:
    CHROM
    BIN_START
    BIN_END
    N_Sites
    N_SNPs
    TajimaD


"""
from docopt import docopt
from utils.vcf import *
from math import isinf
import sys
import os
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


debug = None


class tajima(vcf):

    """
      Subclass of the vcf object
      used for calculating tajima's D
    """

    def __init__(self, filename):
        vcf.__init__(self, filename)

    def calc_tajima(self, window_size, step_size, extra=False):
        # Tajima D Constants
        n = self.n * 2
        a1 = sum([1.0 / i for i in xrange(1, n)])
        a2 = sum([1.0 / i**2 for i in xrange(1, n)])
        b1 = (n + 1.0) / (3.0 * (n - 1))
        b2 = (2.0 * (n**2 + n + 3.0)) / \
            (9.0 * n * (n - 1.0))
        c1 = b1 - (1.0 / a1)
        c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1**2))
        e1 = c1 / a1
        e2 = c2 / (a1**2 + a2)

        if args["--sliding"]:
            shift_method = "POS-Sliding"
            step_size = None
        else:
            shift_method = "POS-Interval"

        for variant_interval in self.window(window_size=window_size, step_size=step_size, shift_method=shift_method):
            pi = 0.0
            S = 0
            n_sites = 0
            for variant in variant_interval:
                n_sites += 1
                AN = variant.INFO.get("AN")  # c;AN : total number of alleles in called genotypes
                # j;AC : allele count in genotypes, for each ALT allele, in the same order as listed
                AC = variant.INFO.get("AC")
                try:
                    # Biallelic sites only!
                    pi += (2.0 * AC * (AN - AC)) / (AN * (AN - 1.0))
                    S += 1
                except:
                    pass
            try:
                CHROM = variant_interval[0].CHROM
                tw = (S / a1)
                var = (e1 * S) + ((e2 * S) * (S - 1.0))
                TajimaD = (pi - tw) / \
                          (var)**(0.5)
                if not isinf(TajimaD) and S > 0:
                    output = [CHROM,
                              variant_interval.lower_bound,
                              variant_interval.upper_bound,
                              n_sites,
                              S,
                              TajimaD]
                    if extra:
                        if step_size is None:
                            step_size = "NA"
                        output += [os.path.split(self.filename)[1],
                                   window_size,
                                   step_size]
                    yield "\t".join(map(str, output))
            except:
                pass

if len(sys.argv) == 1:
    debug = ["tajima", "100000", "10000",
             "~/Dropbox/AndersenLab/wormreagents/Variation/Andersen_VCF/20150731_WI_PASS.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug)
    if args["<vcf>"] == "":
        print(__doc__)
    wz = int(args["<window-size>"].replace(",", ""))
    sz = None
    if not args["--sliding"]:
        sz = int(args["<step-size>"].replace(",", ""))
    if args["--no-header"] is False:
        header_line = ["CHROM",
                       "BIN_START",
                       "BIN_END",
                       "N_Sites",
                       "N_SNPs",
                       "TajimaD"]
        if args["--extra"]:
            header_line += ["filename",
                            "window_size",
                            "step_size"]
        print "\t".join(header_line)
    for i in tajima(args["<vcf>"]).calc_tajima(wz, sz, extra=args["--extra"]):
        print(i)
