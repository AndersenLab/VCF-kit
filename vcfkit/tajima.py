#! /usr/bin/env python
"""
usage:
  vk tajima [--no-header --extra] <window-size> <step-size> <vcf>
  vk tajima [--no-header --extra] <window-size> --sliding <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --window-size               Size of window
  --step-size                 Distance to move window.
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
from vcfkit import __version__
from docopt import docopt
from utils.vcf import *
from math import isinf
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

    def calc_tajima(self, window_size, step_size, sliding, extra=False):
        # Tajima D Constants
        n = self.n * 2
        a1 = sum([1.0 / i for i in xrange(1, n)])
        a2 = sum([1.0 / i**2.0 for i in xrange(1, n)])
        b1 = (n + 1.0) / (3.0 * (n - 1.0))
        b2 = (2.0 * (n**2 + n + 3.0)) / \
            (9.0 * n * (n - 1.0))
        c1 = b1 - (1.0 / a1)
        c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1**2))
        e1 = c1 / a1
        e2 = c2 / ((a1*a1) + a2)

        if sliding:
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
                # Use only diploid, biallelic sites, 0 < AF < 1
                if variant.ploidy == 2 and len(variant.ALT) == 1 and 0 < variant.aaf < 1 and variant.is_snp:
                    AN = variant.num_called*2
                    # j;AC : allele count in genotypes, for each ALT allele, in the same order as listed
                    AC = variant.INFO["AC"]
                    if AC > 0:
                        p = (AC*1.0/AN)*1.0
                        pi += p * (1.0-p)
                        S += 1
            try:
                CHROM = variant_interval[0].CHROM
                pi = 2.0*pi*n/(n-1.0)
                tw = (S*1.0 / a1)
                var = (e1 * S) + ((e2 * S) * (S - 1))
                TajimaD = (pi - tw) / \
                          (var)**(0.5)
                if not isinf(TajimaD) and S > 0 and AC > 0:
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


def large_int(i):
    return int(float(i.replace(",", "")))


def main(debug=None):
    args = docopt(__doc__,
                  version=__version__,
                  argv=debug)

    if args["<vcf>"] == "":
        print(__doc__)

    wz = large_int(args["<window-size>"])
    sz = None
    if not args["--sliding"]:
        sz = large_int(args["<step-size>"])

    if wz < sz:
        exit(puts_err(colored.red("\n\tWindow size must be >= step size.\n")))

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
    for i in tajima(args["<vcf>"]).calc_tajima(wz, sz, args["--sliding"], extra=args["--extra"]):
        print(i)

if __name__ == '__main__':
    main()
