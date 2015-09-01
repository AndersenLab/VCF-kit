#! /usr/bin/env python
"""
usage:
  tb.py tajima [--no-header] <window-size> <step-size> <vcf> 
  tb.py tajima [--no-header] <window-size> --sliding <vcf> 

command:
  tajima        Calculate Tajima's D

options:
  -h --help                   Show this screen.
  --version                   Show version.  
  --window-size               blah
  --step-size                 blash 

output:
    CHROM
    BIN_START
    BIN_END
    N_Sites
    N_SNPs
    TajimaD


"""
from docopt import docopt
from subprocess import call, Popen, PIPE
from itertools import combinations
from utils.vcf import *
from math import isinf
import sys


debug = None



class tajima(vcf):
    """
      Subclass of the vcf object
      used for calculating tajima's D
    """
    def __init__(self, filename):
      vcf.__init__(self,filename)

    def calc_tajima(self, window_size, step_size, sliding = False):
        # Tajima D Constants
        n = self.n*2
        a1 = sum([1.0/i for i in xrange(1,n)])
        a2 = sum([1.0/i**2 for i in xrange(1,n)])
        b1 = (n + 1.0) / (3.0*(n - 1))
        b2 = (2.0 * (n**2 + n + 3.0)) / \
           (9.0 * n * (n -1.0))
        c1 = b1 - (1.0 / a1)
        c2 = b2 - ((n + 2.0) / (a1*n)) + (a2 / (a1**2) )
        e1 = c1 / a1
        e2 = c2 / (a1**2 + a2)

        if args["--sliding"]:
            shift_method = "POS-Sliding"
            step_size = None
        else:
            shift_method = "POS-Interval"

        for variant_interval in self.window(window_size= window_size, step_size = step_size, shift_method=shift_method):
            pi = 0.0
            S = 0
            n_sites = 0
            for variant in variant_interval:
                n_sites += 1
                AN = variant.INFO.get("AN")  # c;AN : total number of alleles in called genotypes
                AC = variant.INFO.get("AC")  # j;AC : allele count in genotypes, for each ALT allele, in the same order as listed
                try:
                    # Biallelic sites only!
                    pi += (2.0 * AC * (AN-AC)) / (AN * (AN-1.0))
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
                    yield CHROM, variant_interval.lower_bound, variant_interval.upper_bound, n_sites, S, TajimaD
            except:
                pass



if __name__ == '__main__':
    print sys.argv
    args = docopt(__doc__, 
                  version='VCF-Toolbox v0.1',
                  argv = debug)
    print args
    if args["<vcf>"] == "":
      print(__doc__)
    wz = int(args["<window-size>"])
    sz = None
    if not args["--sliding"]:
      sz = int(args["<step-size>"])
    if args["--no-header"] == False:
        print("CHROM\tBIN_START\tBIN_END\tN_Sites\tN_SNPs\tTajimaD")
    for i in tajima(args["<vcf>"]).calc_tajima(wz,sz):
        print("\t".join([str(x) for x in i] ))



