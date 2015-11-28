#! /usr/bin/env python
"""
usage:
  tb.py freq <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.



"""
from docopt import docopt
from utils.vcf import *
from utils.fasta import *
from collections import defaultdict
import sys


class freq_vcf(vcf):

    """
        Subclass of vcf that calculates frequency of alleles by strain
    """

    def __init__(self, filename):
        vcf.__init__(self, filename)

    def calc_af(self):
        af_freq = {}
        af_freq = {sample: defaultdict(int) for sample in self.samples}
        for line in self:
            """
            0/0 -> 0
            0/1 -> 1
            ./. -> 2
            1/1 -> 3
            """
            for i in [sample for sample, gt in zip(self.samples, line.gt_types) if gt == 3]:
                af_freq[i][line.num_hom_alt] += 1
        # Output results
        print "\t".join(["sample", "freq_of_gt", "n_gt_at_freq"])
        for sample in af_freq.keys():
            for i in xrange(1, len(self.samples) + 1):
                out = "\t".join(map(str, [sample, i, af_freq[sample][i]]))
                print(out)


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
    freq = freq_vcf(args["<vcf>"])
    freq.calc_af()
