#! /usr/bin/env python
"""
usage:
  tb freq sample_hom_gt <vcf>
  tb freq spectrum (alleles|genotypes) [--rate] <vcf>

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
from utils import autoconvert

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
            for i in xrange(1, len(self.samples)+1):
                out = "\t".join(map(str, [sample, i, af_freq[sample][i]]))
                print(out)

    def calc_frequency(self):
        s_freq = defaultdict(int)
        if args["genotypes"]:
            """
                Calculate number of given genotypes
            """
            for line in self:
                ref = line.num_hom_ref
                alt = line.num_hom_alt
                het = line.num_het
                mis = line.num_unknown
                if args["--rate"]:
                    ref = float(ref)/vcf.n
                    het = float(het)/vcf.n
                    alt = float(alt)/vcf.n
                    mis = float(mis)/vcf.n
                s_site = '__'.join(map(str, [ref, het, alt, mis]))
                s_freq[s_site] += 1
            freq_set = [(int(n), map(autoconvert, x.split("__"))) for n, x in zip(s_freq.values(), s_freq.keys())]
            freq_set = sorted(freq_set, key=lambda s: s[0], reverse=True)
            print("n\tref\thet\talt\tmis")
            for n, freqs in freq_set:
                print(str(n) + '\t' + '\t'.join(map(str, freqs)))
        elif args["alleles"]:
            """
                Calculate ALT Allele Frequency Spectrum
            """
            print("n\talt_allele_freq")
            for line in self:
                s_freq[line.aaf] += 1
            s_freq = sorted(s_freq.items(), key=lambda s: s[0], reverse=True)
            for freq, n in s_freq:
                print(str(n) + "\t" + str(freq))




debug = None
if len(sys.argv) == 1:
    debug = ["freq", "spectrum", "alleles", "../test.vcf.gz"]


if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)

    vcf = freq_vcf(args["<vcf>"])
    if args["sample_hom_gt"]:
        vcf.calc_af()
    elif args["spectrum"]:
        vcf.calc_frequency()

    