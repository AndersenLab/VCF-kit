#! /usr/bin/env python
"""
usage:
  tb hmm [--alt=<alt_sample> --vcf-out] <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
# Suppress the rocket ship!
import matplotlib
matplotlib.use("Agg")
from utils.vcf import *
from utils.fasta import *
from collections import defaultdict
import sys
from utils import autoconvert
from yahmm import *
import itertools
import numpy as np
from pprint import pprint as pp
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

# Setup hmm
model = Model(name="RIL_GT")

ref = State(DiscreteDistribution({'ref': 0.97, 'alt': 0.03}), name = 'ref')
alt = State(DiscreteDistribution({'ref': 0.03, 'alt': 0.97}), name = 'alt')

model.add_transition(model.start, ref, 0.5)
model.add_transition(model.start, alt, 0.5)

# Transition matrix, with 0.05 subtracted from each probability to add to
# the probability of exiting the hmm
model.add_transition(ref, ref, 0.94)
model.add_transition(ref, alt, 0.005)
model.add_transition(alt, ref, 0.005)
model.add_transition(alt, alt, 0.94)

model.add_transition(ref, model.end, 0.01)
model.add_transition(alt, model.end, 0.01)

model.bake(verbose=False)

to_model = {0: 'ref', 1: 'alt'}
to_gt = {0: '0/0', 1: '1/1'}
from_model = {True: 1, False: 0}

debug = None
if len(sys.argv) == 1:
    debug = ["hmm", "--alt=CB4856", "--vcf-out",  "../test.vcf.gz"]

class ranges:
    """
        Class for storing genomic ranges of genotypes.
    """
    def __init__(self):
        self.range = defaultdict(list)

    def add(self, chrom, start, end, sample, gt, site_count, out = False):
        if out is False:
            print('\t'.join(map(str, [chrom, start, end, sample, gt + 1, site_count])))
        else:
            self.range[sample].append([chrom, int(start), int(end), int(gt)])


    def get(self, qchrom, qpos, qsample):
        gt_found = [x for x in self.range[sample] if x[0] == qchrom and qpos >= x[1] and qpos <= x[2]]
        if gt_found:
            return gt_found[0][3]
        else:
            return None


if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False)

    v = vcf(args["<vcf>"])

    # Put genotypes into large array
    """
        0/0 -> 0
        0/1 -> 1
        ./. -> 2
        1/1 -> 3
    """
    gt_list = []
    chromosome = []
    positions = []

    if args["--alt"]:
        # Subset by sample
        alt_sample = args["--alt"]
        alt_sample_index = v.samples.index(alt_sample)
    append_gt = True
    for n, line in enumerate(v):
        if args["--alt"]:
            append_gt = (line.gt_types[alt_sample_index] == 3)
        if append_gt:
            chromosome.append(line.CHROM)
            positions.append(line.POS)
            gt_list.append(list(line.gt_types))
    gt_set = np.vstack([gt_list]).astype("int32")
    gt_set[gt_set == 1] = -9  # Het
    gt_set[gt_set == 2] = -9  # Missing
    gt_set[gt_set == 3] = 1  # Alt
    # Result dict
    vcf_gt_out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    sample_gt_matrix = np.array([])

    chrom_to_factor = {y: x for x, y in enumerate(set(chromosome))}
    sample_to_factor = {y: x for x, y in enumerate(v.samples)}
    hmm_gt_set = []
    gtr = ranges()

    if args["--vcf-out"] is False:
        print("chrom\tstart\tend\tsample\tgt\tsites")

    for sample, column in zip(v.samples, gt_set.T):
        sample_gt = zip(chromosome, positions, column)
        sample_gt = [x for x in sample_gt if x[2] in [0, 1]]
        sequence = [to_model[x[2]] for x in sample_gt]
        if sequence:
            results = model.forward_backward(sequence)[1]
            if model.states[0].name == 'ref':
                result_gt = [from_model[x] for x in np.greater(results[:, 1], results[:, 0])]
            else:
                result_gt = [from_model[x] for x in np.greater(results[:, 0], results[:, 1])]
            results = ((a, b, c, d) for (a, b, c), d in zip(sample_gt, result_gt))
            n = 0
            site_count = 0
            last_contig = "null"
            last_pred = -8
            for chrom, pos, orig, pred in results:
                site_count += 1
                if chrom != last_contig or pred != last_pred:
                    if n > 0:
                        gtr.add(contig, start, end, sample, last_pred, site_count, args["--vcf-out"] )
                        site_count = 0
                    contig = chrom
                    start = pos
                end = pos
                last_pred = pred
                last_contig = chrom
                last_pos = pos
                n += 1
            gtr.add(contig, start, end, sample, last_pred, site_count, args["--vcf-out"] )

    if args["--vcf-out"]:
        v = vcf(args["<vcf>"])
        print(v.insert_header_line("##FORMAT=<ID=GT_ORIG,Number=1,Type=String,Description=\"Original Genotype replaced by HMM\">"))
        for n, line in enumerate(v):
            line = variant_line(line)
            if line.has_gt and line.chrom in chromosome and line.pos in positions:
                for sample_col, sample in enumerate(v.samples):
                    gt_orig = line.fetch_gt_from_index(sample_col)
                    line.modify_gt_format(sample_col, "GT_ORIG", gt_orig)
                    new_gt = gtr.get(line.chrom, line.pos, sample)
                    if new_gt is not None:
                        line.modify_gt_format(sample_col, "GT", to_gt[new_gt])
                print(line)
