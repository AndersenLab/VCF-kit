#! /usr/bin/env python
"""
usage:
  vk hmm [--alt=<alt_sample> --vcf-out] <vcf>

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
from itertools import groupby
from operator import itemgetter
signal(SIGPIPE, SIG_DFL)

# Setup hmm
model = Model(name="RIL_GT")

ref = State(DiscreteDistribution({'ref': 0.97, 'alt': 0.03}), name='ref')
alt = State(DiscreteDistribution({'ref': 0.03, 'alt': 0.97}), name='alt')

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
    debug = ["hmm", "--alt=CB4856", "--vcf-out", "../test.vcf.gz"]


def generate_cigar(arr):
    grouped = [(k, sum(1 for i in g)) for k,g in groupby(arr)]
    return "".join([{0: "R", 1: "A"}[x] + str(y) for x,y in grouped]), len(grouped) - 1


class ranges:
    """
        Class for storing genomic ranges of genotypes.
    """

    def process_results(self, results, out=False):
        if out is False:
            n, self.dp_sum, self.site_count, self.supporting_site_count = 0, 0, 0, 0
            dp_avg = ""
            last_contig = "null"
            last_pred = -8
            for chrom_set in [list(g) for k, g in groupby(results, itemgetter(0))]:
                chrom = chrom_set[0][0]
                for interval in [list(g) for k, g in groupby(chrom_set, itemgetter(3))]:
                    orig = [x[2] for x in interval]
                    pred = [x[3] for x in interval]
                    gt = pred[0]
                    orig_cigar, switches = generate_cigar(orig)
                    supporting_sites = len([x for x in interval if x[2] == x[3]])
                    dp_avg = 0
                    if supporting_sites > 0:
                        dp_avg = sum([x[4] for x in interval if x[2] == x[3]])*1.0/supporting_sites
                    gt = interval[0][3]
                    start = min([x[1] for x in interval])
                    end = max([x[1] for x in interval])
                    sites = len(interval)
                
                    output = [chrom,
                              start,
                              end,
                              sample,
                              gt + 1,
                              supporting_sites,
                              sites,
                              dp_avg, 
                              switches,
                              orig_cigar]
                    out_line = "\t".join(map(str,output))
                    print(out_line)



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
    gt_list, dp_list, dv_list = [], [], []
    chromosome = []
    positions = []

    append_gt = True
    for n, line in enumerate(v):
        if args["--alt"]:
            # Check if gt is 1/1 for alt sample.
            append_gt = (line.gt_types[v.samples.index(args["--alt"])] == 3)
        if append_gt:
            chromosome.append(line.CHROM)
            positions.append(line.POS)
            gt_list.append(np.array(line.gt_types))
            dp_list.append(np.array(line.format("DP", int)))
    gt_set = np.vstack([gt_list]).astype("int8").T

    dp_set = np.vstack([dp_list]).astype("int32").T
    gt_set[gt_set == 1] = -9  # Het
    gt_set[gt_set == 2] = -9  # Missing
    gt_set[gt_set == 3] = 1  # Alt

    dp_set[dp_set < 0] = 0

    # Result dict
    vcf_gt_out = defaultdict(lambda: defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))))
    sample_gt_matrix = np.array([])

    chrom_to_factor = {y: x for x, y in enumerate(set(chromosome))}
    sample_to_factor = {y: x for x, y in enumerate(v.samples)}
    hmm_gt_set = []
    gtr = ranges()

    if args["--vcf-out"] is False:
        print("chrom\tstart\tend\tsample\tgt\tsupporting_sites\tsites\tDP\tswitches\tCIGAR")

    s = 0
    for sample, column in zip(v.samples, gt_set):
        sample_gt = zip(chromosome, positions, column)
        sample_gt = [x for x in sample_gt if x[2] in [0, 1]]
        sequence = [to_model[x[2]] for x in sample_gt]
        gtr.sample = sample
        if sequence:
            results = model.forward_backward(sequence)[1]
            if model.states[0].name == 'ref':
                result_gt = [from_model[x]
                             for x in np.greater(results[:, 1], results[:, 0])]
            else:
                result_gt = [from_model[x]
                             for x in np.greater(results[:, 0], results[:, 1])]
            dp_s = dp_set[0][s]
            dp_s = dp_s[dp_s > 0]
            results = ((a, b, c, d, e)
                       for (a, b, c), d, e in zip(sample_gt, result_gt, dp_s))
            gtr.s = s

            results = list(results)

            gtr.process_results(results, args["--vcf-out"])
        s += 1

    if args["--vcf-out"]:
        v = vcf(args["<vcf>"])
        for n, line in enumerate(v):
            line = variant_line(line)
            if line.has_gt and line.chrom in chromosome and line.pos in positions:
                for sample_col, sample in enumerate(v.samples):
                    gt_orig = line[sample]
                    line.modify_gt_format(sample_col, "GT_ORIG", gt_orig)
                    new_gt = gtr.get(line.chrom, line.pos, sample)
                    if new_gt is not None:
                        line.modify_gt_format(sample_col, "GT", to_gt[new_gt])
                print(line)
