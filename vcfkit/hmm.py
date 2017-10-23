#! /usr/bin/env python
"""
usage:
  vk hmm [options] --A=<parent_1> --B=<parent_2> <vcf>

Example

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --vcf-out                   Output VCF instead of intervals.
  --endfill                   Don't leave gaps at the ends of chromosomes.
  --infill                    Fill in missing portions.
  --state=<state>             State probability [default: 0.97].
  --transition=<transition>   Transition probability [default: 1e-9]


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
from pomegranate import *
import itertools
import numpy as np
from signal import signal, SIGPIPE, SIG_DFL
from itertools import groupby
from operator import itemgetter
from intervaltree import IntervalTree
from vk import __version__
from clint.textui import colored, puts_err, indent
signal(SIGPIPE, SIG_DFL)


def generate_model(state, transition):
    # Setup hmm
    model = HiddenMarkovModel()

    A = State(DiscreteDistribution({'A': state, 'B': 1-state}), name='A')
    B = State(DiscreteDistribution({'A': 1-state, 'B': state}), name='B')

    model.add_transition(model.start, A, 0.5)
    model.add_transition(model.start, B, 0.5)


    model.add_transition(A, A, 1-transition)
    model.add_transition(A, B, transition)
    model.add_transition(B, A, transition)
    model.add_transition(B, B, 1-transition)

    model.add_transition(A, model.end, 0.5)
    model.add_transition(B, model.end, 0.5)

    model.bake(verbose=False)
    return model

to_model = {0: 'A', 1: 'B'}
to_gt = {0: '0/0', 1: '1/1'}
from_model = {True: 1, False: 0}

debug = None
if len(sys.argv) == 1:
    debug = ["hmm","--A=REF", "--B=ALT", "--vcf-out", "../test.vcf.gz"]


def generate_RLE(arr):
    grouped = [(k, sum(1 for i in g)) for k, g in groupby(arr)]
    return "".join([{0: "A", 1: "B"}[x] + str(y) for x, y in grouped]), len(grouped) - 1


if __name__ == '__main__':
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False,
                  version=__version__)
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
    p1_gt_set = []
    p2_gt_set = []

    append_gt = True
    for n, line in enumerate(v):
        ref_alt = ["REF", "ALT"]
        if args["--A"] == "REF":
            p1_gt = 0
        else:
            p1_gt = line.gt_types[v.samples.index(args["--A"])]
        if args["--B"] == "ALT":
            p2_gt = 3
        else:
            p2_gt = line.gt_types[v.samples.index(args["--B"])]

        if p1_gt != p2_gt and p1_gt in [0,3] and p2_gt in [0,3]:
            chromosome.append(line.CHROM)
            positions.append(line.POS)
            gt_list.append(np.array(line.gt_types))
            dp_list.append(np.array(line.format("DP", int)))
            p1_gt_set.append(p1_gt)
            p2_gt_set.append(p2_gt)
    gt_set = np.vstack([gt_list]).astype("int8").T

    dp_set = np.vstack([dp_list]).astype("int32").T
    gt_set[gt_set == 1] = -9  # Het
    gt_set[gt_set == 2] = -9  # Missing
    
    gt_set[gt_set == 3] = 1  # Alt GT
    p1_gt_set = np.array(p1_gt_set)
    p2_gt_set = np.array(p2_gt_set)
    p1_gt_set[p1_gt_set == 3] = 1
    p2_gt_set[p2_gt_set == 3] = 1

    dp_set[dp_set < 0] = 0

    # Result dict
    vcf_gt_out = defaultdict(lambda: defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))))
    sample_gt_matrix = np.array([])

    chrom_to_factor = {y: x for x, y in enumerate(set(chromosome))}
    sample_to_factor = {y: x for x, y in enumerate(v.samples)}
    hmm_gt_set = []

    if args["--vcf-out"] and args["<vcf>"] == "-":
        with indent(2):
            exit(puts_err(colored.blue("Cannot use vcf-out with stdin.")))

    if not args["--vcf-out"]:
        print("chrom\tstart\tend\tsample\tgt\tgt_name\tsupporting_sites\tsites\tDP\tswitches\trle")

    s = 0
    tree = {}
    for sample, column in zip(v.samples, gt_set):
        sample_gt = zip(chromosome, positions, column == p2_gt_set, column)
        sample_gt = [x[:3] for x in sample_gt if x[3] in [0, 1]]
        sequence = [to_model[x[2]] for x in sample_gt]

        model = generate_model(float(args['--state']), float(args['--transition']))
        if sequence:
            results = model.forward_backward(sequence)[1]
            if model.states[0].name == 'A':
                result_gt = [from_model[x]
                             for x in np.greater(results[:, 1], results[:, 0])]
            else:
                result_gt = [from_model[x]
                             for x in np.greater(results[:, 0], results[:, 1])]
            dp_s = dp_set[0][s]
            dp_s = dp_s[dp_s > 0]
            results = ((a, b, c, d, e)
                       for (a, b, c), d, e in zip(sample_gt, result_gt, dp_s))

            results = list(results)

            dp_avg = ""
            chrom_sets = [list(g) for k, g in groupby(results, itemgetter(0))]
            tree[sample] = {}
            for chrom_set in chrom_sets:
                end = 0
                chrom = chrom_set[0][0]
                tree[sample][chrom] = IntervalTree()
                interval_set = [list(g) for k, g in groupby(chrom_set, itemgetter(3))]
                interval_set_len = len(interval_set)
                for interval_n, interval in enumerate(interval_set):
                    orig = [x[2] for x in interval]
                    pred = [x[3] for x in interval]
                    gt = pred[0]
                    orig_RLE, switches = generate_RLE(orig)
                    supporting_sites = len([x for x in interval if x[2] == x[3]])
                    dp_avg = 0
                    if supporting_sites > 0:
                        dp_avg = sum([x[4] for x in interval if x[2] == x[3]])*1.0/supporting_sites
                    gt = interval[0][3]
                    if supporting_sites > 0:
                        if args["--infill"]:
                            start = end + 1
                        else:
                            start = min([x[1] for x in interval])
                        end = max([x[1] for x in interval])
                        sites = len(interval)

                        # End Fill
                        if args["--endfill"] and interval_n == 0:
                            start = 1
                        if args["--endfill"] and interval_n == interval_set_len - 1:
                            end = v.contigs[chrom]
                        gt_name = {0: args['--A'], 1: args['--B']}[gt]
                        output = [chrom,
                                  start,
                                  end,
                                  sample,
                                  gt + 1,
                                  gt_name,
                                  supporting_sites,
                                  sites,
                                  dp_avg,
                                  switches,
                                  orig_RLE]
                        if not args["--vcf-out"]:
                            out_line = "\t".join(map(str, output))
                            print(out_line)
                        else:
                            tree[sample][chrom][start:end+1] = gt
        s += 1

    if args["--vcf-out"]:
        v = vcf(args["<vcf>"])
        chrom_pos = [x[0:2] for x in sample_gt]
        # Add GT_ORIG and print raw header.
        v.add_format_to_header({"ID": "GT_ORIG",
                                "Description": "Original genotype",
                                "Type": "Character",
                                "Number": "1"})
        print(v.raw_header.strip())
        for n, line in enumerate(v):
            line = variant_line(line, v.samples)
            for sample_n, sample in enumerate(v.samples):
                gt_orig = line.get_gt("GT", sample_n)
                try:
                    new_gt = next(iter(tree[sample][line.chrom].search(line.pos))).data
                except:
                    new_gt = None
                if new_gt is not None:
                    line.set_gt("GT_ORIG", sample_n, gt_orig)
                    line.set_gt("GT", sample_n, to_gt[new_gt])
            print(line)
