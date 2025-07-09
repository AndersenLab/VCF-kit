#! /usr/bin/env python
"""
usage:
  vk phylo fasta <vcf> [<region>]
  vk phylo tree (nj|upgma) [--plot] <vcf> [<region>]

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
import os
import sys
import tempfile
import webbrowser
from pkgutil import get_data
from subprocess import PIPE, Popen

import numpy as np
from jinja2 import Template
import biotite.sequence.align as align
import biotite.sequence.phylo as phylo
from biotite.sequence import NucleotideSequence
import networkx as nx
import matplotlib.pyplot as plt

from clint.textui import colored, indent, puts_err
from docopt import docopt
from vcfkit import __version__
from vcfkit.utils import check_program_exists
from vcfkit.utils.vcf import vcf


def main(debug=None):
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False,
                  version=__version__)


    def first(s):
        return s[0].replace(".", "N")

    firstv = np.vectorize(first)

    v = vcf(args["<vcf>"])

    if len(v.samples) <= 1:
        exit(puts_err(colored.red("\n\tVCF must have at least two samples.\n")))
    
    if args["<region>"]:
        variant_set = v(args["<region>"])
    else:
        variant_set = v

    if args["fasta"] or args["tree"]:
        """
            Generate an aligned fasta from a VCF file.
        """
        gt_set = []
        for line in variant_set:
            if line.is_snp:
                gt_set.append(firstv(line.gt_bases))
        if len(gt_set) == 0:
            exit(puts_err("No genotypes"))
        gt_set = np.vstack(gt_set)
        seqs = [''.join(gt_set[:, i]) for i in range(gt_set.shape[1])]

        if args["fasta"]:
            for sample, seq in list(zip(v.samples, seqs)):
                print(">" + sample)
                print(seq)

        elif args["tree"]:
            """
            Generate a phylogenetic tree using an aligned fasta with muscle.
            """

            with indent(4):
                puts_err(colored.blue("\nGenerating Fasta\n"))
            trace = align.Alignment.trace_from_strings(seqs)
            aligned_seqs = align.Alignment([NucleotideSequence(seq) for seq in seqs], trace)
            distances = 1 - align.get_pairwise_sequence_identity(aligned_seqs, mode="all")

            if args["nj"]:
                tree = phylo.neighbor_joining(distances)
            else:
                tree = phylo.upgma(distances)

            # output tree
            print(tree.to_newick(labels=v.samples))
            
            if args["--plot"]:
                # graph = tree.as_graph().to_undirected()
                # fig = plt.figure(figsize=(8.0, 8.0))
                # ax = fig.gca()
                # ax.axis("off")
                # # Calculate position of nodes in the plot
                # pos = nx.kamada_kawai_layout(graph)
                # # Assign the gene names to the nodes that represent a reference index
                # node_labels = {i: name for i, name in enumerate(v.samples)}
                # nx.draw_networkx_edges(graph, pos, ax=ax)
                # nx.draw_networkx_labels(
                #     graph,
                #     pos,
                #     ax=ax,
                #     labels=node_labels,
                #     font_size=7,
                #     # Draw a white background behind the labeled nodes
                #     # for better readability
                #     bbox=dict(pad=0, color="white"),
                # )
                # fig.tight_layout()

                # plt.show()

                prefix = os.path.dirname(os.path.abspath(sys.modules['vcfkit'].__file__)) + "/static"
                template = open(prefix + "/tree.html",'r').read()
                tree_template = Template(template)
                html_out = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
                with html_out as f:
                    tree = tree.to_newick(labels=v.samples)
                    sample_len = len(v.samples)
                    f.write(tree_template.render(tree=tree, prefix=prefix, sample_len=sample_len).encode())
                    webbrowser.open("file://" + html_out.name)

if __name__ == '__main__':
    main()
