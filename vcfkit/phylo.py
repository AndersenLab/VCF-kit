#! /usr/bin/env python
"""
usage:
  vk phylo fasta <vcf>
  vk phylo tree (nj|upgma) [--plot] <vcf>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
import vk
from vcfkit import __version__
from utils.vcf import *
from subprocess import Popen, PIPE
from utils import check_program_exists
from clint.textui import colored, indent, puts_err
import os


def main(debug=None):
    args = docopt(__doc__,
                  argv=debug,
                  options_first=False,
                  version=__version__)
    module_path = os.path.split(os.path.realpath(__file__))[0]
    v = vcf(args["<vcf>"])
    samples = v.samples
    _ROOT = os.path.split(os.path.dirname(vk.__file__))[0]
    if args["fasta"] or args["tree"]:
        """
            Generate an aligned fasta from a VCF file.
        """
        seqs = {}
        for sample in samples:
            seqs[sample] = []
        for line in v:
            if line.is_snp:
                non_missing = [x.replace(".", "-") for x in line.gt_bases]
                sample_gt = zip(samples, [x[-1] for x in non_missing])
            for sample, gt in sample_gt:
                seqs[sample].append(gt)
        if not args["tree"]:
            for sample, seq in seqs.items():
                print(">" + sample)
                print(''.join(seq))
        elif args["tree"]:
            """
            Generate a phylogenetic tree using an aligned fasta with muscle.
            """

            # Check for muscle dependency
            check_program_exists("muscle")
            fasta = ""
            with indent(4):
                puts_err(colored.blue("\nGenerating Fasta\n"))
            for sample, seq in seqs.items():
                fasta += ">" + sample + "\n" + ''.join(seq) + "\n"
            tree_type = "upgma"  # default is upgma
            if args["nj"]:
                tree_type = "neighborjoining"
            with indent(4):
                puts_err(colored.blue("\nGenerating " + tree_type + " Tree\n"))
            comm = ["muscle", "-maketree", "-in", "-", "-cluster", tree_type]
            tree, err = Popen(comm, stdin=PIPE, stdout=PIPE).communicate(input=fasta)
            print(tree)
            if args["--plot"]:
                from jinja2 import Template
                import webbrowser
                import tempfile
                # R code for plotting here!
                prefix = _ROOT + "/static"
                tree_template = Template(open(_ROOT + "/static/tree.html", 'r').read())
                html_out = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
                with html_out as f:
                    tree = tree.replace("\n", "")
                    sample_len = len(samples)
                    f.write(tree_template.render(**locals()))
                    # print html_out.name
                    webbrowser.open("file://" + html_out.name)

if __name__ == '__main__':
    main()
