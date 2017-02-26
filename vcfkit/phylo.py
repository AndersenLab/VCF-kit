#! /usr/bin/env python
"""
usage:
  vk phylo fasta <vcf> [<region>]
  vk phylo tree (nj|upgma) [--plot] <vcf> [<region>]

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from vcfkit import __version__
from utils.vcf import *
from subprocess import Popen, PIPE
from utils import check_program_exists
from clint.textui import colored, indent, puts_err
import os
from pkgutil import get_data
import sys
import numpy as np

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
        gt_set = np.chararray((0,len(v.samples)))
        gt_set = []
        for line in variant_set:
            if line.is_snp:
                gt_set.append(firstv(line.gt_bases))
        if len(gt_set) == 0:
            exit(puts_err("No genotypes"))
        gt_set = np.vstack(gt_set)
        seqs = zip(v.samples, np.transpose(gt_set))
        if args["fasta"]:
            for sample, seq in seqs:
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
            for sample, seq in seqs:
                fasta += ">" + sample + "\n" + ''.join(seq) + "\n"
            tree_type = "upgma"  # default is upgma
            if args["nj"]:
                tree_type = "neighborjoining"
            with indent(4):
                puts_err(colored.blue("\nGenerating " + tree_type + " Tree\n"))
            comm = ["muscle", "-maketree", "-in", "-", "-cluster", tree_type]
            tree, err = Popen(comm, stdin=PIPE, stdout=PIPE).communicate(input=fasta)
            
            # output tree
            print(tree)
            
            if args["--plot"]:
                from jinja2 import Template
                import webbrowser
                import tempfile
                prefix = os.path.dirname(os.path.abspath(sys.modules['vcfkit'].__file__)) + "/static"
                template = open(prefix + "/tree.html",'r').read()
                tree_template = Template(template)
                html_out = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
                with html_out as f:
                    tree = tree.replace("\n", "")
                    sample_len = len(v.samples)
                    f.write(tree_template.render(**locals()))
                    webbrowser.open("file://" + html_out.name)

if __name__ == '__main__':
    main()
