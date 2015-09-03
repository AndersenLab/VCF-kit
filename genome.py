#! /usr/bin/env python
"""
usage:
  tb.py genome <vcf>
  tb.py genome --search=<term>
  tb.py genome --setup=<asm_name>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from utils.vcf import *
import sys
from clint.textui import colored, puts, indent
from Bio.Entrez import esearch, read, efetch, esummary
from time import time
from tabulate import tabulate as tb
import requests
import os


debug = None
if len(sys.argv) == 1:
    debug = ['genome', "test.vcf.gz"]


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)
    print args
    if args["<vcf>"] is not None:
        v = vcf(args["<vcf>"])
        with indent(2):
            puts(colored.blue('\nReference:\n'))
        with indent(4):
            puts(colored.green(v.metadata["reference"] + "\n"))

    elif args["--search"]:
        # Download and cache a list of genomes from NCBI for searching
        genome_file_path = os.path.split(os.path.realpath(__file__))[0] + "/genomes.txt"
        if os.path.isfile(genome_file_path):
            fileTime = os.path.getctime(genome_file_path)
        else:
            fileTime = 0
        if fileTime < (time() - 60*60*24*2):
            with indent(2):
                puts(colored.blue('\nDownloading list of reference genomes\n'))
            r = requests.get("http://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt")
            genome_file = open(genome_file_path, "w")
            with genome_file as f:
                f.write(r.text.encode('utf-8').strip())

        # Cache result
        header = ["assembly_accession", # 0
                 "bioproject", # 1
                 "organism_name", # 7
                 "asm_name", # 15
                 "ftp_path"] # 19

        with indent(2):
            puts(colored.blue('\nSearching...\n'))

        with open(genome_file_path, "r") as f:
            results = []
            for line in f:
                line = line.strip().split("\t")
                line = [x for k,x in enumerate(line) if k in [0,1,7,15,19]]
                if args["--search"].lower() in line[2].lower() and line[4] != "na":
                    results.append(line)
        with indent(4):
            puts(tb(results, headers = header))
        with indent(2):
            puts(colored.blue('\nTo download a genome and setup for use:'))
        with indent(4):
            puts(colored.green("\ntb.py genome --setup=<asm_name>\n"))
