#! /usr/bin/env python
"""
usage:
  tb.py genome <vcf>
  tb.py genome --search=<term>
  tb.py genome --download=<asm_name>

options:
  -h --help                   Show this screen.
  --version                   Show version.

"""
from docopt import docopt
from utils.vcf import *
import sys
from clint.textui import colored, puts, indent, progress
from Bio.Entrez import esearch, read, efetch, esummary
from subprocess import call
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
    module_path = os.path.split(os.path.realpath(__file__))[0]
    genome_file_path = module_path + "/genomes.txt"
    print args
    if args["<vcf>"] is not None:
        v = vcf(args["<vcf>"])
        with indent(2):
            puts(colored.blue('\nReference:\n'))
        with indent(4):
            puts(colored.green(v.metadata["reference"] + "\n"))

    elif args["--search"]:
        # Download and cache a list of genomes from NCBI for searching
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
    elif args["--download"]:
        with open(genome_file_path, "r") as f:
            results = []
            for line in f:
                line = line.strip().split("\t")
                line = [x for k,x in enumerate(line) if k in [0,1,7,15,19]]
                if args["--download"] == line[3]:
                    results.append(line)
        if len(results) == 0:
            with indent(2):
                puts(colored.red('\nError: Genome not found\n'))
        else:
            ref_dl = results[0]
            ref_dir = module_path + "/reference"
            url = ref_dl[4].replace("ftp://", "http://") + "/" + os.path.split(ref_dl[4])[1] + "_genomic.fna.gz"
            if not os.path.exists(ref_dir):
                os.makedirs(ref_dir)
            with indent(2):
                puts(colored.green('\nDownloading: ' + ref_dl[3] + "; " + url + '\n'))

            # stack overflow: 15644964; 
            
            r = requests.get(url, stream = True)
            ref_filename = ref_dir + "/" + ref_dl[3] + ".tmp.fa.gz"
            with open(ref_filename , 'wb') as f:
                total_length = int(r.headers.get('content-length'))
                for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length/1024) + 1): 
                    if chunk:
                        f.write(chunk)
                        f.flush()

            with indent(2):
                puts(colored.green('\nSwitching from gzip to bgzip\n'))
            # Convert to bgzip
            call(["gunzip", "-f", ref_filename])
            comm_bgzip = "bgzip --stdout {ref_filename} > {ref_out}"
            comm_bgzip = comm_bgzip.format(ref_filename = ref_filename.replace(".fa.gz",".fa"),
                              ref_out = ref_filename.replace(".tmp",""))
            call(comm_bgzip, shell = True)
            ref_filename = ref_filename.replace(".tmp","")

            with indent(2):
                puts(colored.green("\nCreating bwa index\n"))

            call(["bwa", "index", ref_filename])

            with indent(2):
                puts(colored.green("\nCreating samtools index\n"))

            call(["samtools", "faidx", ref_filename])

            with indent(2):
                puts(colored.green("\nCreating blast index\n"))

            call(["makeblastdb", "-in", ref_filename, "-dbtype=nucl"])

            # Add error checking here...

            with indent(2):
                puts(colored.green("\nComplete!\n"))

