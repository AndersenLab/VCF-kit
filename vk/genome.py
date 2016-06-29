#! /usr/bin/env python
"""
usage:
  vk genome location [<path>]
  vk genome list
  vk genome --search=<term>
  vk genome [options] --download=<asm_name> [--accession-chrom-names]
  vk genome [options] wormbase (--species)

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --directory=<dir>           Set Genome Directory

"""
from docopt import docopt
from utils import which
from utils.vcf import *
from utils.reference import *
import sys
from clint.textui import colored, puts, puts_err, indent, progress
import gzip
from subprocess import call
from time import time
from tabulate import tabulate as tab
import re
import requests
import os


# ftplib
debug = None
if len(sys.argv) == 1:
    debug = ['genome', "test.vcf.gz"]

if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)

    # Setup Genomes Directory
    if args["location"] and args["<path>"]:
        if args["<path>"] == "-":
            genome_directory = get_genome_directory_file()
            os.remove(genome_directory)
            get_genome_directory_file()
        else:
            with open(get_genome_directory_file(), "w") as f:
                genome_directory=os.path.realpath(args["<path>"])
                with indent(2):
                    puts(colored.blue("\nSet genome location to: " + genome_directory + "/\n"))
                f.write(genome_directory)
                # create directory if not exists
                if not os.path.exists(genome_directory):
                    os.makedirs(genome_directory)
                exit()

    if args["--directory"]:
        genome_directory = os.path.realpath(args["--download"])
    else:
        genome_directory = get_genome_directory()

    with indent(2):
        puts(colored.blue("\nGenome Directory: " + genome_directory + "\n"))
    
    genome_db = get_genome_directory() + "/genomes.db"

    ################
    # List Genomes #
    ################
    if args["list"]:
        genome_list = get_genome_list()
        with indent(2):
            exit(puts(colored.blue("\n".join(genome_list) + "\n")))

    ##################
    # Search Genomes #
    ##################
    elif args["--search"]:
        # Download and cache a list of genomes from NCBI for searching
        if os.path.isfile(genome_db):
            fileTime = os.path.getctime(genome_db)
        else:
            fileTime = 0
        if fileTime < (time() - 60 * 60 * 24 * 2):
            with indent(2):
                puts(colored.blue('\nDownloading list of reference genomes\n'))
            r = requests.get("http://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt")
            genome_file = open(genome_db, "w")
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

        with open(genome_db, "r") as f:
            results = []
            for line in f:
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    line = [x for k, x in enumerate(line) if k in [0, 1, 7, 15, 19]]
                    if args["--search"].lower() in line[2].lower() and line[4] != "na":
                        results.append(line)
        with indent(4):
            puts(tab(results, headers = header))
        with indent(2):
            puts(colored.blue('\nTo download a genome and setup for use:'))
        with indent(4):
            puts(colored.green("\nvk genome --download=<asm_name>\n"))
    elif args["--download"]:
        with open(genome_db, "r") as f:
            results = []
            for line in f:
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    line = [x for k, x in enumerate(line) if k in [0,1,7,15,19]]
                    if args["--download"] == line[3]:
                        results.append(line)
        if len(results) == 0:
            with indent(2):
                puts(colored.red('\nError: Genome not found\n'))
        else:
            reference_download = results[0]
            if not args["--directory"]:
                reference_directory = genome_directory + "/" + args["--download"] + "/"
            else:
                reference_directory = genome_directory + "/"
            if not os.path.exists(reference_directory):
                os.makedirs(reference_directory)
            url = reference_download[4].replace("ftp://", "http://") + "/" + os.path.split(reference_download[4])[1] + "_genomic.fna.gz"
            if not os.path.exists(reference_directory):
                os.makedirs(reference_directory)
            with indent(2):
                puts(colored.green('\nDownloading: ' + reference_download[3] + "; " + url + '\n'))

            # stack overflow: 15644964; 
            r = requests.get(url, stream = True)
            ref_filename = reference_directory + "/" + reference_download[3] + ".tmp.fa.gz"
            with open(ref_filename , 'wb') as f:
                total_length = int(r.headers.get('content-length'))
                for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length/1024) + 1): 
                    if chunk:
                        f.write(chunk)
                        f.flush()


            # Fix chromosome names
            if not args["--accession-chrom-names"]:
                with indent(2):
                    puts(colored.green('\nFixing chromosome names\n'))

                with open(ref_filename.replace(".fa.gz",".fa"), 'w') as outfa:
                    with gzip.open(ref_filename, 'rb') as f:
                        for line in f:
                            outline = line
                            if line.startswith(">"):
                                chrom_name = re.match(".*[c|C]hromosome ([A-Za-z0-9]+)[W]*", line)
                                if chrom_name is not None:
                                    outline = ">" + chrom_name.group(1) + "\n"
                                elif line.lower().find("mitochon") > 0:
                                    outline = ">MtDNA\n"
                                puts(colored.blue(line.strip("\n>")) + " --> " + colored.blue(outline.strip("\n>")))
                            outfa.write(outline)

            if which("bgzip"):
                with indent(2):
                    puts(colored.green('\nSwitching from gzip to bgzip\n'))
                # Convert to bgzip
                if args["--accession-chrom-names"]:
                    call(["gunzip", "-f", ref_filename])
                comm_bgzip = "bgzip -fc {ref_filename} > {ref_out}"
                comm_bgzip = comm_bgzip.format(ref_filename = ref_filename.replace(".fa.gz",".fa"),
                                  ref_out = ref_filename.replace(".tmp",""))
                call(comm_bgzip, shell = True)
                ref_filename = ref_filename.replace(".tmp","")
            else:
                puts_err(colored.red("Please install bgzip."))
                exit()

            if which("bwa"):
                with indent(2):
                    puts(colored.green("\nCreating bwa index\n"))
                call(["bwa", "index", ref_filename])

            if which("samtools"):
                with indent(2):
                    puts(colored.green("\nCreating samtools index\n"))
                call(["samtools", "faidx", ref_filename])

            if which("makeblastdb"):
                with indent(2):
                    puts(colored.green("\nCreating blast index\n"))
                comm = "gunzip -c {ref} | makeblastdb -in - -dbtype=nucl -title={ref} -out={ref}".format(ref=ref_filename)
                call(comm, shell = True)

            # Remove temp files
            if args["--accession-chrom-names"]:
                os.remove(ref_filename.replace(".fa.gz",".tmp.fa.gz"))

            # Remove temporary files
            os.remove(ref_filename.replace(".fa.gz",".tmp.fa.gz"))
            os.remove(ref_filename.replace(".fa.gz",".tmp.fa"))

            # Add error checking here...

            with indent(2):
                puts(colored.green("\nComplete!\n"))

