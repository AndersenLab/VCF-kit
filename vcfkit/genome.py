#! /usr/bin/env python
"""
usage:
  vk genome location [<path>]
  vk genome list
  vk genome --search=<term>
  vk genome ncbi [options] --ref=<asm_name> [--accession-chrom-names]
  vk genome wormbase [options] --ref=<asm_name>

options:
  -h --help                   Show this screen.
  --directory=<dir>           Set Genome Directory

"""
from docopt import docopt
from utils import which
from utils.vcf import *
from utils.reference import *
from clint.textui import colored, puts, puts_err, indent, progress
import gzip
from subprocess import call
from time import time
from tabulate import tabulate as tab
import requests
import os
import urllib
from sys import exit  # Used by exit(); don't remove.
from Bio import Entrez


def fetch_chrom_name(id):
    try:
        if not id.startswith("NC_"):
            return id
        Entrez.email = "vcf-kit@vcf-kit.com"
        chrom = Entrez.read(Entrez.efetch(db="nuccore", id=id, rettype="gb", retmode="xml"))
        gb_feature_quals = chrom[0]["GBSeq_feature-table"][0]["GBFeature_quals"]
        features = dict([x.values() for x in gb_feature_quals])
        if "organelle" in features:
            if features["organelle"] == "mitochondrion":
                return "MtDNA"
        else:
            chrom_name = features["chromosome"]
            return chrom_name
    except:
        return id


def download_genomes(genome_db):
    if os.path.isfile(genome_db):
        fileTime = os.path.getmtime(genome_db)
    else:
        fileTime = 0
    if (time() - fileTime) > (3 * 30 * 24 * 60 * 60) or is_non_zero_file(genome_db) is False:
        with indent(2):
            puts(colored.blue('\nDownloading list of reference genomes\n'))
        r = requests.get("http://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt")
        genome_file = open(genome_db, "w")
        with genome_file as f:
            f.write(r.text.encode('utf-8').strip())


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def main(debug=None):
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  argv=debug,
                  options_first=False)
    # Setup Genomes Directory
    if args["location"] and args["<path>"]:
        if args["<path>"] == "-":
            genome_directory = get_genome_directory_file()
            os.remove(genome_directory)
            return get_genome_directory_file()
        else:
            with open(get_genome_directory_file(), "w") as f:
                genome_directory = os.path.realpath(args["<path>"])
                with indent(2):
                    puts(colored.blue("\nSet genome location to: " + genome_directory + "/\n"))
                f.write(genome_directory)
                # create directory if not exists
                if not os.path.exists(genome_directory):
                    os.makedirs(genome_directory)
                return genome_directory

    if args["--directory"]:
        genome_directory = os.path.realpath(args["--ref"])
    else:
        genome_directory = get_genome_directory()

    with indent(2):
        puts(genome_directory)

    if args["location"] and not args["<path>"]:
        return genome_directory

    genome_db = get_genome_directory() + "/genomes.db"

    ################
    # List Genomes #
    ################
    if args["list"]:
        output_genome_list()

    ##################
    # Search Genomes #
    ##################
    elif args["--search"]:
        # Download and cache a list of genomes from NCBI for searching
        download_genomes(genome_db)

        # Cache result
        header = ["assembly_accession",  # 0
                  "bioproject",  # 1
                  "organism_name",  # 7
                  "asm_name",  # 15
                  "ftp_path"]  # 19

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
            puts(tab(results, headers=header))
        with indent(2):
            puts(colored.blue('\nTo download a genome and setup for use:'))
        with indent(4):
            puts(colored.green("\nvk genome ncbi --ref=<asm_name>\n"))
        return results
    elif args["--ref"]:
        # reference name.
        reference_name = args["--ref"]

        # Ensure genome db is available
        download_genomes(genome_db)

        # reference directory
        if not args["--directory"]:
            reference_directory = genome_directory + "/" + reference_name + "/"
        else:
            reference_directory = genome_directory + "/"
        if not os.path.exists(reference_directory):
            os.makedirs(reference_directory)

        # base reference filename.
        ref_filename = reference_directory + reference_name + ".tmp.fa.gz"

        if args["wormbase"]:
            asm_url = "ftp://ftp.wormbase.org/pub/wormbase/releases/{asm_name}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.{asm_name}.genomic.fa.gz"
            reference_download = asm_url.format(asm_name=args["--ref"])
            comm = "curl {reference_download} > {ref_filename}".format(**locals())
            print(comm)
            call(comm, shell = True)
            # Unzip wormbase genome
            call(["gunzip", "-f", ref_filename])
        else:
            # NCBI
            with open(genome_db, "r") as f:
                results = []
                for line in f:
                    if not line.startswith("#"):
                        line = line.strip().split("\t")
                        line = [x for k, x in enumerate(line) if k in [0, 1, 7, 15, 19]]
                        if args["--ref"] == line[3]:
                            results.append(line)
                            reference_download = results[0]
                            url = reference_download[4].replace("ftp://", "http://") + "/" + os.path.split(reference_download[4])[1] + "_genomic.fna.gz"
            if len(results) == 0:
                with indent(2):
                    puts(colored.red('\nError: Genome ' + args["--ref"] + ' not found\n'))

            with indent(2):
                puts(colored.green('\nDownloading: ' + reference_name + "; " + url + '\n'))

            # stack overflow: 15644964;
            r = requests.get(url, stream=True)

            with open(ref_filename, 'wb') as f:
                total_length = int(r.headers.get('content-length'))
                for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length / 1024) + 1):
                    if chunk:
                        f.write(chunk)
                        f.flush()

        # Fix chromosome names
        if not args["--accession-chrom-names"] and not args['wormbase']:
            with indent(2):
                puts(colored.green('\nFixing chromosome names\n'))

            with open(ref_filename.replace(".fa.gz", ".fa"), 'w') as outfa:
                with gzip.open(ref_filename, 'rb') as f:
                    for line in f:
                        outline = line
                        if line.startswith(">"):
                            acc = line.split(" ")[0].strip(">")
                            chrom_name = fetch_chrom_name(acc)
                            if chrom_name is not None:
                                outline = ">" + chrom_name + "\n"
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
            comm_bgzip = comm_bgzip.format(ref_filename=ref_filename.replace(".fa.gz", ".fa"),
                                           ref_out=ref_filename.replace(".tmp", ""))
            print(comm_bgzip)
            call(comm_bgzip, shell=True)
            ref_filename = ref_filename.replace(".tmp", "")
        else:
            with indent(2):
                puts_err(colored.red("Please install bgzip."))
            exit()

        if which("bwa"):
            with indent(2):
                puts(colored.green("\nCreating bwa index\n"))
            call(["bwa", "index", ref_filename])
        else:
            with indent(2):
                puts(colored.blue("\nSkipping bwa index; bwa not installed\n"))

        if which("samtools"):
            with indent(2):
                puts(colored.green("\nCreating samtools index\n"))
            call(["samtools", "faidx", ref_filename])
        else:
            with indent(2):
                puts(colored.blue("\nSkipping samtools index; Samtools not installed\n"))

        if which("makeblastdb"):
            with indent(2):
                puts(colored.green("\nCreating blast database\n"))
            comm = "gunzip -c {ref} | makeblastdb -in - -dbtype=nucl -title={ref} -out={ref}".format(ref=ref_filename)
            call(comm, shell=True)
        else:
            with indent(2):
                puts(colored.blue("\nSkipping creation of blast database; blast is not installed\n"))

        # Remove temp files
        if args["--accession-chrom-names"]:
            os.remove(ref_filename.replace(".fa.gz", ".tmp.fa.gz"))

        # Remove temporary files
        try:
            os.remove(ref_filename.replace(".fa.gz", ".tmp.fa.gz"))
            os.remove(ref_filename.replace(".fa.gz", ".tmp.fa"))
        except:
            pass

        # Add error checking here...

        with indent(2):
            puts(colored.green("\nComplete!\n"))


if __name__ == '__main__':
    main()
