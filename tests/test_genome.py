import pytest
from subprocess import check_output, Popen, PIPE, call
import hashlib
import os
from os.path import expanduser
from os import listdir
import re
from vcfkit import genome
from test import terminal


def test_genome_location():
    directory = genome.main(["genome","location"])
    assert os.path.exists(os.path.expanduser("~/.genome"))


def test_set_genome_location():
    directory = genome.main(["genome","location", "."])
    assert directory == os.getcwd()
    # Set genome location back to default
    genome.main(["genome","location", "-"])


def test_search_genome():
    results = genome.main(["genome", "--search", "coli"])
    search_result = [x for x in results if x[0] == "GCF_900042795.1"][0][3]
    assert search_result == "F1L3"


def test_download_genome():
    genome.main(["genome", "ncbi", "--ref", "F1L3"])
    genome_files = [x for x in listdir(expanduser("~/.genome/F1L3")) if x != ".DS_Store"]
    assert len(genome_files) == 11


def test_samtools_idx():
    out = check_output("samtools faidx ~/.genome/F1L3/F1L3.fa.gz NZ_FCPC01000094.1:1-20", shell=True)
    seq = out.splitlines()[1]
    assert seq == "CCTCACCGGATAACGCCGGC"


def test_amb_file():
    """
        The amb file records the appearance of non-ATGC characters (e.g. N) 
    """
    amb_file = expanduser("~/.genome/F1L3/") + "F1L3.fa.gz.amb"
    with open(amb_file, "r") as f:
        assert '5252759 94 1 575 10 N' == ' '.join(f.read().splitlines())


def test_ann_file():
    ann = expanduser("~/.genome/F1L3/") + "F1L3.fa.gz.ann"
    with open(ann, "r") as f:
        f = f.read()
        assert f.startswith("5252759 94 11")
        assert len(f.splitlines()) == 189


def test_fai_file():
    """
        Fasta Index
    """
    fai = expanduser("~/.genome/F1L3/") + "F1L3.fa.gz.fai"
    fai = open(fai, 'r').read().splitlines()
    assert len(fai) == 94
    assert fai[0] == "NZ_FCPC01000001.1\t436251\t19\t80\t81"

def test_list_genomes():
    out, err = terminal(["vk", "genome", "list"])
    assert "F1L3" in err

