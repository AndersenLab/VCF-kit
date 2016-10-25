import pytest
from subprocess import check_output, Popen, PIPE, call
import hashlib
import os
from os.path import expanduser
from os import listdir
import re
from vcfkit import genome


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


def test_download_genome_files():
    genome_files = [expanduser("~/.genome/F1L3/") + x for x in listdir(expanduser("~/.genome/F1L3"))]
    # Hash files
    hashed_files = set([hashlib.sha1(open(fname, 'rb').read()).hexdigest() for fname in genome_files if not fname.endswith(".nin")])
    assert '9568184f20821902ac0bbc1f4937691928515488' in hashed_files
    assert '6ccd1454ba12fef067d637bb14300841bde3ed81' in hashed_files
    assert '5c77732ce11c5cc16622de85d322e2bf157c3489' in hashed_files
    assert '5a67b3343c618392bc6e10ae9241961b43509b91' in hashed_files