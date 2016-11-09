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
    hashed_files = set([hashlib.sha256(open(fname, 'rb').read()).hexdigest() for fname in genome_files if not fname.endswith(".nin")])
    assert '96fcce5dc0e9a382d2ee5fd579255855656954f7d76859526f60d5e55a528a87' in hashed_files # F1L3.fa.gz
    assert '50aa468e0f13ba82aa87a364f7966e50315bf1bd71cccec5a69f1885311120d5' in hashed_files # F1L3.fa.gz.amb
    assert 'e98601ee2bc72b633527124c1142cd6211f28765159d592bcfef6b80ece3f390' in hashed_files # F1L3.fa.gz.ann
    assert 'e19e692153c4acfeaaf60e4b6d2c4adf3dc52954080b6195793dc2456b0aac5c' in hashed_files # F1L3.fa.gz.bwt
    assert 'efee15166ca3ba67cc1ae64b0b7425e1822bbc2d44a9f8fbf7c8e1062c3c423d' in hashed_files # F1L3.fa.gz.fai
    assert '6f67f186d73837024ff3b7dcd6799301497a713709d1e93d28231bf223132f9e' in hashed_files # F1L3.fa.gz.gzi
    assert '5fe458b1463bb36536848ac2c72f1f1b2b1dba8e71696a5e1da3a154ff1b9fcc' in hashed_files # F1L3.fa.gz.nhr
    assert '3f7f8e0a72338a7b767c6978bd7e73b5cd94171a5cbd1187363807060dde34a2' in hashed_files # F1L3.fa.gz.nsq
    assert '7f474886d402c7be74650370d6bc49db011d9d139bbfb5bf5cb42f1f65fe0ab7' in hashed_files # F1L3.fa.gz.pac
    assert '6d80cdd3d919c93a5a832d501c1ff63c0918bb24ae65bc9b72663ccf00b14d44' in hashed_files # F1L3.fa.gz.pac
        

