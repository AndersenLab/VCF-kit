import pytest
from subprocess import check_output, Popen, PIPE, call
import hashlib
import os
from os.path import expanduser
from os import listdir
import re

def test_genome_location():
    out, err = Popen(["vk","genome","location"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    directory = out[out.find("Genome Directory"):].split(" ")[2].strip()
    assert os.path.exists(directory)


def test_set_genome_location():
    out, err = Popen(["vk","genome","location", "."], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    print(out[out.find("Set genome location"):].split(" ")[3])
    directory = out[out.find("Set genome location"):].split(" ")[4].strip()
    assert directory == os.getcwd() + "/"
    # Set genome location back to default
    out, err = Popen(["vk","genome","location", "-"], stdout = PIPE).communicate()


def test_search_genome():
    out, err = Popen(["vk", "genome", "--search", "coli"], stdout = PIPE).communicate()
    search_result = [re.split("\W+", x.strip()) for x in out.splitlines() if x.find("GCF_900042795.1") > 0][0][5]
    assert search_result == "F1L3"


def test_download_genome():
    call(["vk", "genome", "--download", "F1L3"])
    genome_files = listdir(expanduser("~/.genome/F1L3"))
    assert len(genome_files) == 11


def test_download_genome_files():
    genome_files = [expanduser("~/.genome/F1L3/") + x for x in listdir(expanduser("~/.genome/F1L3"))]
    # Hash files
    hashed_files = sorted([(os.path.split(fname)[1], hashlib.sha256(open(fname, 'rb').read()).digest()) for fname in genome_files])
    m = hashlib.sha256(str(hashed_files)).hexdigest()
    assert m == "ffef1336ec82212f2b2f10aa9bc5386af92646f53c576f0480270c092d770648"