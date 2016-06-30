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


F1L3_hashed = set(['1d24107a072a588de6b843036b4584938140a1c6', # F1L3.fa.gz.pac
                   '9568184f20821902ac0bbc1f4937691928515488', # F1L3.fa.gz.bwt
                   '85cef94245fdd399b692df36ba0ad7256f49904b', # F1L3.fa.gz.fai
                   'd1467e4a070af05eda87aa57890af5a81609f6b7', # F1L3.fa.gz.nsq
                   'bf0c4dac3b200d44cd567cdbab59a5a5cbdecb6a', # F1L3.fa.gz
                   '5a67b3343c618392bc6e10ae9241961b43509b91', # F1L3.fa.gz.amb
                   '5c77732ce11c5cc16622de85d322e2bf157c3489', # F1L3.fa.gz.nhr
                   '6ccd1454ba12fef067d637bb14300841bde3ed81', # F1L3.fa.gz.ann
                   'f68b57c5633a7fcb0e16f1a8f41c5927ef3f17aa', # F1L3.fa.gz.sa
                   # '8b5aeb9a99e902b2f34a8b2b7447af28b9d1fccc',
                   '29c1333a16e2f643048442e18da55eaf7df0931b']) # F1L3.fa.gz.gzi

def test_download_genome_files():
    genome_files = [expanduser("~/.genome/F1L3/") + x for x in listdir(expanduser("~/.genome/F1L3"))]
    # Hash files
    hashed_files = set([hashlib.sha1(open(fname, 'rb').read()).hexdigest() for fname in genome_files if not fname.endswith(".nin")])
    assert hashed_files == F1L3_hashed