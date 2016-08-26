from vcfkit import phylo
from subprocess import Popen, PIPE
import hashlib
from test import Capturing

def test_phylo_fasta():
    with Capturing() as out:
        phylo.main(["phylo", "fasta", "data/test.vcf.gz"])
    h = hashlib.sha224(str(out)).hexdigest()
    assert h == "610ba5d2b533d730b6e9fb3005453d3d94710c9021fa5ca46fb5d7d8"


def test_phylo_nj():
    with Capturing() as out:
        phylo.main(["phylo", "tree", "nj", "data/test.vcf.gz"])
    h = hashlib.sha224(str(out)).hexdigest()
    assert h == "6c21e0508729c604d403a29cb40dca810455009fec793a9a99da9927"


def test_phylo_upgma():
    with Capturing() as out:
        comm = ["phylo", "tree", "upgma", "data/test.vcf.gz"]
        phylo.main(comm)
    h = hashlib.sha224(str(out)).hexdigest()
    assert h == "cd58f6912c40715ef5b2b1d8a381743454eb42cc6d0ba8a6c847ea83"
