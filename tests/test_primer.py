from vcfkit import primer, genome
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing
import os

def test_template():
    if not os.path.exists(os.path.expanduser("~/.genome/WBcel235/WBcel235.fa.gz")):
        genome.main(["genome","ncbi","--ref=WBcel235"])
    with Capturing() as out:
        primer.main(["primer", "template", "--size=500", "--ref=WBcel235","--region=I:1-100000", "data/test.vcf.gz"])
    
    t = eval(str(out))[1].split("\t")[6][0:20]
    assert t == "ATTTCCCTCTCCGGGAAATG"


def test_snip():
    if not os.path.exists(os.path.expanduser("~/.genome/WBcel235/WBcel235.fa.gz")):
        genome.main(["genome","ncbi","--ref=WBcel235"])
    with Capturing() as out:
        primer.main(["primer", "snip", "--size=500", "--ref=WBcel235","--region=I:542216-543215", "data/test.vcf.gz"])
    
    t = eval(str(out))[1].split("\t")[7:9]
    assert t == ['496:496,184', '429,496:429,67,184']



def test_sanger():
    with Capturing() as out:
        primer.main(["primer", "sanger", "--ref=WBcel235","--region=I:542216-543215", "data/test.vcf.gz"])
    
    t = eval(str(out))[1].split("\t")[8:10]
    assert t == ['GAGAAGGACGGGACCCTTTG', 'aGGCCAGAACCTCGTGAAAC']
