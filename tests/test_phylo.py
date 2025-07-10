from vcfkit import phylo
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing


def test_phylo_fasta():
    with Capturing() as out:
        phylo.main(["phylo", "fasta", "test_data/test.vcf.gz"])
    last_line = out[-1]
    assert 2958 == len(last_line)
    assert last_line.startswith("AANNCTC")


def test_phylo_nj():
    with Capturing() as out:
        phylo.main(["phylo", "tree", "nj", "test_data/test.vcf.gz"])
    line = out[0].split(',')[11]
    assert line.startswith("DL200:0.1099")


def test_phylo_upgma():
    with Capturing() as out:
        phylo.main(["phylo", "tree", "upgma", "test_data/test.vcf.gz"])
    line = out[0].split(',')[11]
    assert line.startswith("DL200:0.1201")
