from vcfkit import phylo
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing


def test_phylo_fasta():
    with Capturing() as out:
        phylo.main(["phylo", "fasta", "test_data/test.vcf.gz"])
    h = hashlib.sha224(str(out).strip()).hexdigest()
    assert h == "3e13b81be200e1128d6ee52ff68139d84ac7d638bc3bf3d272241766"


def test_phylo_nj():
    with Capturing() as out:
        phylo.main(["phylo", "tree", "nj", "test_data/test.vcf.gz"])
    h = hashlib.sha224(str(out).strip()).hexdigest()
    assert h == "5faa8909fab5d5141e8232d1fb05705a01cd608cba9610d740166fff"


def test_phylo_upgma():
    with Capturing() as out:
        comm = ["phylo", "tree", "upgma", "test_data/test.vcf.gz"]
        phylo.main(comm)
    h = hashlib.sha224(str(out).strip()).hexdigest()
    assert h == "4d0deed3dd8421ebf4c3f4458e027c7f2af7c4a699bfd8ce0f3067cf"


def test_phylo_plot():
    comm = ["phylo", "tree", "upgma", "--plot", "test_data/test.vcf.gz"]
    out = phylo.main(comm)
    h = hashlib.sha224(str(out).strip()).hexdigest()