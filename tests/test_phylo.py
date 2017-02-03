from vcfkit import phylo
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing


def test_phylo_fasta():
    with Capturing() as out:
        phylo.main(["phylo", "fasta", "test_data/test.vcf.gz"])
    h = hashlib.sha224(str(out)).hexdigest()
    assert h == "81eda4b80e96eed31749e0a28d1acbd345fcd216023d17735978e999"


def test_phylo_nj():
    with Capturing() as out:
        phylo.main(["phylo", "tree", "nj", "test_data/test.vcf.gz"])
    h = hashlib.sha224(str(out)).hexdigest()
    assert h == "96cdd2d7221bcbea8d86e099fdc1af02594d9dcad4daf8a5f85d0557"


def test_phylo_upgma():
    with Capturing() as out:
        comm = ["phylo", "tree", "upgma", "test_data/test.vcf.gz"]
        phylo.main(comm)
    h = hashlib.sha224(str(out)).hexdigest()
    assert h == "e3ef017c77bca5fd3feef7ed6215a37a5f0bb41101c125b3e493049e"


def test_phylo_plot():
    comm = ["phylo", "tree", "upgma", "--plot", "test_data/test.vcf.gz"]
    out = phylo.main(comm)
    h = hashlib.sha224(str(out)).hexdigest()