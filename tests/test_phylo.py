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
    h = hashlib.sha224(''.join(out)).hexdigest()
    assert h == "707c94f6edc7e5806f80cdb8ff02367cedb94d675ef1aafde19006ed"


def test_phylo_upgma():
    with Capturing() as out:
        comm = ["phylo", "tree", "upgma", "test_data/test.vcf.gz"]
        phylo.main(comm)
    h = hashlib.sha224(''.join(out)).hexdigest()
    assert h == "142ece5d98398583dccb9ca4ecf76cc0dc5689a6beaa54721a72a785"
