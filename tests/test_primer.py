from vcfkit import primer, genome
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing


def test_template():
    if not os.path.exists(os.path.expanduser("~/.genome/WBcel235/WBcel235.fa.gz")):
        genome.main(["genome","ncbi","--ref=WBcel235"])
    with Capturing() as out:
        primer.main(["primer", "template", "--ref=WBcel235","--region=I:1-100000", "data/test.vcf.gz"])
    
    t = eval(str(out))[1].split("\t")[6][0:20]
    assert t == "aacagtaaaaaaatcagtat"
