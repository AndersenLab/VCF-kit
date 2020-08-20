from vcfkit import call, genome
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing
import os

# def test_call_DL238():
#     if not os.path.exists(os.path.expanduser("~/.genome/WBcel235/WBcel235.fa.gz")):
#         genome.main(["genome","ncbi","--ref=WBcel235"])
#     with Capturing() as out:
#         call.main(["call", "test_data/DL238.ab1", "--ref=WBcel235", "--vcf-sites", "test_data/DL238.vcf.gz"])
#     out = eval(str(out))
#     assert len(out) == 5
#     out = out[1].split("\t")
#     assert out[0] == "X"
#     assert int(out[1]) == 14557228
#     assert out[2] == out[3]
#     assert out[5] == out[6]
#     assert out[9] == "TN"
