from vcfkit import freq
from subprocess import Popen, PIPE
import hashlib
from test import Capturing, terminal

def test_sample_hom_gt():
    with Capturing() as out:
        freq.main(["freq", "sample_hom_gt", "data/test.vcf.gz"])
    out = str(out).splitlines()
    assert out[0] == 'sample\tfreq_of_gt\tn_gt_at_freq'
    assert out[10] == 'QG536\t10\t11'
