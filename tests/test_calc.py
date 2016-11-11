from vcfkit import calc
from subprocess import Popen, PIPE
import hashlib
from test import Capturing, terminal

def test_sample_hom_gt():
    with Capturing() as out:
        calc.main(["calc", "sample_hom_gt", "data/test.vcf.gz"])
    assert out[0] == 'sample\tfreq_of_gt\tn_gt_at_freq'
    assert out[10] == 'QG536\t10\t11'


def test_spectrum_genotypes_count():
    with Capturing() as out:
        calc.main(["calc", "genotypes", "data/test.vcf.gz"])
    assert out[0] == 'n\tref\thet\talt\tmis'
    assert out[-1] == '1\t9\t2\t1\t2'


def test_spectrum_genotypes_frequency():
    with Capturing() as out:
        calc.main(["calc", "genotypes", "data/test.vcf.gz"])
    assert out[-1] == '1\t9\t2\t1\t2'


def test_spectrum_alleles_count():
    with Capturing() as out:
        calc.main(["calc", "spectrum", "data/test.vcf.gz"])
    assert out[-2] == '23\t0.0357142857143'

