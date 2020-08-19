from vcfkit import calc
from subprocess import Popen, PIPE
import hashlib
from tests.test_utilities import Capturing

def test_sample_hom_gt():
    with Capturing() as out:
        calc.main(["calc", "sample_hom_gt", "test_data/test.vcf.gz"])
    assert out[0] == 'sample\tfreq_of_gt\tn_gt_at_freq'
    assert out[10] == 'AB1\t10\t8'


def test_spectrum_genotypes_count():
    with Capturing() as out:
        calc.main(["calc", "genotypes", "test_data/test.vcf.gz"])
    assert out[0] == 'n\tref\thet\talt\tmis'
    assert out[-1] == '1\t0\t0\t12\t2'


def test_spectrum_genotypes_frequency():
    with Capturing() as out:
        calc.main(["calc", "genotypes", "test_data/test.vcf.gz"])
    assert out[-1] == '1\t0\t0\t12\t2'


def test_spectrum_alleles_count():
    with Capturing() as out:
        calc.main(["calc", "spectrum", "test_data/test.vcf.gz"])
    assert round(float(out[-2].split("\t")[1]), 5) == 0.03571

