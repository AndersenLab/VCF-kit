from vcfkit import tajima
from subprocess import Popen, PIPE
from tests.test_utilities import Capturing
from vcfkit.utils import autoconvert

def nearly_equal(a,b,sig_fig=5):
    return ( a==b or 
             int(a*10**sig_fig) == int(b*10**sig_fig)
           )


def test_tajima():
    with Capturing() as out:
        tajima.main(["tajima","--no-header","100000", "1000", "data/test.vcf.gz"])
    assert out[0] == 'I\t6000\t106000\t2\t1\t-0.740993868265'
    assert out[-1] == 'X\t17567000\t17667000\t2\t1\t-0.740993868265'


def test_compare_vcftools():
    """
        Compare the Tajima's D Calculation with vcftools output.
    """
    comm = ["vcftools", "--TajimaD", "100000", "--gzvcf", "data/test.vcf.gz"]
    out, err = Popen(comm, stdout=PIPE, stderr=PIPE).communicate()
    vcftools_tajima = open("out.Tajima.D", 'r').read()
    vcftools_tajima = [x.split("\t")[0:4] for x in open("data/out.Tajima.D", 'r').read().splitlines()]
    with Capturing() as out:
        tajima.main(["tajima","--no-header","100000", "100000", "data/test.vcf.gz"])
    vcftools_tajima = dict([(x[0] + ":" + x[1], autoconvert(x[3])) for x in vcftools_tajima if x[3] != "nan"])
    out = [x.split("\t") for x in out]
    vcfkit_tajima = dict([(x[0] + ":" + x[1], autoconvert(x[5])) for x in out])
    for k in vcfkit_tajima.keys():
        assert nearly_equal(vcftools_tajima[k], vcfkit_tajima[k], sig_fig = 3) is True


