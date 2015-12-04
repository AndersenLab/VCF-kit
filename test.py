import pytest
from subprocess import check_output, Popen, PIPE
import hashlib

def test_tajima():
    out, err = Popen(["tb","tajima","--no-header","--print-header","100000", "1000", "test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(out)
    assert m.hexdigest() == "d41d8cd98f00b204e9800998ecf8427e"


def test_vcf2tsv_long():
    out, err = Popen(["tb","vcf2tsv","long","--print-header","test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(out)
    assert m.hexdigest() == "5d36ceaf5298e6c63b76cd13eefac647"

def test_vcf2tsv_wide():
    out, err = Popen(["tb","vcf2tsv","wide","--print-header","test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(out)
    assert m.hexdigest() == "0fed86670ec7d2d2cae9b108b70d39d9"
