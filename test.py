import pytest
from subprocess import check_output, Popen, PIPE
import hashlib

def test_tajima():
    out, err = Popen(["tb","tajima","--no-header","--print-header","100000", "1000", "test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(str(out.splitlines()))
    assert m.hexdigest() == "d751713988987e9331980363e24189ce"


def test_vcf2tsv_long():
    out, err = Popen(["tb","vcf2tsv","long","--print-header","test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(str(out.splitlines()))
    assert m.hexdigest() == "38ec18fad37b9f43e0d5b9e5f6276e79"

def test_vcf2tsv_wide():
    out, err = Popen(["tb","vcf2tsv","wide","--print-header","test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(str(out.splitlines()))
    assert m.hexdigest() == "eca9efc25822aab59dfe29c42c1b01fc"
