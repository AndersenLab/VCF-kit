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
    print out
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(str(out.splitlines()))
    assert m.hexdigest() == "0725e38970302c7389ca630d33cf42a1"

def test_vcf2tsv_wide():
    out, err = Popen(["tb","vcf2tsv","wide","--print-header","test.vcf.gz"], stdout = PIPE).communicate()
    if err:
        raise Exception(err)
    m = hashlib.md5()
    m.update(str(out.splitlines()))
    assert m.hexdigest() == "69b154692c6c85fcd30aa43ce5b5a00a"

#def test_hmm():
#    out, err = Popen(["tb", "hmm"])