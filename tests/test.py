import pytest
import vcfkit
from vcfkit import *
from subprocess import check_output, Popen, PIPE
import hashlib

from cStringIO import StringIO
import sys

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


class Capturing_err(list):
    def __enter__(self):
        self._stderr = sys.stderr
        sys.stderr = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stderr = self._stderr

def terminal(command):
    return Popen(command, stdout=PIPE, stderr=PIPE ).communicate()

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