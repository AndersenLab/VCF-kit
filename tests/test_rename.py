from vcfkit import rename
from cyvcf2 import VCF
from subprocess import Popen, PIPE
import hashlib
from test_utilities import terminal
from cStringIO import StringIO


def test_rename_prefix():
    comm = "vk rename --prefix prefix data/test.vcf.gz | bcftools query --list-samples"
    out, err = Popen(comm, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    assert all([x.startswith("prefix") for x in out.splitlines()])


def test_rename_suffix():
    comm = "vk rename --suffix data/test.vcf.gz | bcftools query --list-samples"
    out, err = Popen(comm, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    assert all([x.endswith("suffix") for x in out.splitlines()])


def test_rename_subst():
    comm = "vk rename --subst N2:TEST --subst WN2001,TEST2 --subst AB1=TEST3 data/test.vcf.gz | bcftools query --list-samples"
    out, err = Popen(comm, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    assert sum([x.startswith("TEST") for x in out.splitlines()]) == 3


def test_rename_combo():    
    comm = "vk rename --prefix PREFIX_ --suffix _SUFFIX --subst N2:TEST --subst WN2001,TEST2 --subst AB1=TEST3 data/test.vcf.gz | bcftools query --list-samples"
    out, err = Popen(comm, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    assert all([x.startswith("PREFIX_") for x in out.splitlines()])
    assert sum(["TEST" in x for x in out.splitlines()]) == 3
    assert all([x.endswith("_SUFFIX") for x in out.splitlines()])
