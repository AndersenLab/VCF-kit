import pytest
from subprocess import check_output


def test_vcf2tsv_long():
    out = check_output("tb vcf2tsv long --print-header test.vcf.gz | head -n 20 | md5", shell = True).strip()
    print(out)
    assert out == "d010d3707c2a1d8ab3613d35d5facc96"

def test_vcf2tsv_wide():
    out = check_output("tb vcf2tsv wide --print-header test.vcf.gz | head -n 20 | md5", shell = True).strip()
    print(out)
    assert out == "99081d87fa3af617693daa8d9d46b82f"
