# Overview

The `vcf2tsv` command can be used to convert a VCF to a tab-delimited file. 

```
usage:
  vk vcf2tsv (wide|long) [--print-header --ANN] <vcf>
```

# Options

* `wide` - Prints a row for every CHROM:POS x # annotations (if `--ANN` is specified). A row will be output for every annotation.
* `long` - Prints a row for every CHROM:POS x sample x annotation (if `--ANN` is specified) combination.

For example, if a SNP at CHROM=II, POS=100 has two annotations and five samples, there will be (CHROM:POS[1]) x (annotations[2]) x (samples[5]) = 10 lines output.

* `--print-header` - Prints a header row.
* `--ANN` - Expands ANN fields added by variant prediction programs (_e.g._ SnpEff). When multiple annotations are present, a row is output for each annotation. This is for both wide and long format.

