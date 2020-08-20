# Overview


The `vcf2tsv` command can be used to convert a VCF to a tab-delimited file. 

Converts a VCF into a tsv - in wide or long format, and taking into account annotation fields (ANN) added by programs such as SNPeff.

```
usage:
  vk vcf2tsv (wide|long) [--print-header --ANN] <vcf>
```

* __wide | long__ - Select one of __wide__ or __long__ to set the output format. 
* __--print-header__ - Print a header row specifying column names.
* __--ANN__ - Parse annotation fields. Results in 1 row x (n)annotations x (n)samples when long or 1 row x (n)annotations x (n)variants when wide.

* `wide` - Prints a row for every CHROM:POS x # annotations (if `--ANN` is specified). A row will be output for every annotation.
* `long` - Prints a row for every CHROM:POS x sample x annotation (if `--ANN` is specified) combination.

For example, if a SNP at CHROM=II, POS=100 has two annotations and five samples, there will be (CHROM:POS[1]) x (annotations[2]) x (samples[5]) = 10 lines output.
* `--ANN` - Expands ANN fields added by variant prediction programs (_e.g._ SnpEff). When multiple annotations are present, a row is output for each annotation. This is for both wide and long format.

