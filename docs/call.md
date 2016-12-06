# Overview

The `call` command can be used to compare variants identified from Sanger sequencing with those present within a VCF.

### call

```
    vk call <seq> --ref=<reference> [--all-sites --vcf-targets <vcf>]
    vk call alignments <seq>  [--ref=<reference>]
```

Perform variant calling using blast. Useful for validating variants using sanger sequencing.  

* __<seq>__ - Fasta (fa), Fastq (fq, fastq), or ab1 format, determined by file extension containing sanger reads.
* __--ref__

* [ ] Add option to output INFO and FORMAT data for every variant.
