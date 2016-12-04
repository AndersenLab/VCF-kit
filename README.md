[![Build Status](https://travis-ci.org/AndersenLab/VCF-kit.svg?branch=master)](https://travis-ci.org/AndersenLab/VCF-kit) [![Coverage Status](https://coveralls.io/repos/github/AndersenLab/vcf-kit/badge.svg?branch=master)](https://coveralls.io/github/AndersenLab/vcf-kit?branch=master) [![Documentation Status](https://readthedocs.org/projects/vcf-kit/badge/?version=latest)](http://vcf-kit.readthedocs.io/en/latest/?badge=latest)
         

VCF-kit - [Documentation](http://vcf-kit.readthedocs.io/en/latest/?badge=latest)
===========

VCF-kit is a command-line based collection of utilities for performing analysis on Variant Call Format (VCF) files. A summary of the commands is provided below.

| Command |Description                                                                                 |
|:----------|------------------------------------------------------------------------------------------|
| __calc__ | Obtain  frequency/count of  genotypes and alleles.                                               |
| __call__ | Compare variants  identified  from  sequences obtained  through alternative methods against a VCF. |
| __filter__ | Filter  variants  with  a minimum or  maximum number  of  REF,  HET,  ALT,  or  missing calls.         |
| __geno__ | Various operations  at  the genotype  level.                                                      |
| __genome__ | Reference genome  processing  and management.                                                  |
| __hmm__ | Hidden-markov model for use in  imputing  genotypes from  parental  genotypes in  linkage studies.   |
| __phylo__ | Generate  dendrograms from  a VCF.                                                              |
| __primer__ | Generate  primers for variant validation.                                                     |
| __rename__ | Add a prefix, suffix, or  substitute  a string  in  sample  names.                                |
| __tajima__ | Calculate Tajima’s  D.                                                                        |
| __vcf2tsv__ | Convert a VCF to  TSV.                                                                       |

## Installation

```
pip install https://github.com/AndersenLab/vcf-kit/archive/0.0.2.tar.gz
```

__Installing Dependencies:__

In addition to python, VCF-kit requires that a number of additional programs be installed. `vk setup` can be used to install these dependencies. `vk setup` requires [homebrew](http://brew.sh/) (or if on linux, [linux-brew](http://linuxbrew.sh/)) to install programs used by VCF-kit. The programs are listed below followed by the versions they have been tested with.

* bwa (v 0.7.12)
* samtools (v 1.3)
* bcftools (v 1.3)
* blast (v 2.2.31+)
* muscle (v 3.8.31)