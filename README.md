[![Build Status](https://travis-ci.org/AndersenLab/VCF-kit.svg?branch=master)](https://travis-ci.org/AndersenLab/VCF-kit) [![Coverage Status](https://coveralls.io/repos/github/AndersenLab/vcf-kit/badge.svg?branch=master)](https://coveralls.io/github/AndersenLab/vcf-kit?branch=master) [![Documentation Status](https://readthedocs.org/projects/vcf-kit/badge/?version=latest)](http://vcf-kit.readthedocs.io/en/latest/?badge=latest)
         

VCF-kit
===========

VCF-kit is a collection of utilities for performing analysis on variant caller format (VCF) files. A summary of the commands is provided below.

| Command |Description                                                                                 |
|:----------|------------------------------------------------------------------------------------------|
| calc | Obtain  frequency/count of  genotypes and alleles.                                               |
| call | Compare variants  identified  from  sequences obtained  through alternative methods against a VCF. |
| filter | Filter  variants  with  a minimum or  maximum number  of  REF,  HET,  ALT,  or  missing calls.         |
| geno | Various operations  at  the genotype  level.                                                      |
| genome | Reference genome  processing  and management.                                                  |
| hmm | Hidden-markov model for use in  imputing  genotypes from  parental  genotypes in  linkage studies.   |
| phylo | Generate  dendrograms from  a VCF.                                                              |
| primer | Generate  primers for variant validation.                                                     |
| rename | Add a prefix, suffix, or  substitute  a string  in  sample  names.                                |
| tajima | Calculate Tajima’s  D.                                                                        |
| vcf2sql | Generate  a tab-separated values  (TSV) and schema  for loading a VCF into  a database.         |
| vcf2tsv | Convert a VCF to  TSV.                                                                       |

## Usage

```
usage:
  vk <command> [<args>...] 
  vk setup
  vk -h | --help
  vk --version

commands:
  calc
  call
  filter
  geno
  genome
  hmm
  phylo
  primer
  rename
  tajima
  vcf2sql
  vcf2tsv
```

## Installation

```
pip install https://github.com/AndersenLab/vcf-kit/archive/0.0.2.tar.gz
```

__Installing Dependencies:__

```
vk setup
```

`vk setup` uses [homebrew](http://brew.sh/) (or if on linux, [linux-brew](http://linuxbrew.sh/)) to install programs used by vcf-kit. Versions listed have been tested:

* bwa (v 0.7.12)
* samtools (v 1.3)
* bcftools (v 1.3)
* blast (v 2.2.31+)
* muscle (v 3.8.31)

## To Do (Broad)

* [ ] Handle missing '-' input if stdin
* [ ] Error function 
* [ ] Replace `insert_header_line()` vcf function with cyvcf2 `add_info_to_header`
* [ ] Replace iterators with for line in vcf("I:1-1000"); replace `variant_line`

