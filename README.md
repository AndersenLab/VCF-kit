[![Build Status](https://travis-ci.org/AndersenLab/VCF-kit.svg?branch=master)](https://travis-ci.org/AndersenLab/VCF-kit) [![Coverage Status](https://coveralls.io/repos/github/AndersenLab/vcf-kit/badge.svg?branch=master)](https://coveralls.io/github/AndersenLab/vcf-kit?branch=master) [![Documentation Status](https://readthedocs.org/projects/vcf-kit/badge/?version=latest)](http://vcf-kit.readthedocs.io/en/latest/?badge=latest)
         

vcf-kit
===========

## Usage

```
usage:
  vk <command> [<args>...] 
  vk setup
  vk -h | --help
  vk --version

commands:
  call
  filter
  freq
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

