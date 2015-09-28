[![Build Status](https://travis-ci.org/AndersenLab/vcf-toolbox.svg?branch=v2)](https://travis-ci.org/AndersenLab/vcf-toolbox)

vcf-toolbox
===========

### Usage

	VCF-Toolbox 0.1

	usage:
	  tb.py <command> [<args>...]
	  tb.py -h | --help
	  tb.py --version

	commands:
	  tajima
	  primer
	  genome
	  phylo
	  freq

### Commands

#### tajima

#### primer

#### genome

#### phylo

#### freq

Calculates the frequency of homozygous genotypes by sample (e.g. number of singletons, doubletons, tripletons, etc. by sample)

```
	tb.py freq <vcf>
```

__Results__

In the table below - the first line indicates that the sample __ED3052__ has 2 singletons (private alleles). The second line indicates there are 3 doubltons.

| sample | freq_of_gt | n_gt_at_freq |
|--------|------------|--------------|
| ED3052 | 1          | 2            |
| ED3052 | 2          | 3            |
| ED3052 | 3          | 2            |
| ED3052 | 4          | 2            |
| ED3052 | 5          | 2            |
| ED3052 | 6          | 1            |
| ED3052 | 7          | 2            |
| ED3052 | 8          | 2            |
| ED3052 | 9          | 1            |

### Roadmap

#### Tajima-D

* [ ] Tajima-D
	* [ ] Documentation

#### Genotyping

* [ ] Primer Design 
	* [ ] snip-SNP
	* [ ] Indels
	* [ ] IDT - PCR order form

#### Genome Manager

* [ ] Check that tools (bwa/samtools/blast) are available.
	* [ ] Skip bwa if not available
* [ ] Error Checking