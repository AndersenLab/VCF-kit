<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [vcf-kit](#vcf-kit)
  - [Usage](#usage)
  - [Commands](#commands)
    - [call](#call)
    - [tajima](#tajima)
          - [Calculate Tajima's D using a sliding window](#calculate-tajimas-d-using-a-sliding-window)
          - [Calculate Tajima's D using a continous sliding window](#calculate-tajimas-d-using-a-continous-sliding-window)
          - [Calculate Tajima's D using a bins](#calculate-tajimas-d-using-a-bins)
    - [primer](#primer)
    - [genome](#genome)
    - [hmm](#hmm)
    - [phylo](#phylo)
          - [Generate a fasta-alignment for variant calls](#generate-a-fasta-alignment-for-variant-calls)
          - [Generate a phylogenetic tree (newick format)](#generate-a-phylogenetic-tree-newick-format)
          - [Plot a phylogeny from a VCF file](#plot-a-phylogeny-from-a-vcf-file)
    - [geno](#geno)
          - [transfer-filter](#transfer-filter)
          - [het-polarization](#het-polarization)
    - [freq](#freq)
    - [vcf2tsv](#vcf2tsv)
    - [vcfcompare](#vcfcompare)
    - [vcf2sql](#vcf2sql)
    - [vcf2bigquery](#vcf2bigquery)
    - [google datastore / amazon](#google-datastore--amazon)
    - [Additional Features](#additional-features)
      - [Possible additions?](#possible-additions)
  - [To Do (Broad)](#to-do-broad)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[![Build Status](https://travis-ci.org/AndersenLab/vcf-kit.svg?branch=master)](https://travis-ci.org/AndersenLab/vcf-kit)

[![Coverage Status](https://coveralls.io/repos/github/AndersenLab/vcf-kit/badge.svg?branch=master)](https://coveralls.io/github/AndersenLab/vcf-kit?branch=master)

[![Documentation Status](https://readthedocs.org/projects/vcf-kit/badge/?version=latest)](http://vcf-kit.readthedocs.io/en/latest/?badge=latest)
         

vcf-kit
===========

## Usage

```
  VCF-Kit 0.1

  usage:
    vk <command> [<args>...]
    vk -h | --help
    vk --version

  commands:
  tajima
  hmm
  filter
  call
  primer
  genome
  rename
  phylo
  freq
  geno
  vcf2tsv
  vcf2sql
```

## Installation

```
pip install https://github.com/AndersenLab/vcf-kit/archive/0.0.2.tar.gz
```

__Installing Dependencies:__

```
vk setup
```

`vk setup` uses [homebrew](http://brew.sh/) (or if on linux, [linux-brew](http://linuxbrew.sh/)) to install programs used by vcf-kit:

* bwa
* samtools
* bcftools
* blast
* muscle

## Disclaimer

This isn't done yet.

## Commands

### call

```
    vk call <seq> --ref=<reference> [--all-sites --vcf-targets <vcf>]
    vk call alignments <seq>  [--ref=<reference>]
```

Perform variant calling using blast. Useful for validating variants using sanger sequencing.  

* __<seq>__ - Fasta (fa), Fastq (fq, fastq), or ab1 format, determined by file extension containing sanger reads.
* __--ref__

* [ ] Add option to output INFO and FORMAT data for every variant.


### primer

* [ ] Spike in primers with variants.

Suite of tools for genotyping: via sanger sequencing, using snip-SNPs, and indels. Generates appropriate primers and predicts band sizes for indels and snip-SNPs.

* [ ] Primer Design 
	* [ ] snip-SNP
	* [ ] Indels
  * [ ] Use bcftools consensus

###### transfer-filter

```
vk geno transfer-filter <vcf>
```

Generate a FORMAT field filter column (FT) and transfer filters applied within the FILTER column to the FT field. This utility is useful if variants are called individually (producing a single VCF) for samples and then later merged into a multi-sample VCF. It enables you to track which variant filters were applied to a specific sample rather than across a population. The following header line is added:

`##FORMAT=<ID=GF,Number=1,Type=String,Description="Genotype Filter">`


###### het-polarization

```
vk geno het-polarization <vcf>
```

Creates a new FORMAT field (HP) and "polarizes" or switches heterozygous genotypes based on genotype likelyhoods (GL) or Phred-scaled genotype likelihoods (PL). For example, a variant with 7 reads supporting a reference call and 1 read supporting an alternate allele might be called as a heterozygous genotype. However, in hermaphroditic species (_e.g. C. elegans_) where limited or no heterozygosity exists, it is much more likely that the genotype is homozygous reference. The HP FORMAT field will list one of the following depending on the change made:

* __AA__ - A heterozygous genotype was switched to a homozygous reference genotype.
* __AB__ - A heterozygous genotype remains a heterozygous genotype.
* __BB__ - A heterozygous genotype was switched to a homozygous alternate genotype.

The following header line is added:

`##FORMAT=<ID=HP,Number=1,Type=String,Description="Flag used to mark whether a variant was polarized">`


### freq

Calculates the frequency of homozygous genotypes by sample (e.g. number of singletons, doubletons, tripletons, etc. by sample)

```
vk freq <vcf>
```

Output from this utility appears as the table below. The first line indicates that the sample __ED3052__ has 2 singletons (private alleles). The second line indicates there are 3 doubltons.

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

### vcf2tsv

Converts a VCF into a tsv - in wide or long format, and taking into account annotation fields (ANN) added by programs such as SNPeff.

```
vk vcf2tsv (wide|long) [--print-header --snpeff] <vcf>
```

* __wide | long__ - Select one of __wide__ or __long__ to set the output format. 
* __--print-header__ - Print a header row specifying column names.
* __--ANN__ - Parse annotation fields. Results in 1 row x (n)annotations x (n)samples when long or 1 row x (n)annotations x (n)variants when wide.

* [X] Wide
* [X] Long
* [X] Parse ANN Fields (e.g. snpeff)
* [ ] Read from stdin

### vcfcompare

Tool for comparing concordance across genome.

### vcf2sql

Prepares a VCF for import into an sql database. Bigquery, Postgres, Mysql, and sqlite will be supported. Data is loaded denormalized and
indices are added to enable easy querying.

```
  vk vcf2sql (bigquery|postgres|mysql|sqlite) <vcf>
```

* [X] Bigquery
* [X] Mysql
* [X] Postgres
* [X] Sqlite
* [X] Automatically load
* [X] Support for SNPeff annotations
* [ ] Reorder columns
* [ ] Support for multi-column types
* [ ] Break load job for bigquery into multiple files if filesize > 4GB
* [ ] Output TGT (bases), allele 1, allele 2 for genotypes.

### vcf2bigquery

* [ ] Split bigquery into its own tool

### google datastore / amazon 

* [ ] Split out into separate tools.

### Additional Features

* [X] Check for BWA, Blast, Primer3, other req'd CLI tools.
* [ ] tstv
* [ ] GC

#### Possible additions?

* [ ] Upload to UCSC
* [ ] Annotate variants (track-based search)


## To Do (Broad)

* [ ] Handle missing '-' input if stdin
* [ ] Error function 
* [ ] Replace `insert_header_line()` vcf function with cyvcf2 `add_info_to_header`
* [ ] Replace iterators with for line in vcf("I:1-1000"); replace `variant_line`

