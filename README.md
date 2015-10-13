[![Build Status](https://travis-ci.org/AndersenLab/vcf-toolbox.svg?branch=master)](https://travis-ci.org/AndersenLab/vcf-toolbox)

vcf-toolbox
===========

## Usage

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
	  geno
	  vcf2tsv

## Commands

## tajima

Generate a Tajima's D statistic using a sliding window or across bins. 

__Calculate Tajima's D using a sliding window__

```
	tb.py tajima 100000 1000 <vcf>
```

When run, the code above would calculate Tajima's D across a 100,000 bp sliding window that moves 1,000 bp with each iteratino.

__Calculate Tajima's D using a continous sliding window__
```
	tb.py tajima 100000 --sliding <vcf>
```

When run, the code above would calculate Tajima's D across a 100,000 bp sliding window that captures every unique bin of variants that fall within 100,000 bp of one another.

__Calculate Tajima's D using a bins__
```
	tb.py tajima 1000 1000 <vcf>
```

The code above will calculate Tajima's D within 1,000 bp bins across the genome.

## primer

## genome


The __genome__ utility can be used to download genomes from NCBI, <wormbase?>, etc. Downloaded genomes are indexed with bwa, samtools (faidx) and blast. 

__Output reference filename for a given VCF__

```
  tb.py genome <vcf>
```

__Search for genomes__

```
  tb.py genome --search=<term>
```

Search NCBI, etc. for genomes.

__Download genomes__

```
  tb.py genome --download=<asm_name> [--fix-chrom-names]
```

To download a genome, specify its assembly name (`asm_name`) as provided from search results. Use `--fix-chrom-names` to replace NCBI chromosome names with more appropriate roman numeral or numeric chromosome names.

## phylo

## geno

The `geno` command offers functions for modifying genotypes.

###### transfer-filter

```
	tb.py geno transfer-filter <vcf>
```

Generate a FORMAT field filter column (GF) and transfer filters applied within the FILTER column to the GF field. This utility is useful if variants are called individually (producing a single VCF) for samples and then later merged into a multi-sample VCF. It enables you to track which variant filters were applied to a specific sample rather than across a population.

###### het-polarization

```
	tb.py geno het-polarization <vcf>
```

Creates a new FORMAT field (HP) and "polarizes" or switches heterozygous genotypes based on genotype likelyhoods (GL) or Phred-scaled genotype likelihoods (PL). For example, a variant with 7 reads supporting a reference call and 1 read supporting an alternate allele might be called as a heterozygous genotype. However, in hermaphroditic species (_e.g. C. elegans_) where limited or no heterozygosity exists, it is much more likely that the genotype is homozygous reference. The HP FORMAT field will list one of the following depending on the change made:

* __AA__ - A heterozygous genotype was switched to a homozygous reference genotype.
* __AB__ - A heterozygous genotype remains a heterozygous genotype.
* __BB__ - A heterozygous genotype was switched to a homozygous alternate genotype.



## freq

Calculates the frequency of homozygous genotypes by sample (e.g. number of singletons, doubletons, tripletons, etc. by sample)

```
	tb.py freq <vcf>
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

## vcf2tsv

Generate long or wide format tab-separated-values (tsv) from a VCF.

```
	tb.py vcf2tsv (wide|long) [--print-header --snpeff] <vcf>
```

* __wide | long__ - Select one of __wide__ or __long__ to set the output format. 
* __--print-header__ - Print a header row specifying column names.
* __--ANN__ - Parse annotation fields. Results in 1 row x (n)annotations x (n)samples when long or 1 row x (n)annotations x (n)variants when wide.

### Roadmap

* [ ] Check for BWA, Blast, Primer3, other req'd CLI tools.

#### Tajima

* [ ] Tajima-D
	* [ ] Documentation

#### Primer

Suite of tools for genotyping: via sanger sequencing, using snip-SNPs, and indels. Generates appropriate primers and predicts band sizes for indels and snip-SNPs.

* [ ] Primer Design 
	* [ ] snip-SNP
	* [ ] Indels

#### Genome 

Manages genomes used for generating primers and other tasks. Performs indexing for all necessary tools.

* [ ] Check that tools (bwa/samtools/blast) are available.
	* [ ] Skip bwa if not available
* [ ] Error Checking

#### phylo

Tools for producing files useful in constructing phylogies.

#### freq

Generates the allele frequency on a per-sample basis. 

* [X] Calc. freq of individual samples.
* [ ] Heterozygous genotypes

#### tstv

* [ ] Calc. tstv ratios and stratify by other variables (e.g. depth)
* [ ] Visualization in R?

#### geno

Utilities for working with genotypes.

* [X] - Transfer filter
* [X] - Heterozygous polarization.

#### vcf2tsv

Converts a VCF into a tsv - in wide or long format, and taking into account annotation fields (ANN) added by programs such as SNPeff.

* [X] Wide
* [X] Long
* [X] Parse ANN Fields (e.g. snpeff)
* [ ] Read from stdin
* [ ] set fields (INFO / FORMAT)
* [ ] generate bigquery schema and output script

#### Possible additions?

* [ ] Upload to UCSC
* [ ] Annotate variants (track-based search)

