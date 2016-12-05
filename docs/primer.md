# Overview

The `primer` command can be used to generate primers for the purpose of verifying genotype calls present within a VCF. Primers can be generated for verifying genotypes using snip-SNPs, indels (based on product size), and for sanger sequencing.

Output is in tab-seperated value (TSV) format, and the columns below are output regardless of method used:

* __CHROM__ - chromosome/contig
* __POS__ - chromosomal position
* __region__ - The region used as a template to generate primers. Note that the amplicons will be smaller
* __REF__ - Reference allele
* __ALT__ - Alternative allele
* __primer_left__ - Left PCR primer
* __primer_right__ - Right PCR primer
* __amplicon_length__ - Sequence of the amplified piece of DNA
* __amplicon_region__ - Genomic region (Chromsome:Start-End) of the amplified region.
* __amplicon_sequence__ - Sequence of the amplified region.
* __0/0__ - Comma-delimited samples with homozygous reference genotypes
* __1/1__ - Comma-delimited samples with homozygous alternative genotypes
* __polymorphic__ - True when there is at least one 0/0 and one 1/1 genotype.

All of the `primer` subcommands require you to specify a reference using `--ref`. A reference can be obtained using the `vk genome` command.

## Options

* __--ref__ - Set the reference genome to use. Obtain genomes using `vk genome`.
* __--region__ - Restrict primer generation to a specific region.
* __--samples__ - Output genotypes for a sample or set of samples.
* __--size__ - Sets the target amplicon size, specified as a range (e.g. 400-600).
* __--polymorphic__ - Only outputs variants where there is at least one 0/0 and one 1/1 genotype.
* __--nprimers__ - Specify the number of primer sets desired for each variant.


### snip-SNP

[snip-SNPs](http://www.wormbook.org/wli/wbg16.1p28/) are single-nucleotide polymorphisms that modify a restriction site, resulting in a restriction fragment length polymorphism (RFLP). They are useful because they can be used to determine the genotype of a sample using only PCR and restriction enzymes. VCF-kit takes the context of SNPs and determines whether a snip-SNP is present. Then it generates primers and estimates product sizes with and without the restriction site present. snip-SNP primers can be generated using `vk primer snip`.

__--enzymes__

Restriction enzymes can be set using the `--enzymes` option as a tab-delimited list or by specifying one of the groups below:

* `ALL` - All restriction enzymes available in BioPython.
* `Common` - A common set of restriction enzymes.
* `HF` - [High fidelity enzymes as denoted by NEB](https://international.neb.com/products/restriction-endonucleases/hf-nicking-master-mix-time-saver-other/high-fidelity-restriction-enzymes).

Lists of restriction enzymes to use can be specified using a comma delimited list:

Specify a group of restriction enzymes to use:

```
vk primer snip --enzymes=HF ...
```

Specify a list of specific restriction enzymes to use:

```
vk primer snip --enzymes=DraI,BssAI ...
```

__Examples__

```
vk primer snip --ref=WBcel235 --enzymes=DraI data/test.vcf.gz
```

##### snip-SNP output

In addition to the common columns listed above, snip-SNP output includes the following:

* __variant_count__ - Number of variants within interval
* __ref_sites__ - Cut Position : Product Sizes for Reference genotype
* __alt_sites__ - Cut positions : Product Sizes for Alternative genotype
* __restriction_enzyme__ - The restriction enzyme to use.
* __restriction_site__ - The sequence/motif used by the restriction enzyme
* __restriction_site_length__ - Length of the restriction site

## Indels

##### Indel output

* __indel_size__ - basepair size of interval
* __indel_type__ - Deletion or Insertion relative to the reference
* __REF_product_size__ - PCR product size with homozygous reference genotype
* __ALT_product_size__ PCR product size with homozygous alternative genotype

## Sanger

The `sanger` command will generate PCR primers that can be used to amplify a region of interest. The left primer can be used to perform Sanger sequencing.

__variant_distance__ - Indicates the distance from the start of the amplicon to the variant position. 

!!! tip

    Set the amplicon size (`--size`) of the region between 500-800 bp. Your variant should be 50-600 bp upstream of the start of the amplicon as indicated by `variant_distance`.


