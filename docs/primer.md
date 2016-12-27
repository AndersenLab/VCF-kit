# Overview

The `primer` command can be used to generate primers for the purpose of verifying genotype calls present within a VCF. Primers can be generated for verifying genotypes using snip-SNPs, indels (based on product size), and using sanger sequencing.

All of the `primer` subcommands require you to specify a reference using `--ref`. A reference can be obtained using the `vk genome` command.

## Options

* __--ref__ - Set the reference genome to use. Obtain genomes using `vk genome`.
* __--region__ - Restrict primer generation to a specific region.
* __--samples__ - Output genotypes for a sample or set of samples.
* __--template__ - The template to use for generating primers.
* __--size__ - For the `template` command, sets the size of the region output, specified as an integer (e.g. 300). For  `sanger`, sets the target amplicon size, specified as a range (e.g. 400-600).
* __--polymorphic__ - Only outputs variants where there is at least one 0/0 and one 1/1 genotype.
* __--nprimers__ - Specify the number of primer sets desired for each variant.

## Output

Output is in tab-seperated value (TSV) format. The following columns are always output:

* __CHROM__ - chromosome/contig
* __POS__ - chromosomal position
* __region__ - The region used as a template to generate primers. Note that the amplicons will be smaller
* __REF__ - Reference allele
* __ALT__ - Alternative allele
* __template_sample__ - Template used for generating primers. Default is to use ALT calls.

These additional columns are output for `snip`, `indel`, and `sanger` options:

* __variant_count__ - Number of variable sites within interval __among the samples specified__. Includes missing/heterozygous sites. A high variant count should be avoided. 
* __primer_left__ - Left PCR primer
* __primer_right__ - Right PCR primer
* __melting_temperature__ - Melting temperature for the left and right primers (TM), separated by a comma.
* __amplicon_length__ - Sequence of the amplified piece of DNA
* __amplicon_region__ - Genomic region (Chromsome:Start-End) of the amplified region.
* __amplicon_sequence__ - Sequence of the amplified region.
* __0/0__ - Comma-delimited samples with homozygous reference genotypes
* __1/1__ - Comma-delimited samples with homozygous alternative genotypes
* __polymorphic__ - True when there is at least one 0/0 and one 1/1 genotype.

Additionally, there are columns specific to the type of primers being generated.

### snip output columns

* __ref_sites__ - Cut Position : Product Sizes for Reference genotype (see example below)
* __alt_sites__ - Cut positions : Product Sizes for Alternative genotype (see example below)
* __restriction_enzyme__ - The restriction enzyme to use.
* __restriction_site__ - The sequence/motif used by the restriction enzyme
* __restriction_site_length__ - Length of the restriction site

__ref_sites/alt_sites__ are denoted as the following example illustrates:

ref_sites: `183:183,472` - There is one restriction site at position 183; Cutting that site produces two products. One is `183` bp, the other is `472` bp.

The alternative genotypes produce a different set of restriction sites and product sizes:

`183,258:183,75,397` - Cut at `183` and `258`; Product sizes are `183`, `75`, and `397`.

### indel output columns

* __indel_size__ - basepair size of interval
* __indel_type__ - Deletion or Insertion relative to the reference
* __REF_product_size__ - PCR product size with homozygous reference genotype
* __ALT_product_size__ - PCR product size with homozygous alternative genotype

### sanger output columns

__variant_distance__ - Indicates the distance from the start of the amplicon to the variant position. 

## Primer3 Records

The section below describes how primers are generated using Primer3. When __VCF-kit__ generates primers, it creates a primer configuration file in the [Boulder-IO](http://primer3.sourceforge.net/primer3_manual.htm) format. The options are set as follows. 

The options below are always used:
```
PRIMER_GC_CLAMP=1
PRIMER_MAX_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_TASK=pick_pcr_primers
```

The options below depend on specified options or the type of primers being generated.
```
PRIMER_NUM_RETURN=5 # Can be set using -nprimers
PRIMER_PRODUCT_SIZE_RANGE=600-800 # Can be set for template and sanger options.
```

```
PRIMER_THERMODYNAMIC_PARAMETERS_PATH # set to either /usr/local/share/primer3_config/ or /usr/local/share/primer3/primer3_config/
```

Changes for each variant; The region down- and upstream of the variant:
```
SEQUENCE_TEMPLATE
```

## Commands

## template

Template can be used to fetch the region surrounding variants. Set the size of the desired region using `--size=<int>`.

## snip

[snip-SNPs](http://www.wormbook.org/wli/wbg16.1p28/) are single-nucleotide polymorphisms that modify a restriction site, resulting in a restriction fragment length polymorphism (RFLP). They are useful because they can be used to determine the genotype of a sample using only PCR and restriction enzymes. VCF-kit takes the context of SNPs and determines whether a snip-SNP is present. Then it generates primers and estimates product sizes with and without the restriction site present. snip-SNP primers can be generated using `vk primer snip`.

__--enzymes__

Restriction enzymes can be set using the `--enzymes` option as a tab-delimited list or by specifying one of the groups below:

* `ALL` - All restriction enzymes available in BioPython.
* `Common` - A common set of restriction enzymes.
* `HF` - [High fidelity enzymes as denoted by NEB](https://international.neb.com/products/restriction-endonucleases/hf-nicking-master-mix-time-saver-other/high-fidelity-restriction-enzymes).

```
vk primer snip --ref=WBcel235 --enzymes=DraI data/test.vcf.gz # Specify a single restriction enzyme
vk primer snip --ref=WBcel235 --enzymes=DraI,BssAI <vcf> # Specify a list by delimiting with a comma
vk primer snip --ref=WBcel235 --enzymes=HF <vcf> # Specify a group of enzymes; All/Common/HF
```

__Examples__

```
vk primer snip --ref=WBcel235 data/test.vcf.gz
```

The table has been formatted for easier reading. Normal output is tab-seperated. 

```
 CHROM | POS     | region            | REF | ALT | template_sample | variant_count | ref_sites           | alt_sites                   | restriction_enzyme | restriction_site | restriction_site_length | primer_left          | primer_right         | melting_temperature | amplicon_length | amplicon_region   | amplicon_sequence | 0/0 | 1-Jan | polymorphic |
-------|---------|-------------------|-----|-----|-----------------|---------------|---------------------|-----------------------------|--------------------|------------------|-------------------------|----------------------|----------------------|---------------------|-----------------|-------------------|-------------------|-----|-------|-------------|
 I     | 1198228 | I:1197728-1198727 | G   | A   | ALT             | 1             | 93,283:93,190,471   | 93:93,661                   | NdeI               | CATATG           | 6                       | gtattcagtgggcaagcagc | GGATTAGGCCACCATCCGAG | 59.547,59.965       | 754             | I:1197943-1198697 | gtatt…            |     |       | FALSE       |
 I     | 1487691 | I:1487191-1488190 | A   | C   | ALT             | 1             | 548,654:548,106,100 | 356,548,654:356,192,106,100 | BsuRI              | GGCC             | 4                       | TCAAAGCTGTTTTTGGCGGG | CTTCCCGACAACTTTGCTGC | 59.896,60.04        | 754             | I:1487335-1488089 | TCAAA…            |     |       | FALSE       |
 I     | 1487691 | I:1487191-1488190 | A   | C   | ALT             | 1             | 548,654:548,106,100 | 356,548,654:356,192,106,100 | BsnI               | GGCC             | 4                       | TCAAAGCTGTTTTTGGCGGG | CTTCCCGACAACTTTGCTGC | 59.896,60.04        | 754             | I:1487335-1488089 | TCAAA…            |     |       | FALSE       |
```

__Note__: Amplicon sequences are truncated in the above output.

## Indels

## Sanger

The `sanger` command will generate PCR primers that can be used to amplify a region of interest. The left primer can then also be used to initiate Sanger sequencing. `sanger` can be used to verify both indels and snps.

!!! tip

    Set the amplicon size (`--size`) of the region between 500-800 bp. Your variant should be 50-600 bp upstream of the start of the amplicon as indicated by `variant_distance`.


