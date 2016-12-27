# Overview

The `call` command can be used to compare variants identified from Sanger sequencing with those present within a VCF. The `call` command uses BLAST to perform variant calling. This is a not a sophisticated method of variant calling, and is merely meant to facilitate in comparing Sanger sequencing with variants called in a VCF. Only homozygous calls are supported. Sanger sequences can be provided as AB1/ABI, FASTQ, or FASTA files. Sequences must be annotated with the sample name if you want to compare the genotype calls.

## Annotating Samples

VCF-kit first tries to identify the sample from sequences using their filename. Because `ABI/AB1` files can only contain one sequence, you must specify their sample using the filename. To identify the sample, VCF-kit searches for the sample name after splitting the string using the `\W` regular expression. For example:

The following Sanger sequence files (in AB1, fa, and fastq format) are all acceptable ways of specifying that the sample from the VCF is `DL238`:
```
DL238.AB1
DL238.fa
DL238.fastq.gz
DL238.20161205.fq.gz
DL238_20161205.fq.gz
DL238-20161205.fq.gz
```

In the first line, `DL238.AB1` is split into `DL238` and `AB1`. Among the split pieces, the name of a sample must match exactly.

The examples below do not work:

```
testDL238.AB1
testDL23820161205.AB1
```

You may have multiple samples within a single fastq or fasta file. In this case, you can use the description line to annotate each record:

```
> DL238
CGTCAAGGGGCACCGGATCG...
```

```
> @FCC1GJUACXX:6:1101:19686:24985#AAGGTCTA/2 |DL238|
CGTCAAGGGGCACCGGATCG...
```

And multiple samples may be present within the same file:

```
> DL238
CGTCAAGGGGCACCGGATCG...
> N2
CGTCAAGGCGCACCGGATCG...
> MY25
CGTCAAGGCGCACCGGATCG...
```

## Usage

```
    vk call <seq> --ref=<reference> [--all-sites --vcf-targets <vcf>]
```

Perform variant calling using blast. Useful for validating variants using sanger sequencing.  

* __&lt;seq&gt;__ - Fasta (fa), Fastq (fq, fastq), or ab1 format, determined by file extension containing sanger reads.
* __--ref__ - The reference genome used to perform calling.
* __--all-sites__ - Outputs comparison of every site that aligns within the alignment.
* __--vcf-sites__ - Only outputs sites present within the VCF.

## Output

The following columns are output when using the `call` command:

* __CHROM__ - chromosome
* __POS__ - position
* __reference__ - The base at the genomic position within the reference genome. Should always be identical to `REF`, when `REF` is listed.
* __REF__ - The reference genotype within the VCF. Only listed when a variant is present.
* __ALT__ - The alternative genotype within the VCF. Only listed when a variant is present.
* __seq_gt__ - The Sanger genotype for the sample listed.
* __vcf_gt__ - The VCF genotype for the sample listed. Homozygous calls are listed as a single base.
* __sample__ - The sample being compared. Must be matched to the Sanger sequence and present within the VCF.
* __variant_type__ - The type of variant observed as compared to the Sanger sequence: SNV, Insertion, or Deletion.
* __classification__ - True/False Positive/Negative, assuming the Sanger sequence is correct.
* __index__ - The index of the Sanger sequence, as aligned.
* __alignment_start__ - The start of the blast alignment.
* __alignment_end__ - The end of the blast alignment.
* __strand__ - The strand aligned to (`+` or `-`).
* __gaps__ - Total number of gaps in the alignment.
* __mismatch__ - Total number of mistmatches in the alignment.
* __bitscore__ - The bitscore of the alignment.


## Examples

<small>__Note:__ Some columns removed from examples for conciseness.</small>

__--all-sites__

Shows all calls from Sanger and VCF. Will output every position that was aligned to reference from Sanger sequence and corresponding VCF calls.
Notice the absence of REF/ALT for sites not called in the VCF below.

```
vk call DL238_sanger.AB1 --ref=WBcel235 --all-sites illumina_sequencing.vcf.gz
```

| CHROM   |      POS | reference   | REF   | ALT   | seq_gt   | vcf_gt   | sample   | variant_type   | classification   |   index |   alignment_start |   alignment_end |
|:--------|---------:|:------------|:------|:------|:---------|:---------|:---------|:---------------|:-----------------|--------:|------------------:|----------------:|
| X       | 14557504 | C           |       |       | C        |          | DL238    |                |                  |     545 |          14556961 |        14557595 |
| X       | 14557505 | T           |       |       | T        |          | DL238    |                |                  |     546 |          14556961 |        14557595 |
| X       | 14557506 | G           | A     | G     | G        | G        | DL238    | snp            | TP               |     547 |          14556961 |        14557595 |
| X       | 14557507 | C           |       |       | C        |          | DL238    |                |                  |     548 |          14556961 |        14557595 |

__--vcf-targets__

Only show variants present in the VCF.

```
vk call DL238_sanger.AB1 --ref=WBcel235 --vcf-sites illumina_resequencing.vcf.gz
```

| CHROM   |      POS | reference   | REF   | ALT   | seq_gt   | vcf_gt   | sample   | variant_type   | classification   |   index |   alignment_start |   alignment_end |
|:--------|---------:|:------------|:------|:------|:---------|:---------|:---------|:---------------|:-----------------|--------:|------------------:|----------------:|
| X       | 14557228 | G           | G     | A     | G        | G        | DL238    | snp            | TN               |     269 |          14556961 |        14557595 |
| X       | 14557388 | T           | T     | C     | T        | T        | DL238    | snp            | TN               |     429 |          14556961 |        14557595 |
| X       | 14557506 | A           | A     | G     | G        | G        | DL238    | snp            | TP               |     547 |          14556961 |        14557595 |
| X       | 14557521 | T           | T     | A     | T        | T        | DL238    | snp            | TN               |     562 |          14556961 |        14557595 |

