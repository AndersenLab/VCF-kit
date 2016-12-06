# Overview

The `call` command can be used to compare variants identified from Sanger sequencing with those present within a VCF. The `call` command uses BLAST to perform variant calling. This is a not a sophisticated method of variant calling, and is merely meant to facilitate in comparing Sanger sequencing with variants called in a VCF. Only homozygous calls are supported. Sanger sequences can be provided as AB1/ABI, FASTQ, or FASTA files. Sequences must be annotated with the sample name if you want to compare the genotype calls.

## Annotating Samples

VCF-kit first tries to identify the sample from sequences using their filename. Because `ABI/AB1` files can only contain one sequence, you must specify their sample using the filename. To identify the sample, VCF-kit searches for the sample name after splitting the string using the `\W` regular expression. For example:

The following all work if the `DL238` sample is located in your VCF:
```
DL238.AB1
DL238.fa
DL238.fastq.gz
DL238.20161205.fq.gz
DL238_20161205.fq.gz
DL238-20161205.fq.gz
```

In the first line, `DL238.AB1` is split into `DL238` and `AB1`. Note that among these split pieces, the name of a sample must match exactly.

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
* __--vcf-targets__ - Only outputs sites present within the VCF.

## Output

The following columns are output when using the `call` command:

* __CHROM__ - chromosome
* __POS__ - position
* __REF__ - The alternative genotype within the VCF.
* __ALT__ - The reference genotype within the VCF.
* __seq_gt__ - The Sanger genotype for the sample listed.
* __vcf_gt__ - The VCF genotype for the sample listed.
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

```
vk call DL238_sanger.AB1 --ref=WBcel235 --all-sites illumina_sequencing.vcf.gz
```

| CHROM   |      POS | REF   | ALT   | seq_gt   | vcf_gt   | sample   | variant_type   | classification   |   index |   alignment_start |   alignment_end | strand   | context           |   gaps |   mismatch |   evalue |   bitscore |
|:--------|---------:|:------|:------|:---------|:---------|:---------|:---------------|:-----------------|--------:|------------------:|----------------:|:---------|:------------------|-------:|-----------:|---------:|-----------:|
| X       | 14556961 | A     | None  | A        |          | DL238    |                |                  |       3 |          14556961 |        14557595 | +        | [A]TTCACTTTT      |      6 |          5 |        0 |       1112 |
| X       | 14556962 | T     | None  | T        |          | DL238    |                |                  |       4 |          14556961 |        14557595 | +        | A[T]TCACTTTTA     |      6 |          5 |        0 |       1112 |
| X       | 14556963 | T     | None  | TC       |          | DL238    | insertion      | FN               |       5 |          14556961 |        14557595 | +        | ATT[TC]ACTTTTATT  |      6 |          5 |        0 |       1112 |
| X       | 14556964 | A     | None  | A        |          | DL238    |                |                  |       7 |          14556961 |        14557595 | +        | ATTC[A]CTTTTATTC  |      6 |          5 |        0 |       1112 |
| X       | 14556965 | C     | None  | C        |          | DL238    |                |                  |       8 |          14556961 |        14557595 | +        | ATTCA[C]TTTTATTCT |      6 |          5 |        0 |       1112 |                                                                                                                 |

```
CHROM   POS REF ALT seq_gt  vcf_gt  sample  variant_type    classification  index   alignment_start alignment_end   strand  context gaps    mismatch    evalue  bitscore    phred_quality   phred_quality_window    description
X   14556963    T   None    TC      DL238   insertion   FN  5   14556961    14557595    +   ATT[TC]ACTTTTATT    6   5   0.0 1112
X   14556984    AA  None    A       DL238   deletion    FN  25  14556961    14557595    +   TTTTTAGGAA[A]TTGCAAGT   6   5   0.0 1112
X   14556989    CA  None    C       DL238   deletion    FN  29  14556961    14557595    +   AGGAATTGC[C]AAGTCTTAC   6   5   0.0 1112
X   14557228    G   None    G   G   DL238       TN  269 14556961    14557595    +   TTGTTGGTTC[G]CCAGGATCG  6   5   0.0 1112
X   14557388    T   None    T   T   DL238       TN  429 14556961    14557595    +   ATTCGCTGGG[T]CCGGCCTCC  6   5   0.0 1112
X   14557506    A   None    G   G   DL238   snp TP  547 14556961    14557595    +   TATACTAGCT[G]CATAGACAA  6   5   0.0 1112
X   14557521    T   None    T   T   DL238       TN  562 14556961    14557595    +   GACAACTGAC[T]GTGTATGTG  6   5   0.0 1112
X   14557572    T   None    AA      DL238   insertion   FN  613 14556961    14557595    +   ATTGATCTCA[AA]GAAAAAATA 6   5   0.0 1112
X   14557576    C   None    A       DL238   snp FN  618 14556961    14557595    +   ATCTCAAGAA[A]AAATAACCG  6   5   0.0 1112
X   14557578    T   None    AA      DL238   insertion   FN  620 14556961    14557595    +   TCAAGAAAAA[AA]TAACCGCAC 6   5   0.0 1112
X   14557586    T   None    A       DL238   snp FN  629 14556961    14557595    +   AAATAACCGC[A]CGCTGGAA   6   5   0.0 1112
X   14557591    TT  None    T       DL238   deletion    FN  632 14556961    14557595    +   ACCGCACGCT[T]GGAA   6   5   0.0 1112
```
