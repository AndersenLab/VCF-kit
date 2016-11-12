# Overview

```
  vk calc sample_hom_gt <vcf>
  vk calc genotypes [--frequency] <vcf>
  vk calc spectrum <vcf>
```

The `calc` command can be used examine the frequency or count of genotypes/alleles from a VCF in different ways.

## Calculate shared homozygous genotypes

```
    vk calc sample_hom_gt <vcf>
```

The above command calculates the frequency of homozygous genotypes (e.g. number of singletons, doubletons, tripletons, etc.) by sample.

An example of the output from this utility appears as the table below. The first column indicates the name of the sample. The second column indicates how often the genotype occurs in the population (within the VCF). Finally, the third column indicates the number of cases that exist for that sample. Therefore, the first row indicates that __ED3052__ has two singletons. The second row indicates that __ED3052__ shares a single genotype with another sample for three variants.


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


## Calculate genotype distributions

The following command will output a table that looks like this:

```
    vk calc genotypes <vcf>
```

Providing the following output:

|   n |   ref |   het |   alt |   mis |
|----:|------:|------:|------:|------:|
| 937 |    14 |     0 |     0 |     0 |
| 328 |    13 |     0 |     0 |     1 |
| 242 |    13 |     0 |     1 |     0 |
| 168 |    12 |     0 |     0 |     2 |
| 101 |    11 |     0 |     0 |     3 |
|  94 |    12 |     0 |     1 |     1 |
|  89 |    12 |     0 |     2 |     0 |
|  73 |    10 |     0 |     0 |     4 |

Where:

* __n__ = number of variants/records with the genotype distribution shown.
* __ref__ - number of homozygous reference genotypes
* __het__ - number of heterozygous genotypes
* __alt__ - number of homozygous alternative genotypes (biallelic variants only)
* __mis__ - number of missing genotypes

In the example above, the first row would represent 937 records where there are 14 reference genotypes and 0 heterozygous, 0 homozygous alternative, and 0 missing genotypes.

Optionally, you can specify the `--frequency` flag to calculate the frequencies of each genotype.

```
    vk calc genotypes --frequency <vcf>
```

|   n |      ref |   het |       alt |       mis |
|----:|---------:|------:|----------:|----------:|
| 937 | 1        |     0 | 0         | 0         |
| 328 | 0.928571 |     0 | 0         | 0.0714286 |
| 242 | 0.928571 |     0 | 0.0714286 | 0         |
| 168 | 0.857143 |     0 | 0         | 0.142857  |
| 101 | 0.785714 |     0 | 0         | 0.214286  |
|  94 | 0.857143 |     0 | 0.0714286 | 0.0714286 |
|  89 | 0.857143 |     0 | 0.142857  | 0         |
|  73 | 0.714286 |     0 | 0         | 0.285714  |
|  62 | 0.785714 |     0 | 0.214286  | 0         |
|  46 | 0.714286 |     0 | 0.285714  | 0         |


## Calculate the frequency of alleles

```
    vk calc spectrum <vcf>
```

The above command will generate the output below:

|   n |   alt_allele_freq |
|----:|------------------:|
|   3 |          1        |
|  10 |          0.928571 |
|   1 |          0.923077 |
|   1 |          0.9      |
|   7 |          0.857143 |
|   1 |          0.846154 |
|   4 |          0.833333 |
|   2 |          0.818182 |
|   5 |          0.785714 |
|   1 |          0.777778 |
|   1 |          0.769231 |
|   7 |          0.75     |
|   7 |          0.714286 |

This provides the number of records and their corresponding allele frequency for alternative alleles only. All alleles that are non-reference are counted.


