# Overview

The calc command can be used to calculate or count genotypes/alleles from a VCF.

Calculates the frequency of homozygous genotypes (e.g. number of singletons, doubletons, tripletons, etc.) by sample.

```
  vk freq sample_hom_gt <vcf>
  vk freq genotypes [--frequency] <vcf>
  vk freq spectrum <vcf>
```

Output from this utility appears as the table below. The first column indicates the name of the sample. The second column indicates how often the genotype occurs in the population (within the VCF). Finally, the third column indicates the number of cases that exist for that sample. Therefore, the first row indicates that __ED3052__ has two singletons. The second row indicates that __ED3052__ shares a single genotype with another sample for three variants.

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
