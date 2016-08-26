# Overview

The `phylo fasta` command can be used to generate a fasta file of variants and enables users to generate a tree file (in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html)) using [UPGMA](https://en.wikipedia.org/wiki/UPGMA) or [neighbor-joining](https://en.wikipedia.org/wiki/Neighbor_joining). It can also generate a plot of the tree/phylogeny.

`phylo fasta` can read from stdin by using `-`.

## Generate fasta-alignment from variant calls

The `phylo fasta` command generates a fasta file by concatenating all single-nucleotide variants from a VCF for each sample. Missing values are replaced with a `-`. The fasta sequence is constructed in memory - so it can be somewhat resource intensive for larger files.

```
vk phylo fasta <vcf>
```

__Output__:

```
>QG536 <-- Sample 1
AGGGATCCT-GGG...
>GXW1 <-- Sample 2
AGAGATCCCTGGG...
>DL200 <-- Sample 3
AGAGA-CCCTGG-...
```

## Generate fasta-alignment from stdin

You may be interested in filtering variants prior to generating a fasta file. The command below reads from stdin.

```
bcftools filter --set-GTs . --exclude 'FMT/DP < 20'  data/test.vcf.gz | vk phylo fasta -
```

## Generating a tree/phylogeny

`vk phylo tree` can be used to generate a tree/phylogeny from a vcf file. This command uses a fasta file (identical to what is produced using `vk phylo fasta`), and uses [MUSCLE](https://en.wikipedia.org/wiki/MUSCLE_(alignment_software) to produce a tree file. 

### Generate a upgma tree

An unweighted pair group method with arithmetic-mean (UPGMA) tree can be constructed using the following command. Output is in newick format.

```
vk phylo tree upgma <vcf>
```

### Generate a neighbor-joining tree

An neighbor-joining tree can be constructed using the following command. Output is in newick format.

```
vk phylo tree nj <vcf>
```

Generate fasta sequences from variant data. This is useful for generating phylogenetic trees from VCF files.

### Output format

__Output - Newick format__

The `phylo tree` command sends output to stdout  is newick format. Newick format can be used 

```
(((N2:0.0250154,PX179:0.02262):0.00270637,(((((EG4946:0.035835,AB1:0.0349638):0.00435886,GXW1:0.0490124):0.00222221,(((WN2001:0.0850733,CB4856:0.130009))...
```

### Plot a phylogeny from a VCF file

`phylo tree` can be used to generate a plot of a phylogeny by adding the `--plot` flag. 

```
vk phylo tree nj --plot <vcf>
```

![phylogeny example](https://github.com/AndersenLab/vcf-toolbox/raw/img/tb_phylo.png)
