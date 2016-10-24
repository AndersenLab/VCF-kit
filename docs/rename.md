# Overview

The `rename` command can be used to add a prefix or suffix to sample names, or to substitute a string within a sample name. Note that because it is difficult to edit VCF files in place - this command requires rewriting the entire file. Therefore, it is best to use it within a pipe prior to output.

```
usage:
  vk rename [--prefix=<prefix> --suffix=<suffix> --subst=<subst>] <vcf>
```

# Options

* `--prefix`
* `--suffix`
* `--subst`

# Examples

When comparing VCFs, you may want to make clear the source of the samples - hence the need to add a prefix or suffix. `vk rename` can be used to add a prefix using:

```
vk rename --prefix="GATK_" <vcf/stdin>
```

A suffix can be added similarly:

```
vk rename --suffix="_BCFTOOLS" <vcf/stdin>
```

Finally, you can substitute the name of a single sample (or part of the name) using `--subst`. Use a `:` to separate what you are substituting from and to:

```
vk rename --subst=N2:REFERENCE_SAMPLE
```
