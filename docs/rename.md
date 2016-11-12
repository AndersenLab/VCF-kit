# Overview

The `rename` command can be used to add a prefix or suffix to sample names, or to substitute a sample names. Note that the command outputs an entirely new VCF to stdout. Therefore, it is best to use it within a pipe prior to output.

```
usage:
  vk rename [--prefix=<prefix> --suffix=<suffix> --subst=<subst>] <vcf>
```

# Options

* `--prefix`
* `--suffix`
* `--subst`

# Examples

When comparing VCFs, you may want to make clear the source of the samples - hence the need to add a prefix or suffix. `vk rename` can be used to add a prefix to all sample names using:

```
vk rename --prefix="GATK_" <vcf/stdin>
```

A suffix can be added to all sample names similarly:

```
vk rename --suffix="_BCFTOOLS" <vcf/stdin>
```

You can substitute the name of a single sample using `--subst`. Use a `:`, `,`, or `=` to separate what you are substituting from and to:

```
vk rename --subst=OLD_SAMPLE_NAME:NEW_SAMPLE_NAME
```

The above command would rename the `N2` sample as `REFERENCE_SAMPLE`.

Multiple sample names can be substituted. as well.

```
vk rename --subst=N2:REFERENCE_SAMPLE --subst WN2001=ANOTHER_RENAME
```

Rename commands can be combined as well. Substitutions are performed first, followed by adding a prefix, and finally adding suffixes.

```
vk rename --prefix PREFIX_ --suffix _SUFFIX --subst=N2:REFERENCE_SAMPLE
```

