# Overview

The filter command can be used to hard or soft-filter variants based on the count or frequency of homozygous reference, heterozygous, homozygous alternative, or missing variants from a  VCF.

```
    vk filter (REF|HET|ALT|MISSING) [--min=<min> --max=<max> --soft-filter=<soft> --mode=<mode>] <vcf>
```

### Options

* __--min=<min>__ - Set a minimum threshold for the specified filter. An integer is interpretted as a minimum count. A float is interpretted as a minimum frequency.
* __--max=<max>__ - Set a maximum threshold for the specified filter. An integer is interpretted as a maximum count. A float is interpretted as a maximum frequency.
* __--soft-filter=<soft>__ - When sepcified, the filter is added to the `FILTER` column and variants are __NOT__ removed.
* __--mode=(+|x)__ - Specifies the mode for soft-filtering. When set to `+`, the filter is added to existing filters in the FILTER column. When set to `x` the soft-filter replaces existing filters in the FILTER column.

### Examples

The command below will return only variants with at least one homozygous REF and one homozygous ALT call:

```
    vk filter REF --min=1 <vcf> | vk filter ALT --min=1 - | bcftools view -O z > one_homozygous.vcf.gz
```

The command below will soft-filter variants with greater than 10% missing calls.

```
    vk filter MISSING --max=0.90 --soft-filter=HIGH-MISSING --mode=x <vcf>
```

Return variants with a no more than 3 heterozygous calls:

```
    vk filter HET --max=3 <vcf>
```