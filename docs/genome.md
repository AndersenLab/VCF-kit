# Overview

The `genome` command can be used to download and prepare reference genomes for use with other tools within `vcf-kit`. `vcf-kit` will do the following when downloading a reference genome:


1. Download the reference genome.
1. If downloaded from NCBI, `vcf-kit` will attempt to replace chromosomes names with shorthand descriptors of if possible (_e.g._ Chromosome I).
1. Unzip the gzipped reference and re-compress using `bgzip`.
1. Create a `bwa` index.
1. Create a `samtools` index.
1. Create a `blast` database.
1. Remove temporary file names.


# Usage

## View/set genome directory

By default, genomes are stored within your home directory in a `.genome` folder. The location of this directory can be viewed by typing:

```
vk genome location
```

!!! tip

    You can use a bash alias to access the genome currently set with vcf-kit. Add this to your __.bash_profile__:
    ```
    alias GENOME=`vk genome location`
    ```
    And you can access the currently set genome using `GENOME`.

Additionally, you can change the location by adding a path:

```bash
vk genome /path-to-my-new-genome-directory
```

## List genomes

A list of downloaded genomes can be viewed using:

```bash
vk genome list
```

## Download Genomes

### Search NCBI

`vcf-kit` makes it easy to obtain and prepare genomes from the [NCBI genome database](http://www.ncbi.nlm.nih.gov/genome/). To do this, it downloads a [text file](http://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt) containing a list of all available genomes and uses this for searching purposes. To search for a genome, you can type:

```
vk genome --search="Human cyclovirus"
```

The results of the search will be output in a table:

```shell

  Genome Directory: /Users/dancook/.genome


  Searching...

    assembly_accession    bioproject    organism_name               asm_name         ...
    --------------------  ------------  --------------------------  ---------------  ...
    GCF_000908835.1       PRJNA209365   Human cyclovirus VS5700009  ViralProj209365  ...
    GCF_000918035.1       PRJNA243497   Human cyclovirus            ViralProj243497  ...

  To download a genome and setup for use:

    vk genome --ref=<asm_name>
```

As the instructions illustrate, you can download the genome by providing the `asm_name` (assembly name).

### Download from NCBI

Set `--ref` to an `asm_name` from the search results table to download a genome.

```bash
vk genome --ref=ViralProj209365
```

### Wormbase

Reference genomes can also be obtained from wormbase. 

```bash
vk genome wormbase --ref=WS245
```

# Custom Directories

It is possible to set the directory to download a genome using the `--directory` parameter.

```bash
vk genome ncbi --directory="." --ref=ViralProj15089
```
