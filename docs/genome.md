# Usage

The `vk genome` command can be used to download and prepare reference genomes for use with other tools within `vcf-kit`. `vcf-kit` will do the following when downloading a reference genome:


1. Download the reference genome.
1. If downloaded from NCBI, `vcf-kit` will attempt to replace chromosomes names with shorthand descriptors of if possible (_e.g._ Chromosome I).
1. Unzip the gzipped reference and re-compress using `bgzip`.
1. Create a `bwa` index.
1. Create a `samtools` index.
1. Create a `blast` database.
1. Remove temporary file names.


## Setup Genomes

### View/set genome directory

By default, genomes are stored within your home directory in a `.genome` folder. The location of this directory can be viewed by typing:

```
vk genome location
```

Similarly, you can change the location by adding a path:

```
vk genome /path-to-my-new-genome-directory
```

### List genomes

A list of downloaded genomes can be viewed using:

```bash
vk genome list
```

## NCBI Genomes

`vcf-kit` makes it easy to obtain and prepare genomes from the [NCBI genome database](http://www.ncbi.nlm.nih.gov/genome/). To do this, it downloads a text file containing a list of all available genomes and uses this for searching purposes. To search for a genome, you can type:

```
vk genome --search="Human cyclovirus"
```

The results of the search will be output in a table:

```shell

  Genome Directory: /Users/dancook/.genome


  Searching...

    assembly_accession    bioproject    organism_name               asm_name         ftp_path
    --------------------  ------------  --------------------------  ---------------  ----------------------------------------------------------------------
    GCF_000908835.1       PRJNA209365   Human cyclovirus VS5700009  ViralProj209365  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000908835.1_ViralProj209365
    GCF_000918035.1       PRJNA243497   Human cyclovirus            ViralProj243497  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000918035.1_ViralProj243497

  To download a genome and setup for use:

    vk genome --ref=<asm_name>
```

As the instructions illustrate, you can download the genome by providing the `asm_name` (assembly name).

### Downloading Genomes from RefSeq

```
vk genome --ref=ViralProj209365
```

## Wormbase Genomes

Reference genomes can also be obtained from wormbase. 

```
vk genome wormbase --ref=WS245
```