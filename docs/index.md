# VCF-kit

![logo.png](logo.png)

## About

The `Variant Call Format Kit` is a collection of tools useful for performing a variety of analyses and operations on variant data stored using the [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format). VCF-kit is open source and community contributions are encouraged.

VCF-kit is open source and is licensed under [the MIT License](https://raw.githubusercontent.com/AndersenLab/VCF-kit/master/LICENSE). We welcome community contributions.


## Installation

__VCF-Kit has been upgraded to Python 3__

VCF-kit has been tested with Python 3.6. VCF-kit makes use of additional software for a variety of tasks:

* bwa (v 0.7.12)
* samtools (v 1.3)
* bcftools (v 1.3)
* blast (v 2.2.31+)
* muscle (v 3.8.31)
* primer3 (v 2.5.0)

You can install these dependencies and VCF-kit using conda, or you can use a Docker image.

#### Conda

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n vcf-kit \
  danielecook::vcf-kit=0.2.6 \
  "bwa>=0.7.17" \
  "samtools>=1.10" \
  "bcftools>=1.10" \
  "blast>=2.2.31" \
  "muscle>=3.8.31" \
  "primer3>=2.5.0"

conda activate vcf-kit
```

#### Docker

You can also run VCF-kit with all installed dependencies using docker: 

```bash
docker run -it andersenlab/vcf-kit vk
```