FROM mambaorg/micromamba:1.5.8-alpine3.20

RUN micromamba create -n vcf-kit -y -c bioconda -c conda-forge -c anaconda \
    biopython==1.77

RUN micromamba install -n vcf-kit -y -c bioconda -c conda-forge -c anaconda \
    "bwa>=0.7.17" \
    "samtools>=1.10" \
    "bcftools>=1.10" \
    "blast>=2.2.31" \
    "muscle>=3.8.31" \
    "primer3>=2.5.0" \
    setuptools

RUN micromamba install -n vcf-kit -y -c bioconda -c conda-forge -c anaconda \
    awesome-slugify \
    matplotlib \
    scipy \
    numpy \
    cython \
    cyvcf2 \
    docopt

RUN micromamba install -n vcf-kit -y -c bioconda -c conda-forge -c anaconda \
    logzero \
    pomegranate \
    clint \
    requests \
    networkx \
    intervaltree \
    tabulate \
    jinja2 \
    pytest \
    pytest-runner

RUN micromamba install -n vcf-kit -y -c bioconda -c conda-forge -c anaconda pomegranate
                    
RUN micromamba clean -a

ENV PATH=/opt/conda/envs/vcf-kit/bin:${PATH}
RUN micromamba env export --name vcf-kit > vcf-kit.yml

RUN pip install https://github.com/AndersenLab/VCF-kit/archive/refs/tags/0.3.0.tar.gz

LABEL Name="vcf-kit" Author="Daniel Cook"
USER root
RUN apk add --no-interactive procps
