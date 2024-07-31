FROM continuumio/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --add channels danielecook \
    && conda create -n "bwa>=0.7.17" \
                       "samtools>=1.10" \
                       "bcftools>=1.10" \
                       "blast>=2.2.31" \
                       "muscle>=3.8.31" \
                       "primer3>=2.5.0" \
    && conda clean -a

RUN pip install https://github.com/AndersenLab/VCF-kit/archive/0.3.0.tar.gz

LABEL Name="vcf-kit" Author="Daniel Cook"
RUN apt-get update && apt-get install -y procps
