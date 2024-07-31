FROM continuumio/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --add channels danielecook \
    && conda create -n vcf-kit \
                       danielecook::vcf-kit=0.2.6 \
                       "bwa>=0.7.17" \
                       "samtools>=1.10" \
                       "bcftools>=1.10" \
                       "blast>=2.2.31" \
                       "muscle>=3.8.31" \
                       "primer3>=2.5.0" \
    && conda clean -a

ENV PATH /opt/conda/envs/vcf-kit/bin:${PATH}
RUN conda env export --name vcf-kit > vcf-kit.yml
LABEL Name="vcf-kit" Author="Daniel Cook"
RUN apt-get update && apt-get install -y procps
