FROM continuumio/miniconda3:4.8.2

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --add channels danielecook \
    && conda create -n vcf-kit danielecook::vcf-kit=2.6 \
    && conda clean -a

ENV PATH /opt/conda/envs/vcf-kit/bin:${PATH}
RUN conda env export --name vcf-kit > vcf-kit.yml
LABEL Name="vcf-kit" Author="Daniel Cook"
