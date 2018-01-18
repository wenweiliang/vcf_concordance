FROM biocontainers/biocontainers:latest

LABEL base.image="biocontainers:latest"
LABEL description="Comparison of different VCF files"
LABEL tags="Genomics"

MAINTAINER Wen-Wei Liang <wenwiliang@gmail.com>

USER root

RUN conda install pyvcf tabix
