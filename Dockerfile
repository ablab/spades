FROM python:3.10.14-slim-bullseye

LABEL    software="spades" \ 
    author="Sebastian Bassi <sebastian@toyoko.io>" \
    base_image="python:3.10.14-slim-bullseye" \ 
    container="spades" \ 
    about.summary="genome assembler for single-cell and isolates data sets" \ 
    about.home="https://ablab.github.io/spades/" \ 
    software.version="4.0.0" \ 
    extra.identifiers.biotools="spades" \ 
    about.paper="Prjibelski A, Antipov D, Meleshko D, Lapidus A, Korobeynikov A (2020) Using SPAdes de novo assembler. Current Protocols in Bioinformatics 70(1): e102. https://doi.org/10.1002/cpbi.102" \ 
    about.license="GPL-2+" \ 
    about.tags="biology::nucleic-acids, field::biology, field::biology:bioinformatics,:c++,:python, interface::commandline,:program, scope::utility, works-with::biological-sequence,:file" 

USER root

RUN apt-get update && apt-get install -y build-essential cmake wget libbz2-dev

RUN wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0.tar.gz

RUN tar xfz SPAdes-4.0.0.tar.gz && rm SPAdes-4.0.0.tar.gz && cd SPAdes-4.0.0/ && ./spades_compile.sh

# ./SPAdes-4.0.0/bin/spades.py --test
