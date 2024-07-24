FROM python:3.10.14-slim-bullseye

LABEL    software="spades" \ 
    author="Sebastian Bassi <sebastian@toyoko.io>" \
    base_image="python:3.10.14-slim-bullseye" \ 
    container="spades" \ 
    about.summary="genome assembler for single-cell and isolates data sets" \ 
    about.home="https://ablab.github.io/spades/" \ 
    software.version="4.0.0" \ 
    upstream.version="4.0.0" \ 
    version="1" \ 
    extra.identifiers.biotools="spades" \ 
    about.copyright="2011-2014 Saint-Petersburg Academic University" \ 
    about.license="GPL-2+" \ 
    about.license_file="/usr/share/doc/spades/copyright" \ 
    extra.binaries="/usr/bin/metaspades,/usr/bin/metaspades.py,/usr/bin/plasmidspades,/usr/bin/plasmidspades.py,/usr/bin/rnaspades,/usr/bin/rnaspades.py,/usr/bin/spades,/usr/bin/spades.py,/usr/bin/truspades,/usr/bin/truspades.py" \ 
    about.tags="biology::nucleic-acids, field::biology, field::biology:bioinformatics,:c++,:python, interface::commandline,:program, scope::utility, works-with::biological-sequence,:file" 

USER root

RUN apt-get update && apt-get install -y build-essential cmake wget libbz2-dev

RUN wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0.tar.gz

RUN tar xfz SPAdes-4.0.0.tar.gz && rm SPAdes-4.0.0.tar.gz && cd SPAdes-4.0.0/ && ./spades_compile.sh

# ./SPAdes-4.0.0/bin/spades.py --test
