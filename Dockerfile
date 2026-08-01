FROM ubuntu:24.04 AS builder

ARG BUILD_JOBS=8
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        cmake \
        libbz2-dev \
        python3 \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN PREFIX=/opt/spades ./spades_compile.sh \
        -j "${BUILD_JOBS}" \
        -DSPADES_ENABLE_PROJECTS=release \
        -DSPADES_USE_NCBISDK=ON


FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="SPAdes" \
      org.opencontainers.image.description="Genome assembler for isolate, single-cell, metagenomic, and transcriptomic sequencing data" \
      org.opencontainers.image.url="https://ablab.github.io/spades/" \
      org.opencontainers.image.source="https://github.com/ablab/spades" \
      org.opencontainers.image.licenses="GPL-2.0-only"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        libbz2-1.0 \
        libgomp1 \
        pigz \
        python3 \
        zlib1g \
    && groupadd --gid 1000 spades \
    && useradd --uid 1000 --gid spades --create-home --shell /bin/bash spades \
    && install --directory --owner spades --group spades /work \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/spades /opt/spades

ENV PATH="/opt/spades/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

WORKDIR /work
USER spades

CMD ["/bin/bash"]
