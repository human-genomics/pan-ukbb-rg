# =============================================================================
# Pan-UKBB EUR One-vs-All Genetic Correlations - Docker Image
# =============================================================================
# Builds a reproducible Linux x86_64 image with the repo scripts, a pinned
# Neale LDSC checkout, and the Python 2.7 LDSC conda environment pre-installed.
# Public Pan-UKBB data downloads and generated sumstats happen at runtime.
#
# Build:
#   docker build -t pan-ukbb-rg .
#
# Run:
#   docker run --rm -v $(pwd)/pan-ukbb-rg-work:/app/pipeline-output pan-ukbb-rg help
# =============================================================================

FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/human-genomics/pan-ukbb-rg"
LABEL org.opencontainers.image.description="Pan-UKBB EUR LDSC one-vs-all genetic correlation pipeline"
LABEL org.opencontainers.image.licenses="MIT"

ARG MAMBA_VERSION=1.5.10
ARG LDSC_REPO=https://github.com/astheeggeggs/ldsc.git
ARG LDSC_COMMIT=a4ee4c8aa065a1c9a586c3b678e9b3040bbebafc

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV LDSC_ENV_PREFIX=/opt/ldsc-neale
ENV LDSC_PYTHON=/opt/ldsc-neale/bin/python
ENV LDSC_DIR=/opt/ldsc-neale-src
ENV ENV_MANAGER=micromamba
ENV PYTHON=python3
ENV JOBS=16

RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        bzip2 \
        ca-certificates \
        coreutils \
        curl \
        findutils \
        gawk \
        git \
        gzip \
        make \
        procps \
        python3 \
        tar \
        util-linux \
    && rm -rf /var/lib/apt/lists/*

RUN curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}" \
    | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

WORKDIR /app

COPY envs/ldsc-neale.yml envs/ldsc-neale.yml
COPY scripts/check_ldsc_env.sh scripts/check_ldsc_env.sh

RUN micromamba create -y -p "${LDSC_ENV_PREFIX}" -f envs/ldsc-neale.yml \
    && micromamba clean --all --yes \
    && bash scripts/check_ldsc_env.sh "${LDSC_PYTHON}"

RUN git clone "${LDSC_REPO}" "${LDSC_DIR}" \
    && git -C "${LDSC_DIR}" checkout "${LDSC_COMMIT}"

COPY . .

RUN chmod +x scripts/docker_entrypoint.sh scripts/check_ldsc_env.sh scripts/setup_neale_ldsc.sh \
    && python3 -m py_compile scripts/*.py \
    && bash scripts/check_ldsc_env.sh "${LDSC_PYTHON}" \
    && git -C "${LDSC_DIR}" rev-parse --verify "${LDSC_COMMIT}^{commit}" >/dev/null

VOLUME ["/app/pipeline-output"]

ENTRYPOINT ["scripts/docker_entrypoint.sh"]
CMD ["help"]
