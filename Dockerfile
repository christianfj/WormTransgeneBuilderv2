# WormTransgeneBuilder: reproducible image.
#
# Everything that could change between two builds is pinned:
#   * the base image (R 4.5.3 on Ubuntu 24.04), by its content digest;
#   * ViennaRNA, by version and by the checksum of its source code;
#   * every R package, by renv.lock, installed from a dated snapshot of CRAN.
# Build:  docker build -t wormtransgenebuilder .
# Run:    docker run --rm -p 3838:3838 wormtransgenebuilder     (then open http://localhost:3838)

ARG R_IMAGE=rocker/r-ver:4.5.3@sha256:c3f39b365d1077fe24f8e9ab2742e352b6d3950897f51af1624a5bb5550c21c0


# ---- Stage 1: compile ViennaRNA (provides the RNAfold program used for the ribosome binding site) ----
FROM ${R_IMAGE} AS viennarna
ARG VIENNARNA_VERSION=2.7.2
ARG VIENNARNA_SHA256=1ab5f4a4f76fc85a2243546088e45f5d85f2d7a56cc656e969b005cce9bfab5f
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates pkg-config \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL -o /tmp/ViennaRNA.tar.gz \
      "https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-${VIENNARNA_VERSION}.tar.gz" \
 && echo "${VIENNARNA_SHA256}  /tmp/ViennaRNA.tar.gz" | sha256sum -c - \
 && tar -xzf /tmp/ViennaRNA.tar.gz -C /tmp \
 && cd "/tmp/ViennaRNA-${VIENNARNA_VERSION}" \
 && ./configure --prefix=/opt/viennarna --without-perl --without-python --without-forester \
      --without-kinfold --without-rnalocmin --without-doc --disable-static --disable-lto \
 && make -j"$(nproc)" \
 && make install


# ---- Stage 2: the app ----
FROM ${R_IMAGE}
ARG CRAN_SNAPSHOT=https://p3m.dev/cran/__linux__/noble/2026-09-30

# ViennaRNA's run-time library (libgomp1), libuv (used by the R package fs), and curl for the health check.
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 libuv1t64 curl \
 && rm -rf /var/lib/apt/lists/*
# RNAfold is a single program (its library is built into it), so only that file is copied.
COPY --from=viennarna /opt/viennarna/bin/RNAfold /usr/local/bin/RNAfold

# R packages, at the versions recorded in renv.lock (only those the app needs at run time, with what
# they depend on; the test and development packages in the lock file are not installed).
ENV RENV_CONFIG_REPOS_OVERRIDE=${CRAN_SNAPSHOT}
COPY renv.lock /build/renv.lock
RUN R -q -e "options(repos = c(CRAN = Sys.getenv('RENV_CONFIG_REPOS_OVERRIDE'))); install.packages('renv')" \
 && R -q -e "renv::restore(lockfile = '/build/renv.lock', library = .Library.site[1], prompt = FALSE, \
      packages = c('bslib', 'config', 'golem', 'shiny', 'shinyjs', 'shinyWidgets', 'stringdist'))" \
 && rm -rf /build /tmp/*

# The app itself (only the files that make up the package).
COPY DESCRIPTION NAMESPACE /app/
COPY R /app/R
COPY data /app/data
COPY inst /app/inst
COPY man /app/man
RUN R CMD INSTALL --no-test-load /app && rm -rf /app
COPY docker/start.R /usr/local/share/wormtransgenebuilder/start.R

# Build-time checks: RNAfold must give the known result for a known sequence, and the app must load.
RUN test "$(echo GGGGAAAACCCC | RNAfold --noPS | tail -1)" = "((((....)))) ( -5.40)" \
 && R -q -e "library(WormTransgeneBuilder); stopifnot(nzchar(system.file('extdata', 'Promoters.csv', package = 'WormTransgeneBuilder')))"

# Run as an unprivileged user.
RUN useradd --create-home --uid 1000 app
USER app
WORKDIR /home/app

EXPOSE 3838
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD curl -fsS http://localhost:3838/ > /dev/null || exit 1

CMD ["Rscript", "/usr/local/share/wormtransgenebuilder/start.R"]
