# syntax=docker/dockerfile:1

# ---- build stage: compile/install heavy deps ----
FROM rocker/r-ver:4.3.3 AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        cmake \
    && rm -rf /var/lib/apt/lists/*

# Posit Package Manager binary snapshot pinned to a date near R 4.3.3. rocker/r-ver:4.3.3 is built
# on Ubuntu 22.04 (jammy), NOT the WSL host's 24.04 (noble) -- verified via /etc/os-release inside
# the image; a noble-targeted binary repo silently installs libstdc++-incompatible binaries
# (GLIBCXX version mismatch at dyn.load time) that fail only when a package actually loads.
# renv.lock records this same repo (minus the platform suffix, which renv re-derives from the
# platform it's restoring on), so this only needs to match at first-snapshot time, not every build.
ENV CRAN_REPO=https://packagemanager.posit.co/cran/__linux__/jammy/2024-06-01
ENV RENV_CONFIG_REPOS_OVERRIDE=${CRAN_REPO}
# renv normally installs into a global cache and leaves symlinks in renv/library -- those symlinks
# would point outside this build stage's filesystem and break when copied into the runtime stage.
# Force real copies into renv/library so it's self-contained and safe to COPY across stages.
ENV RENV_CONFIG_CACHE_SYMLINKS=FALSE

WORKDIR /project

# renv.lock is the source of truth for exact package versions (see dev/bootstrap_renv.R for how it
# was generated). Only copy what renv::restore() needs -- not the whole repo -- so R source
# changes don't invalidate this layer's cache.
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/settings.json renv/
RUN R -e "install.packages('renv', repos='${CRAN_REPO}')"
RUN R -e "renv::restore()"

# ---- runtime stage: slim image, non-root ----
FROM rocker/r-ver:4.3.3 AS runtime

# Quarto CLI (for analysis/explore.qmd) is pinned by version + verified download URL, the same
# reproducibility standard as everything else in this image -- not "curl the latest installer".
ENV QUARTO_VERSION=1.10.18

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 libssl3 libxml2 libfontconfig1 libharfbuzz0b libfribidi0 \
        libfreetype6 libpng16-16 libtiff5 libjpeg-turbo8 \
        curl ca-certificates \
    && curl -fsSL -o /tmp/quarto.deb \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && rm /tmp/quarto.deb \
    && apt-get purge -y curl && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 1000 sim

WORKDIR /project
COPY --from=build --chown=sim:sim /project/renv ./renv
COPY --chown=sim:sim . .

USER sim
ENTRYPOINT ["Rscript"]
CMD ["run.R"]
