# Reproducibility

## What is fixed

| Part | How it is pinned |
|---|---|
| Operating system and R | The base image `rocker/r-ver:4.5.3` (R 4.5.3, Ubuntu 24.04), by its content digest (see the first lines of the `Dockerfile`) |
| R packages | `renv.lock` records the exact version of every package; the image installs them from a dated snapshot of CRAN (the `CRAN_SNAPSHOT` line in the `Dockerfile`) |
| ViennaRNA | Version 2.7.2, compiled from the source archive, checked against its SHA-256 checksum |
| The app | The release tag (for example `v2.0.0`) |

The version of the app, R and ViennaRNA is printed when the image starts and shown in the **About** tab.
Quote them with any design.

## Running an archived version

Each release carries a saved copy of its image (with a checksum) as a file attached to the GitHub release. Zenodo's automatic GitHub archiving stores only the *source code* of a release, so the saved image is uploaded to Zenodo separately, by hand, and listed in the same record or linked from it (see the DOI in the README).

```sh
sha256sum -c wormtransgenebuilder-v2.0.0-linux-amd64.tar.gz.sha256   # check the download
docker load < wormtransgenebuilder-v2.0.0-linux-amd64.tar.gz
docker run --rm -p 3838:3838 ghcr.io/christianfj/wormtransgenebuilderv2:v2.0.0
```

Open <http://localhost:3838>. The image is built for Intel/AMD processors (amd64); on Apple Silicon Docker
runs it by emulation, or you can rebuild it (below).

## Rebuilding the image

```sh
docker build -t wormtransgenebuilder .
```

This takes about 10 minutes on an amd64 machine. It fails if `RNAfold` does not return the expected energy
for a test sequence.

## Running the tests

```r
renv::restore()        # exact packages
devtools::test()       # about 3,650 checks, about 3 minutes
```

`RNAfold` (ViennaRNA 2.7.2) is needed for the ribosome-binding-site tests, which skip themselves otherwise.
One test drives the app in headless Chrome and is skipped if Chrome is absent.

## What results can be reproduced exactly

* **High expression** and all restriction-site/piRNA choices under it: no random numbers; the same input always gives the same sequence.
* **No codon optimization:** fixed seed; the same input gives the same output.
* **Ubiquitous, Random, GLO:** random choices; each run differs. Tests therefore check properties (same protein, sites removed, thresholds met) rather than exact sequences. For a fixed seed in a script, `set.seed()` before calling `optimize_transgene()` gives repeatable results with the same software versions.

During development the speed-ups (piRNA matching, GLO scoring, restriction-site search) were each checked
to give outputs identical to the earlier implementations, and the same fixed run gave identical output on
macOS (arm64), Linux arm64 and Linux amd64.
