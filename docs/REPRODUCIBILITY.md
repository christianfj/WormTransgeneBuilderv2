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
# Intel/AMD computers, Windows and most servers (amd64); use "linux-arm64" for Apple Silicon Macs (M1 and later):
shasum -a 256 -c wormtransgenebuilder-v2.0.1-linux-amd64.tar.gz.sha256      # check the download (on Linux: sha256sum -c ...)
docker load < wormtransgenebuilder-v2.0.1-linux-amd64.tar.gz
docker run --rm -p 3838:3838 ghcr.io/christianfj/wormtransgenebuilderv2:v2.0.1
```

Open <http://localhost:3838>. From version 2.0.1 on, each release has **two** images: `amd64` (Intel/AMD processors, Windows and
most servers) and `arm64` (Apple Silicon Macs and arm64 servers). The name `ghcr.io/christianfj/wormtransgenebuilderv2:<version>` is
one multi-architecture image: `docker pull` and `docker run` choose the right one for your computer automatically. The saved files
(`...-linux-amd64.tar.gz` and `...-linux-arm64.tar.gz`) are for loading without a registry; use the one that matches your computer
(`uname -m` prints `x86_64` for amd64 and `arm64` for Apple Silicon). Version 2.0.0 has an amd64 image only; on Apple Silicon Docker
runs it by emulation, which is slower. Each release also has a small file `wormtransgenebuilder-<version>-image-digest.txt` with the
exact digests of its images.

## Rebuilding the image

```sh
docker build -t wormtransgenebuilder .
```

This takes about 10 minutes on an amd64 machine. It fails if `RNAfold` does not return the expected energy
for a test sequence.

## Running the tests

```r
renv::restore()        # exact packages
devtools::test()       # about 4,180 checks, about 3 minutes
```

`RNAfold` (ViennaRNA 2.7.2) is needed for the ribosome-binding-site tests, which skip themselves otherwise.
One test drives the app in headless Chrome and is skipped if Chrome is absent.

## What results can be reproduced exactly

* **High expression** and all restriction-site/piRNA choices under it: no random numbers; the same input always gives the same sequence. (With the ribosome-binding-site step the start is drawn at random, so the app fixes the seed and the result is still repeatable.)
* **No codon optimization:** fixed seed; the same input gives the same output.
* **Ubiquitous, Random, GLO:** random choices; each run differs. Tests therefore check properties (same protein, sites removed, thresholds met) rather than exact sequences. For a fixed seed in a script, `set.seed()` before calling `optimize_transgene()` gives repeatable results with the same software versions.

During development the speed-ups (piRNA matching, GLO scoring, restriction-site search) were each checked
to give outputs identical to the earlier implementations, and the same fixed run gave identical output on
macOS (arm64), Linux arm64 and Linux amd64.
