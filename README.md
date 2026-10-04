# WormTransgeneBuilder

An online app for comprehensive *C. elegans* transgene design.

> **Beta release.** The tool is being tested by the *C. elegans* community, and the manuscript describing it
> is in preparation (not yet submitted). Results are meant to be checked by the user before ordering
> synthesis. Please report problems as [GitHub issues](../../issues). The app is expected to be maintained
> only until about mid-2027; the archived release (see [Citing](#citing-and-archived-versions)) will stay
> usable after that.

**Use it online:** <https://www.wormbuilder.org/transgenebuilder/>

## What it does

Paste a coding sequence (DNA or protein, up to 20 kb) and choose what to do with it:

| Step | Options |
|---|---|
| **Codon optimization** | Ubiquitous expression · High expression · Germline-optimized (GLO) · Random recoding (synonymous codons equally likely) · none |
| **Molecular tags** | N- and C-terminal tags (for example SV40 NLS, *egl-13* NLS, linkers) |
| **Biological shields** | Minimize piRNA homology (choose how many mismatches to known piRNAs are required) · remove restriction-enzyme sites (234 NEB enzymes, grouped as *commonly used*, *Golden Gate* and *all other*) |
| **Gene architecture** | Ribosome-binding-site optimization (ViennaRNA) · consensus start (aaaaATG) · native or synthetic introns · 5' and 3' UTRs · promoters (38 validated 300 bp promoters, grouped by neuron and tissue) |
| **Synthesis checks** | Twist Bioscience guidelines (homopolymers, global and local GC content) |

The result is shown in a sequence viewer (promoter, UTRs, tags, introns, piRNA and enzyme sites, non-optimal
codons) and can be downloaded as an annotated GenBank file (opens in ApE, Benchling and similar programs).
The amino-acid sequence is never changed by the optimization steps. How each step works, with references,
is in [docs/METHODS.md](docs/METHODS.md).

## Run it yourself

### With Docker (recommended: the same software as the published version)

```sh
docker run --rm -p 3838:3838 ghcr.io/christianfj/wormtransgenebuilderv2:latest
```

Then open <http://localhost:3838>. For a specific, archived version use its tag instead of `latest`
(for example `:v2.0.0`). A saved copy of every release is also attached to the release on GitHub (and
deposited on Zenodo); see [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

### From R

You need R 4.5, and, for the ribosome-binding-site step, [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)
(the program `RNAfold` on your `PATH`; without it that step is skipped with a warning).

```r
# in the project folder
install.packages("renv")
renv::restore()          # installs the exact package versions recorded in renv.lock
devtools::load_all()
run_app()
```

## Documentation

* [docs/METHODS.md](docs/METHODS.md): what each step does and the references behind it.
* [docs/DATA.md](docs/DATA.md): the data tables the app uses and where they come from.
* [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md): exact versions, how to run an archived image, how to run the tests.
* [docs/LIMITATIONS.md](docs/LIMITATIONS.md): what the tool does not do, and known limits.
* [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md): running the image as a public service.

## Testing

`devtools::test()` runs about 3,650 automated checks (the amino-acid sequence is unchanged, requested
restriction sites are gone, piRNA matches stay under the threshold, GenBank files are well formed, and
so on). They also run on GitHub on every change.

## Citing and archived versions

A manuscript is in preparation. Until it is available, please cite the software as in
[CITATION.cff](CITATION.cff) (GitHub shows a "Cite this repository" button). Every tagged release is
archived on Zenodo with a DOI: *DOI to be added after the first release*.

## Authors and licence

Amhed M. Vargas-Velazquez, Sonia El Mouridi, Sarah AlHarbi, Khlifa Alnaim, Henrik Bringmann, Daniel J.
Dickinson and Christian Frøkjær-Jensen (corresponding author: cfjensen@kaust.edu.sa), Laboratory of
Synthetic Genome Biology, King Abdullah University of Science and Technology (KAUST), with the Technical
University Dresden and the University of Texas at Austin. Released under the
[GNU General Public License v3](LICENSE). The data tables have their own sources; see
[docs/DATA.md](docs/DATA.md).
