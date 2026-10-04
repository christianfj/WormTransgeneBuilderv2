# Benchmarks

How long each step of WormTransgeneBuilder takes, and how much memory it uses, on fixed
test genes. The point is to tell whether a change made the app faster, and to notice if it
also changed the output.

## Run it

Needs the `bench`, `devtools` and `rlang` packages (`install.packages("bench")`).

```sh
Rscript dev/benchmarks/benchmark.R                      # standard: 300 bp and 3 kb genes, about 15-25 minutes
Rscript dev/benchmarks/benchmark.R --suite quick        # 300 bp only, a few minutes
Rscript dev/benchmarks/benchmark.R --suite full         # adds a 20 kb gene (the app's limit); slow
Rscript dev/benchmarks/benchmark.R --label before-x     # name the results file
Rscript dev/benchmarks/benchmark.R --only "pirna"       # only cases whose name matches
Rscript dev/benchmarks/benchmark.R --iterations 10      # repeat each case more often
```

Results are written to `dev/benchmarks/results/bench_<label>.csv`.

## Compare two runs

```sh
Rscript dev/benchmarks/compare.R dev/benchmarks/results/bench_before.csv dev/benchmarks/results/bench_after.csv
```

Shows the median time of each case in both runs, the speed-up (before / after, so `2.00x`
means twice as fast), the change in memory, and whether the output changed.

## Workflow for a speed-up

1. Before changing anything, run the benchmark with a clear label: `--label before-<change>`.
2. Make the change and run the tests (`devtools::test()`).
3. Run the benchmark again: `--label after-<change>`.
4. Compare. Every row should say `same` in the `output` column. `CHANGED` means the code now
   returns something different under the same random seed. That is expected only when you
   changed an algorithm on purpose (and the tests should say so); for a pure speed-up it is
   a red flag.

## What is measured

- **Start-up**, in fresh R sessions: `devtools::load_all()`, and the installed package (the
  way the app runs in Docker): the `library()` call, and the first use of the 54 MB GLO
  dictionary, which is loaded lazily.
- **Core steps** at each gene size: GC, CAI and GLO scores, synthesis checks, codon recoding,
  GLO optimization, piRNA distances and both piRNA shields, restriction-site removal (one
  enzyme and all 21, random and High-expression methods), intron insertion, RBS optimization,
  the sequence viewer and the GenBank export.
- **The whole pipeline**, `optimize_transgene()`, for each method, and once with every option on.

Each case runs under a fixed random seed. The output column is a hash of what the code
returned, so identical output gives an identical hash.

## The test genes

`make_benchmark_cds()` builds a reproducible gene of a given length (300 bp, 3 kb, 20 kb):
random codons drawn with the Ubiq usage weights, with real piRNA targets (about one per
900 bp) and BsaI / XhoI sites (about one per 1500 bp) planted in frame, so the shields have
real work to do. `make_best_codon_benchmark_cds()` does the same with every codon set to the
best codon, for the High-expression shields, which start from a recoded sequence.

## Reading the numbers

- Timings depend on the computer and on what else it is doing. Only compare runs made on the
  same machine, with nothing else busy. `compare.R` warns if the platform, R version or number
  of cores differ.
- `git_state` in the results says whether the code had uncommitted changes when it was run.
- The piRNA steps use several CPU cores (`stringdist` with all cores but one), so they are
  sensitive to other programs running at the same time.
- 20 kb runs take a long time; the `full` suite repeats them only once.
- Whole-pipeline timings depend on how much work the shields find. Recoding can destroy a
  planted piRNA target or restriction site by chance, so a method that recodes more can look
  faster (for example GLO, whose recoding happens to remove the planted target in the 300 bp
  gene). That is real behaviour, not noise, and it is the same every run because the seed is
  fixed, but read the pipeline rows together with the single-step rows.
- `--outdir` writes the results file somewhere other than `dev/benchmarks/results/`.

## Known cost centres (first baseline, 300 bp gene, Apple Silicon)

Measured before any optimization, for orientation only. See the results files for the
current numbers.

- `get_pirna_distances` takes about 1.2 s on a 300 bp gene, because every 20-mer window is
  compared with every piRNA in the database (18,150 when this was measured; 17,849 since the
  301 entries without a sequence were removed). The whole pipeline calls it several times.
- The random piRNA shield takes about 5 s and the High-expression one about 2 s.
- High-expression restriction-site removal with all 21 enzymes takes about 2.7 s.
- On an installed package, the first use of the old 54 MB GLO dictionary took about 9 s (now
  0.23 s with the compact score table; the start-up row of the benchmark times the new loader).

## Comparing one function old against new

Whole-suite timings are noisy when other programs use the CPU (a VPN extension alone can use
half a core). To compare one function reliably, alternate the old and new code, repeat, and
take the best time. Two scripts do this and also check the results are exactly identical:

- `compare_pirna_distances.R`: the original `stringdist` comparison against the bit-encoded,
  multi-core piRNA engine.
- `compare_glo.R`: the previous `calc_glo_score()` / `optimize_glo()` (taken from git history)
  against the score-table versions, including the memory the score data uses.
