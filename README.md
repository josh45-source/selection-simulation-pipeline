# Selection simulation MVP

A long-horizon selection simulation pipeline built on AlphaSimR, with recurrent mutation,
checkpoint/resume, and a Parquet-plus-DuckDB storage layer queried through a Quarto dashboard.
The design goal is a lean MVP: the smallest infrastructure that proves the full loop works and
gets the scientific details right, not a scaled-up production system. See
`claude_code_prompt_selection_mvp.md` for the original design brief.

**Status.** The full loop runs end to end: simulate, write Parquet, query with DuckDB, render a
dashboard. Checkpoint/resume and independent RNG streams per replicate are implemented and
verified bit-for-bit. Four validation checks are in place and passing: neutral heterozygosity
decay, the breeder's equation, F_ROH against pedigree F, and inbreeding depression against a
literature benchmark. Read the Limitations section before trusting numbers from a run: two of the
project's central scientific claims (reproducibility, inbreeding depression) come with real
caveats that are easy to miss if you only skim the validation suite's PASS output.

## Setup

Everything runs in Docker; there is no supported way to run this outside a container, since the
package versions are pinned through `renv.lock` and the image is the unit of reproducibility. You
need Docker and enough disk space for the image (a few GB, mostly `arrow` and `duckdb`).

```bash
docker build -t selection:latest .
```

The build has two stages. The first restores the exact package set from `renv.lock` using
`renv::restore()`, pointed at a dated Posit Package Manager snapshot rather than "whatever CRAN
has today." The second copies the restored library and the project source into a slim,
non-root runtime image. Rebuilding after a source change is fast: only the second stage reruns,
since Docker caches the package-restore layer as long as `renv.lock` hasn't changed.

## Running a scenario

The image runs as a non-root user, so the bind-mounted `data/` output directory must be owned by
whatever UID actually runs the container. Otherwise Docker auto-creates it as root on first mount
and the container can't write to it. `--user "$(id -u):$(id -g)"` runs the container as your own
host UID and GID instead of the image's baked-in default, so this works regardless of what UID
you are on your machine, not just uid 1000:

```bash
mkdir -p data
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd)/data:/project/data" \
  selection:latest run.R config/scenario.slice.yaml   # fast smoke test, seconds
  # or: config/scenario.example.yaml for the full-scale scenario, much longer
```

`run.R` runs every replicate of every selection scenario named in the config, writes the results
to partitioned Parquet under `data/`, and produces a summary plot. Re-running the same command
after an interruption resumes automatically: a checkpoint left on disk from a prior partial run
of that exact scenario, replicate, and selection intensity is picked up rather than restarting
from generation 1. See `R/checkpoint.R` and `R/sim_loop.R` for how checkpointing works, and
`validation/test_resume_equals_uninterrupted.R` for the proof that a resumed run is bit-identical
to an uninterrupted one.

## Rendering the dashboard

The dashboard reads whatever is currently in `data/`, so run a scenario first. Quarto needs its
own entrypoint, since the image's default entrypoint is `Rscript`:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --entrypoint quarto \
  -v "$(pwd)/data:/project/data" \
  -v "$(pwd)/analysis:/project/analysis" \
  selection:latest render analysis/explore.qmd
```

This writes `analysis/explore.html`, self-contained with embedded images, ready to open directly
in a browser. It shows additive genetic variance and fixation progress as distributions across
replicates rather than single trajectories, since a lone run of a stochastic process is close to
uninformative on its own, and a validation summary with citations resolved against
`references.bib`.

## Running the validation suite

Each validation script is self-contained and can be run directly against the image:

```bash
docker run --rm --user "$(id -u):$(id -g)" selection:latest validation/test_froh.R
```

Swap in any of `test_neutral_heterozygosity.R`, `test_breeders_equation.R`,
`test_replicate_independence.R`, `test_resume_equals_uninterrupted.R`, or
`test_inbreeding_depression.R`. Each one prints the actual numbers it produced, not just a
pass/fail line: realized versus nominal effective population size, where the breeder's equation
starts overpredicting, the correlation between F_ROH and pedigree F, and the inbreeding
depression slope. `dev/` holds the scratch scripts that found and diagnosed the bugs documented
below; `dev/README.md` explains which ones are load-bearing evidence and which were one-off
debugging.

## Data model

Three Parquet datasets, partitioned by `scenario_id` and `replicate_id`, under `data/`:

- `summary_metrics`: narrow, one row per generation and replicate. Additive genetic variance and
  mean genetic value.
- `per_locus_freq`: per-locus allele frequency, sampled every `sampling.locus_cadence`
  generations and restricted to segregating sites by default, to keep the table small rather than
  writing one row per locus per generation per replicate.
- `founders/<scenario>_rep<n>.rds`: not Parquet, and not yet a formal provenance table. One
  serialized founder population per replicate, written once, immediately after it's built. See
  Limitations for why this exists.

`data/checkpoints/` holds the rolling checkpoint used for resume; unlike the founders directory,
these files are overwritten as a run progresses and aren't meant to be kept as a historical
record.

## Project structure

```
.
├── Dockerfile                  # multi-stage build; renv::restore(), then a slim runtime image
├── renv.lock                   # pinned package versions, the source of truth for reproducibility
├── references.bib              # bibliography, cited from code comments and the dashboard
├── config/
│   ├── scenario.example.yaml   # full-scale scenario
│   └── scenario.slice.yaml     # small, fast scenario for smoke-testing the pipeline
├── R/
│   ├── engine.R                # founder population and genome setup
│   ├── mutation.R              # recurrent mutation via AlphaSimR's mutate()
│   ├── sim_loop.R              # the serial per-generation loop, with checkpointing
│   ├── checkpoint.R            # save/restore pop, SimParam, RNG state, and accumulated output
│   ├── rng.R                   # independent L'Ecuyer-CMRG streams per replicate
│   ├── metrics.R               # additive variance, allele frequency, fixation, heterozygosity
│   ├── froh.R                  # F_ROH via detectRUNS
│   ├── pedigree.R              # pedigree-based F via a direct tabular-method implementation
│   └── io.R                    # partitioned Parquet writer
├── run.R                       # entrypoint: reads a scenario config, runs it, writes output
├── validation/                 # the five validation checks, each self-contained and runnable alone
├── analysis/
│   └── explore.qmd             # the dashboard
└── dev/                        # diagnostic scripts that found real bugs; see dev/README.md
```

## Limitations

**Founder populations are not reproducible from a seed.** `AlphaSimR:::MaCS()`, the compiled
coalescent function that `runMacs()` calls to build founder haplotypes, does not produce
identical output across calls given identical literal seeds. This was confirmed directly,
bypassing R's own RNG entirely (`dev/verify_macs_direct.R`), and it is a limitation in AlphaSimR's
compiled code, not something fixable from this project. It is also undocumented: AlphaSimR's
changelog shows the maintainers forced single-threading in their own doc examples (v1.5.1), which
independently confirms that multi-threaded execution isn't reproducible, but nowhere describes
fixing MaCS's seed handling itself. Everything *after* founders exist, including mutation,
selection, crossing, and checkpoint/resume, is reproducible and proven bit-identical. The founder
population itself is not. State this plainly: **the qualitative pattern of a run is robust, but
the exact numbers for a specific run come from the persisted founder object, not the seed.** A
scenario config and a seed will get you a run with the same statistical character, not the same
run. `run.R` persists each replicate's founder population once, to
`data/founders/<scenario>_rep<n>.rds`, specifically so a given run's founders can be recovered
exactly. A future Provenance dataset (git SHA, container digest, full parameter set) should
reference this file, not the seed, as the record of how to reproduce a run.

**Inbreeding depression requires a trait architecture the main pipeline doesn't use.** A purely
additive trait, which is what every scenario config in this repository uses, cannot show mean
inbreeding depression by construction: depression requires directional dominance, where the
higher-value allele tends to be dominant, so rising homozygosity from inbreeding unmasks more of
the lower-value allele's effect. `validation/test_inbreeding_depression.R` builds its own
directional-dominance population with `SP$addTraitAD` specifically to demonstrate this mechanism,
separately from `build_founder_pop()` and every other scenario in the project. Extending dominance
to the main pipeline, so that scenarios themselves show inbreeding depression rather than only the
validation test, is a real feature and has not been decided or built.

**The Doekes benchmark comparison is qualitative, not a literal reproduction.** The simulated
trait has an arbitrary, mean-zero scale rather than a natural biological unit like birth weight or
litter size, so a percent-of-mean comparison to the Doekes et al. (2021) meta-analysis benchmark
of roughly 0.13% decline per 1% rise in F is not well defined here: dividing by a founder-generation
mean phenotype close to zero makes the ratio numerically unstable, which was confirmed
empirically before this was caught. `test_inbreeding_depression.R` instead reports the slope in
phenotypic standard deviation units per unit F, which is stable and directionally comparable to
Doekes' finding of small but real depression, without claiming numeric equivalence to it.

**`optiSel` is not installed.** It's unresolvable from the pinned package repository for R 4.3.3,
for reasons unrelated to its own version floor (which is only R >= 3.5.0). `detectRUNS` covers
F_ROH directly and is what this project actually uses; `optiSel` would only matter if optimal
contribution selection features are added later.

**MoBPS remains the documented fallback**, per the original design brief, if AlphaSimR's founder
non-determinism or its other constraints become limiting, particularly for recurrent mutation and
optimal contribution selection.
