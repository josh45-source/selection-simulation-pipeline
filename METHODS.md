# Methods

This document consolidates the simulation model, the reproducibility design, and the validation
results into one narrative. It draws only on work already in the repository: the R source, the
scenario config, and the validation suite. Where a number appears below, it comes from an actual
run of the corresponding validation test, not a projection.

## Simulation model

The engine is AlphaSimR [@gaynor2021alphasimr], run through a thin wrapper (`R/engine.R`,
`R/mutation.R`, `R/sim_loop.R`) rather than called directly, so the engine can be swapped later
without touching the rest of the pipeline. MoBPS [@pook2020mobps] is the documented fallback if
AlphaSimR's constraints, particularly around recurrent mutation, become limiting.

**Genome.** Each scenario defines three site classes per chromosome, all designated at founder
initialization from an oversampled coalescent sample (`AlphaSimR::runMacs`, `species = "GENERIC"`):
a QTL panel that carries the trait, a mutating neutral-marker panel used for general
allele-frequency and fixation tracking, and a validation-marker panel held at a fixed mutation
rate of zero. The QTL and mutating-marker panels are each split into "active" sites, which carry
nearly all of the founder population's additive variance, and "reservoir" sites, drawn from the
same oversampled pool at near-boundary founder frequency but still assigned effects at
initialization, so that when recurrent mutation later activates one, it already has an effect
size rather than needing one assigned dynamically. The validation-marker panel exists specifically
so the neutral heterozygosity-decay theory check is not confounded by mutation pushing allele
frequencies toward a mutation-drift equilibrium rather than following the pure-drift curve the
theory assumes.

**Recurrent mutation.** A fixed founder haplotype set exhausts standing variation over long
horizons, and selection response flatlines as a modeling artifact rather than a biological one
[@robertson1960limits; @hill1986recurrent]. Mutation is applied each generation with AlphaSimR's
`mutate()` function, which flips individual haplotype copies rather than forcing an individual
homozygous, confirmed empirically before relying on it (`dev/verify_mutate_heterozygosity.R`):
after a `mutate()` call, converted homozygous sites showed real heterozygous outcomes matching
what a single-copy mutation event should produce, not the artificial homozygosity that
AlphaSimR's `editGenome()` would have introduced instead.

The mutation rate is calibrated from a target mutational heritability rather than a molecular
per-bp rate, since each site here is an effective locus abstracting many base pairs and a
genomic per-bp rate of roughly 1e-8 under-counts mutational input by orders of magnitude
[@houle1996mutational; @lynch1998qtl]. With founder additive variance $V_A = 1$ and heritability
$h^2 = 0.3$: $V_P = V_A / h^2 = 3.333$, $V_E = V_P - V_A = 2.333$. Targeting a mutational
heritability of $h^2_m = V_M / V_E = 10^{-3}$ gives $V_M = 2.333 \times 10^{-3}$ per generation.
Using $V_M = 2 \mu n E[a^2]$ with $n = 50{,}000$ mutable QTL sites and $E[a^2] \approx V_A /
n_{\text{active}} = 10^{-4}$ (the approximation that founder additive variance is carried almost
entirely by the active sites) gives $\mu \approx 2.3 \times 10^{-4}$ per site per generation. This
is a calibrated starting point, not a derived-exact value: the actual pass/fail signal is whether
additive variance empirically stabilizes at a non-trivial level over a few hundred generations
rather than decaying to zero, which is what the mutation-selection-balance check
(`dev/verify_va_stability.R`) tests directly, and $\mu$ is meant to be retuned from that result if
needed.

**Generation loop.** Each generation applies mutation, assigns phenotypes at the configured
heritability, selects parents by phenotype at the configured intensity, and produces the next
generation by random crossing. The loop is serial within a replicate, since generation $t+1$
depends on $t$, but replicates and scenarios are independent and safe to parallelize.

## Reproducibility design and its limits

**Single-threaded execution.** `SimParam$nThreads` and `runMacs()`'s own separate `nThreads`
argument both default to all available cores via OpenMP. AlphaSimR's own documentation examples
all carry a hidden `SP$nThreads = 1L`, an implicit acknowledgment from the maintainers that
multi-threaded execution is not reproducible from a seed. Both are forced to 1 in this project.
This was not a precaution taken on faith: multi-threading was confirmed empirically to break
reproducibility before it was disabled, and disabling it was confirmed to fix everything
downstream of founder generation (see below).

**Independent RNG streams.** Each replicate draws from its own L'Ecuyer-CMRG stream
[@lecuyer2002streams; @lecuyer1999cmrg], generated via `parallel::nextRNGStream()` rather than
distinct integer seeds, so that replicates cannot be correlated with each other regardless of how
many are run. `validation/test_replicate_independence.R` confirms `make_replicate_streams()`
hands out distinct seed vectors per replicate and that two replicates in the same run produce
different founders and different trajectories, while a replicate resumed from its own checkpoint
reproduces itself exactly.

**Checkpoint and resume.** A checkpoint saves the population object, its `SimParam`, the chip
index map, every output row accumulated so far, and the full RNG state, `.Random.seed`, not a
seed integer, together (`R/checkpoint.R`). Restoring that state in a fresh process needs one
non-obvious step: `RNGkind("L'Ecuyer-CMRG")` has to be set before `.Random.seed` is assigned, or
the restored state silently fails to take effect in a session that never switched to that
generator, a bug that was caught, not assumed absent, by a test that runs the resumed half in a
genuinely separate Rscript process
(`validation/test_resume_equals_uninterrupted.R`). With that fix in place, a resumed run
reproduces an uninterrupted one bit-for-bit, verified on the actual genotype matrices and the
full recorded output series, not just summary statistics.

**Founder generation is not reproducible from a seed on the AlphaSimR release we pin.**
`AlphaSimR:::MaCS()`, the compiled coalescent function `runMacs()` calls, does not produce
identical output across calls given identical literal seeds, confirmed directly and bypassing R's
own RNG entirely (`dev/verify_macs_direct.R`). This holds even single-threaded, so it is not the
same issue as the threading problem above; it is a separate defect in AlphaSimR's compiled code.

It is a documented upstream bug with an upstream fix that has not yet reached CRAN. Reported as
[issue #228](https://github.com/gaynorr/AlphaSimR/issues/228) and diagnosed in [issue
#266](https://github.com/gaynorr/AlphaSimR/issues/266) by AlphaSimR co-author Gregor Gorjanc: the
compiled code draws from Armadillo's RNG, whose state is separate from R's, so `set.seed()` never
reaches it. [PR #265](https://github.com/gaynorr/AlphaSimR/pull/265) fixes this by routing seeds
through R's RNG (`runMacs()` now draws them with `sample.int()` and passes them numerically, where
1.5.3 built them independently), and was merged to the `devel` branch on 2026-04-03.

We verified the fix rather than assuming it from the thread
(`dev/verify_macs_direct_current_version.R`). On `devel` (2.1.0.9004, containing PR #265) both the
direct `MaCS()` call and issue #228's own `set.seed()` + `runMacs()` repro are reproducible,
single- and multi-threaded, with a negative control confirming distinct seeds still diverge. Both
released versions still fail the same tests: the pinned 1.5.3, and current CRAN 2.1.0, which is
built from `master` and predates the `devel` merge. Upgrading the pin to current CRAN would
therefore not resolve it, and the founder-persistence workaround remains necessary until the fix
appears in a CRAN release.

The practical consequence, for as long as we pin an affected release: a scenario config and a
seed reproduce the statistical character of a run, everything after founders exist, not the
specific run. `run.R` persists each replicate's founder population once, immediately after it is
built, to
`data/founders/<scenario>_rep<n>.rds`, specifically so a given run's founders can be recovered
exactly when needed. State this plainly, since it is easy to miss if only the validation suite's
PASS output is read: **the pattern is robust across runs; the exact numbers for a specific run
come from the persisted founder object, not the seed.**

## Validation results

Four checks, each run separately from the main pipeline and each reporting the actual numbers it
produced, not only pass or fail.

**Neutral heterozygosity decay** [@falconer1996qtg; @lynch1998qtl]. Expected heterozygosity on the
$\mu = 0$ validation-marker panel should decay at $(1 - 1/(2N_e))$ per generation under pure
drift. Fitting a log-linear decay to five replicates under no selection gives an implied $N_e$ of
104.0 (range 89.4 to 125.0) against a nominal $N_e$ of 100, close agreement. Under low-intensity
selection the same fit gives an implied $N_e$ of 14.1 (range 9.9 to 17.6): only a fifth of the
population contributes parents each generation, and that concentration compounds drift far beyond
what census size alone would suggest. The gap between the two conditions is not an error; it is
the point of running the selected condition at all.

**Breeder's equation**, $R = h^2 S$ [@falconer1996qtg]. Predicting response from a fixed nominal
heritability tracks the actual response well early and increasingly overpredicts as additive
variance depletes. The cumulative gap between predicted and actual response was 2.5% by generation
5, 20.6% by generation 10, 94.1% by generation 25, and 476.5% by generation 100. Realized
heritability, $V_A / V_P$ measured directly rather than assumed, fell from a mean of 0.237 in the
first ten generations to 0.021 in the last ten, against the nominal 0.300 the prediction used
throughout: the equation does not fail because the underlying biology changes, it fails because a
single early heritability estimate stops describing the population.

**$F_{ROH}$ against pedigree $F$** [@curik2014roh; @keller2011roh; @kardos2015genomic]. Under 40
generations of deliberately intense inbreeding (a founder population of 40, 8 parents per
generation), both measures rise together, from a mean of 0.27 ($F_{ROH}$) and 0.22 (pedigree) at
generation 5 to 0.94 and 0.90 by generation 40. At the final generation the two correlate at
$r \approx 0.34$ to $0.37$ across runs but are not identical: $F_{ROH}$ sits consistently above
pedigree $F$ by roughly 0.04 to 0.07, consistent with realized recombination and Mendelian
sampling variance pushing genomic inbreeding past the population-average expectation that pedigree
$F$ represents.

**Inbreeding depression** [@doekes2021meta]. A purely additive trait cannot show mean inbreeding
depression by construction, since depression requires directional dominance: the higher-value
allele tending to be dominant, so rising homozygosity unmasks more of the lower-value allele's
effect. This was tested, not assumed: pooling ten replicates (12,500 observations) of a purely
additive population and regressing phenotype on pedigree $F$, controlling for a replicate by
generation cohort effect to remove genetic drift's own correlation with generation number, gives a
non-significant slope (0.472, SE 0.372, $p = 0.205$). The same design with a directional-dominance
architecture (`SP\$addTraitAD`, mean dominance degree 0.6) gives a highly significant negative
slope (-2.369, SE 0.364, $p < 0.0001$), about -1.2 phenotypic standard deviations per unit rise in
$F$. The comparison to the Doekes et al. benchmark of roughly 0.13% decline per 1% rise in $F$ is
qualitative, not literal: the simulated trait has an arbitrary mean-zero scale, not a natural
biological unit like birth weight, so a percent-of-mean comparison is not well defined here, and
was confirmed numerically unstable before the standard-deviation framing was adopted instead. The
directional-dominance architecture used for this test is not the one the main pipeline's scenarios
run under; extending it there is a separate decision, not yet made.

One methodological point worth stating directly, since it shaped the final test design: pooling
generations and regressing phenotype on $F$ without controlling for generation is confounded,
because $F$ and generation number are nearly perfectly correlated in an inbreeding design, and
genetic drift's own random walk in mean breeding value over generations leaks into the estimated
$F$ effect. This was caught, not avoided in advance, when the purely additive condition first
showed a large, significant slope, which is not mechanistically possible for a purely additive
trait and was the signal that something in the regression, not the biology, was wrong.
