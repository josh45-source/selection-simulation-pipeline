# Proposed extensions

A menu, not a plan. Nothing here is implemented, and nothing here should be implemented on the
strength of this document alone — it exists so a domain expert can react to it and say which
directions are scientifically worth the cost.

Each entry states what it would add, what *new* validation it would require, and the tradeoff.
The validation column is the important one: this project's standard is that a scientific claim
ships with an empirical check that reports real numbers, so an extension's true cost is the
validation it forces, not the API calls it takes. Feasibility is assessed against the AlphaSimR
version this project pins (1.5.3); where an extension needs something AlphaSimR does not provide,
that is stated rather than glossed.

Two constraints apply across everything below:

- **Founder non-determinism.** Until the upstream MaCS fix reaches CRAN (see README Limitations),
  any extension that needs a *specific* run reproduced must go through the persisted founder
  objects. This makes cross-scenario comparisons more awkward than they look.
- **Single-threaded execution.** `SP$nThreads = 1L` is required for reproducibility. Extensions
  that multiply the compute cost per generation cannot be bought back with threads.

---

## 1. Dominance in the main pipeline

**What it would add.** The pipeline's scenarios use a purely additive trait (`SP$addTraitA`), which
cannot show mean inbreeding depression by construction. `validation/test_inbreeding_depression.R`
already builds its own directional-dominance population with `SP$addTraitAD(meanDD=)` to
demonstrate the mechanism. Promoting that architecture to the main scenarios would let inbreeding
depression appear in the pipeline's own output rather than only in a test, making the
depression-versus-response tradeoff a first-class result.

**New validation required.** The existing depression test would need re-anchoring against the main
pipeline's trait rather than its own. `compute_va()` currently reports `varG()`, which is exactly
equal to `varA()` for a purely additive trait but *not* under dominance — a dominance pipeline
needs additive and dominance variance separated (`varA()`, `varD()`, and the `genicVar*`
counterparts) and a check that the split behaves as theory predicts as inbreeding rises.
Heritability bookkeeping also changes: `setPheno(h2=)` is narrow-sense and computed from founder
variances, so `H2` versus `h2` becomes a real choice that needs stating.

**Tradeoff.** This is the highest-value extension on the list and also the one that most disturbs
existing validated numbers: every result currently reported would shift, and the neutral-decay and
breeder's-equation tests would need re-running and re-interpreting. It is a re-baselining, not an
addition.

## 2. Genomic selection

**What it would add.** Selection on estimated breeding values instead of phenotypes. AlphaSimR
supports this directly: `RRBLUP()` / `fastRRBLUP()` to train, `setEBV()` to attach EBVs, and
`selectInd(use="ebv")` to select on them. This would let the project compare phenotypic and
genomic selection under the same drift, mutation, and inbreeding machinery it has already
validated — a genuinely interesting contrast, since genomic selection is known to accelerate both
response *and* inbreeding.

**New validation required.** Prediction accuracy (correlation of EBV with true breeding value)
must be measured and reported per generation, not assumed; accuracy decay as training data ages
is the failure mode that makes naive genomic-selection simulations wrong. The existing
breeder's-equation check would need a genomic analogue, and the F_ROH machinery becomes
load-bearing rather than confirmatory, since the headline claim would be about inbreeding cost.
The mutating neutral-marker chip becomes a real training panel, which means the SNP-chip
allocation issue noted in the audit would have to be settled first.

**Tradeoff.** Substantial compute per generation (a model fit per cycle, single-threaded), and it
changes the project's identity from "long-term selection dynamics" to "breeding-scheme
comparison." Scientifically strong, but it is a second project sharing a codebase.

## 3. Multi-trait selection and index weights

**What it would add.** Correlated traits via `addTraitA(corA=)`, selected on an index
(`selIndex()`, `smithHazel()`). This makes antagonistic correlations expressible — the classic
case where selecting hard on one trait drags another down — which single-trait simulations cannot
represent at all.

**New validation required.** A correlated-response check: realized correlated response against the
multivariate breeder's equation, with the realized genetic correlation tracked over time, since
selection erodes it. The neutral-decay and F_ROH tests carry over unchanged.

**Tradeoff.** Moderate cost, well-supported by the API, and it composes cleanly with the existing
tests. The risk is interpretive rather than technical: index weights are an arbitrary choice, and
results become a function of that choice unless the weights are justified externally.

## 4. IBD-based inbreeding as a third estimator

**What it would add.** This project already computes pedigree F and F_ROH and reports that they
correlate but are not identical. AlphaSimR can track true identity-by-descent directly
(`SP$setTrackRec(TRUE)`, then `pullIbdHaplo()`), which gives the *true* IBD inbreeding the other
two are estimating. That converts an observed discrepancy into a measured bias: which estimator is
closer to truth, and how does the gap behave as inbreeding accumulates.

**New validation required.** Comparatively little — this *is* validation. It would need a check
that pedigree F and F_ROH both track IBD F in the expected direction, and an honest statement of
where each fails. The existing F_ROH test would gain a ground-truth reference rather than being
replaced.

**Tradeoff.** The cheapest scientifically meaningful item on this list and the best fit with what
the project already does. Cost is memory and runtime: recombination tracking is not free over 100
generations, and would need a slice-scale feasibility check before committing.

## 5. Optimal contribution selection

**What it would add.** Selection that maximises genetic gain subject to an explicit constraint on
the rate of inbreeding — the standard tool for the exact tradeoff this project measures. It is the
natural endpoint of the inbreeding-versus-response story.

**New validation required.** That realised ΔF actually tracks the constraint the optimiser was
given; that is the whole claim, and it is easy to get silently wrong. It would also need a
constrained-versus-unconstrained comparison to show the method costs something in response.

**Tradeoff.** The honest blocker: **AlphaSimR has no native OCS.** `selectCross()` is a
convenience wrapper around `selectInd()` + `randCross()`, not an optimiser. Real OCS needs
`optiSel`, which the README already records as unresolvable from the pinned package repository for
R 4.3.3, or a hand-rolled optimiser over the numerator relationship matrix — which is a
substantial piece of numerical work that would itself need validating against a published
benchmark before any result using it could be trusted. Highest scientific value per result,
highest implementation risk on this list.

## 6. Other species histories and real genomes

**What it would add.** `runMacs()` ships CATTLE, WHEAT, and MAIZE demographic histories alongside
GENERIC, and `runMacs2()` allows a custom effective-population-size history. Beyond that,
`importGenMap()`, `importHaplo()`, and `importInbredGeno()` accept real marker maps and
haplotypes, which would let the pipeline run on an actual species' genome architecture rather than
a generic one.

**New validation required.** Every existing test is calibrated against GENERIC's parameters and
would need re-running per history — in particular the F_ROH configuration, whose `maxGap` is
scaled to this genome's marker spacing and whose 1 Morgan : 100 Mb physical conversion is a
GENERIC-specific assumption that is simply wrong for a real map. Marker density changes ROH
detection behaviour more than it changes anything else here.

**Tradeoff.** Low implementation cost, high credibility gain for a domain audience — a cattle
history makes the inbreeding-depression comparison to the Doekes benchmark far less abstract. The
manual explicitly cautions that the WHEAT and MAIZE histories do not faithfully represent their
species, so those two should be treated as legacy settings rather than realism.

## 7. GxE and epistasis

**What it would add.** `addTraitAG()` / `addTraitADG()` add genotype-by-environment interaction
with an environmental covariate sampled per generation (`setPheno(p=)`); the `TraitADE` family adds
epistatic variance (`varAA()`, `genicVarAA()`). Both make the genetic architecture more realistic
and both change how variance depletes under long-term selection.

**New validation required.** For GxE, that the realised environmental variance component matches
what was specified, and a response comparison across environments. For epistasis, that additive and
epistatic variance are separated correctly and that the well-known conversion of epistatic to
additive variance under drift actually appears — which is a genuinely demanding check to design.

**Tradeoff.** Epistasis in particular buys realism at a steep interpretive cost: results become
strongly dependent on architecture parameters that are hard to justify from data, and the project's
current strength is that its claims are simple enough to check. Recommended only if a specific
question requires it.

## 8. Breeding-scheme structure

**What it would add.** Doubled haploids (`makeDH()`), selfing (`self()`), hybrid crossing with
general and specific combining ability (`hybridCross()`, `calcGCA()`, `setPhenoGCA()`), family
structure (`selectFam()`, `selectWithinFam()`), and the usefulness criterion (`usefulness()`).
Together these would let the pipeline represent an actual breeding programme rather than
generation-on-generation mass selection.

**New validation required.** Scheme-specific expectations for inbreeding accumulation, which
differ sharply from the random-mating baseline the current tests assume. The pedigree F
implementation would need extending: it is a tabular-method calculator written specifically
because AlphaSimR is monoecious, and selfing and DH lines stress exactly that code path.

**Tradeoff.** Each piece is individually cheap, but they only mean something in combination, and
the combination is a different project. Worth listing so it can be explicitly declined.

---

## If only one thing

**Item 4 (IBD-based inbreeding)** is the recommendation: it is the cheapest, it strengthens a claim
the project already makes rather than opening a new front, and it turns the existing "F_ROH and
pedigree F correlate but differ" result into a statement about which one is right. **Item 1
(dominance)** is the highest-value if a re-baselining is acceptable. **Item 5 (OCS)** is the most
scientifically attractive and the one most likely to consume the project.
