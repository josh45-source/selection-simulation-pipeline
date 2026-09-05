# Per-generation metrics: additive genetic variance and per-locus allele frequency/fixation.
#
# F_ROH lives in R/froh.R (detectRUNS-based) and the inbreeding-depression regression in
# validation/test_inbreeding_depression.R, rather than here. Their citations are wired into
# analysis/explore.qmd against references.bib: curik2014roh, keller2011roh, and kardos2015genomic
# for ROH-based genomic inbreeding, doekes2021meta as the inbreeding-depression benchmark.

library(AlphaSimR)

# Additive genetic variance. Uses varA(), not varG(): for the purely additive traits this
# pipeline runs (SP$addTraitA) the two are numerically identical -- verified equal to within
# floating-point error across 100 generations of selection and drift -- so this returns exactly
# what varG() did and no reported number changes. They diverge the moment a trait carries
# dominance or epistasis, where varG() is total genetic variance and varA() is the additive
# component this column claims to be. validation/test_inbreeding_depression.R already builds an
# addTraitAD population (it does not call this function), and PROPOSED_EXTENSIONS.md item 1
# proposes dominance in the main pipeline, so the wrong one here is a latent bug rather than a
# hypothetical.
#
# varA() needs the SimParam; varG() did not. That is the only reason this takes a second argument.
compute_va <- function(pop, simParam) {
  # varA() returns a Trait x Trait var-covar matrix (with dimnames like "Trait1") even for a
  # single trait -- strip that down to a plain unnamed scalar so it doesn't leak AlphaSimR's
  # internal trait naming into the Parquet schema (data.frame(va = <named matrix>) would otherwise
  # name the column after the matrix's own dimname, not the "va" argument name).
  as.numeric(varA(pop, simParam = simParam))[1]
}

# geno: individuals x sites matrix, values in {0,1,2} (e.g. from pullQtlGeno/pullSnpGeno)
allele_freq <- function(geno) {
  colMeans(geno) / 2
}

is_fixed <- function(freq, tol = 0) {
  freq <= tol | freq >= (1 - tol)
}

# Expected heterozygosity / gene diversity (2pq averaged across sites), the quantity the
# H_t = H_0 * (1 - 1/(2*Ne))^t neutral decay theory predicts -- not observed (per-genotype)
# heterozygosity, which conflates with within-population inbreeding.
# ref: falconer1996qtg, lynch1998qtl
heterozygosity <- function(geno) {
  freq <- allele_freq(geno)
  mean(2 * freq * (1 - freq))
}
