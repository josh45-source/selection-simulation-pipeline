# Metrics for the smallest end-to-end slice: additive genetic variance and per-locus allele
# frequency/fixation. F_ROH and inbreeding-depression regression are layered on in a later step
# once the slice (sim -> Parquet -> DuckDB -> plot) is proven end-to-end.
#
# TODO(cite): when F_ROH is implemented here (via detectRUNS), cite curik2014roh, keller2011roh
# for ROH as the genomic inbreeding estimator, and kardos2015genomic, doekes2019recent for using
# genomic/ROH-based F over pedigree-only F. If a GRM-based F is added, cite vanraden2008grm.
# If the inbreeding-depression regression is added, cite doekes2021meta as the calibration
# benchmark (~0.13% trait decline per 1% rise in F).

library(AlphaSimR)

compute_va <- function(pop) {
  # varG() returns a Trait x Trait var-covar matrix (with dimnames like "Trait1") even for a
  # single trait -- strip that down to a plain unnamed scalar so it doesn't leak AlphaSimR's
  # internal trait naming into the Parquet schema (data.frame(va = <named matrix>) would otherwise
  # name the column after the matrix's own dimname, not the "va" argument name).
  as.numeric(varG(pop))[1]
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
