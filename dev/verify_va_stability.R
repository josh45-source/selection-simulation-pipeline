library(AlphaSimR)

n_founders   <- 100
n_chr        <- 10
n_qtl_total  <- 5000   # per chr: 1000 active + 4000 reservoir, per scenario.example.yaml
n_gen        <- 300
sel_fraction <- 0.20    # "low" intensity scenario
h2           <- 0.3
mu_qtl       <- 2.3e-4  # calibrated in scenario.example.yaml from target h2_m = 1e-3

run_scenario <- function(mutRate, seed, report_freq_spectrum = FALSE) {
  set.seed(seed)
  founderPop <- runMacs(nInd = n_founders, nChr = n_chr, segSites = n_qtl_total, species = "GENERIC")
  SP <- SimParam$new(founderPop)
  SP$addTraitA(n_qtl_total, mean = 0, var = 1)
  SP$setTrackPed(FALSE)

  pop <- newPop(founderPop, simParam = SP)

  if (report_freq_spectrum) {
    geno0 <- pullQtlGeno(pop, simParam = SP)
    p0 <- colMeans(geno0) / 2
    cat("Founder QTL allele-frequency spectrum (runMacs, GENERIC):\n")
    print(summary(p0))
    cat("Sites with freq < 0.02 or > 0.98 (near-monomorphic, natural reservoir):",
        sprintf("%.1f%%\n", 100 * mean(p0 < 0.02 | p0 > 0.98)))
  }

  va_trace <- numeric(n_gen)
  for (gen in seq_len(n_gen)) {
    if (mutRate > 0) {
      pop <- mutate(pop, mutRate = mutRate, simParam = SP)
    }
    pop <- setPheno(pop, h2 = h2, simParam = SP)
    parents <- selectInd(pop, nInd = round(n_founders * sel_fraction), use = "pheno", simParam = SP)
    pop <- randCross(parents, nCrosses = n_founders, simParam = SP)
    va_trace[gen] <- varG(pop)
  }
  va_trace
}

cat("Running WITHOUT mutation (mu = 0) ...\n")
va_no_mut <- run_scenario(0, seed = 1)

cat("Running WITH mutation (mu =", mu_qtl, ") ...\n")
va_mut <- run_scenario(mu_qtl, seed = 1, report_freq_spectrum = TRUE)

report_gens <- seq(20, n_gen, by = 20)
cat("\ngen\tVa_no_mutation\tVa_with_mutation\n")
for (g in report_gens) {
  cat(g, "\t", round(va_no_mut[g], 5), "\t\t", round(va_mut[g], 5), "\n", sep = "")
}

va0 <- mean(va_no_mut[1:10])
cat("\nFounder-era Va (mean of gen 1-10):", round(va0, 5), "\n")
cat("No-mutation Va at gen", n_gen, ":", round(va_no_mut[n_gen], 5),
    sprintf(" (%.1f%% of founder Va)\n", 100 * va_no_mut[n_gen] / va0))
cat("With-mutation Va at gen", n_gen, ":", round(va_mut[n_gen], 5),
    sprintf(" (%.1f%% of founder Va)\n", 100 * va_mut[n_gen] / va0))

plateaued <- mean(va_mut[(n_gen - 49):n_gen])
decayed   <- mean(va_no_mut[(n_gen - 49):n_gen])
cat("\nMean Va, last 50 gens -- no mutation:", round(decayed, 5),
    " | with mutation:", round(plateaued, 5), "\n")
if (plateaued > 5 * decayed && plateaued > 0.05 * va0) {
  cat("CONCLUSION: mutation model plateaus Va at a non-trivial level vs. decay -- PASS.\n")
} else {
  cat("CONCLUSION: mutation does not meaningfully rescue Va from decay at this rate/site count -- FAIL, needs recalibration.\n")
}
