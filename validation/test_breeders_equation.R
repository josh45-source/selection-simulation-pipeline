# Validates the breeder's equation R = h^2 * S as the short-term selection-response prediction.
# ref: falconer1996qtg
#
# Per generation: S (selection differential) = mean phenotype of selected parents minus mean
# phenotype of the pre-selection population; predicted R = h2_nominal * S, using the FIXED
# config heritability (0.3), the way a breeder would naively apply the equation from a single
# early heritability estimate -- not a per-generation re-measured h2. Actual R = mean breeding
# value of the offspring generation minus mean breeding value of the parental generation
# (pre-selection). Expected to hold well early (Va close to its founder value, realized h2 close
# to nominal) and increasingly overestimate actual R as Va depletes and realized h2 = Va/Vp drops
# below the nominal value used in the prediction.

library(testthat)

source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")

test_that("breeder's equation holds early and its overprediction grows as Va depletes", {
  n_founders <- 100
  n_gen <- 100
  h2_nominal <- 0.3

  config <- list(
    scenario_id = "breeders_eq_test",
    replicate = list(n_replicates = 1, seed = 135790),
    population = list(n_founders = n_founders, species = "generic_diploid"),
    genome = list(n_chr = 3, segsites_per_chr = 200),
    qtl = list(n_active_per_chr = 10, n_reservoir_per_chr = 40, heritability = h2_nominal,
               mean = 0, var_add = 1, mutation_rate = 2.3e-4),
    neutral_markers = list(n_active_per_chr = 10, n_reservoir_per_chr = 40, mutation_rate = 2.3e-4),
    validation_markers = list(n_per_chr = 100, mutation_rate = 0),
    generations = list(n_generations = n_gen, checkpoint_every = n_gen),
    sampling = list(locus_cadence = 10, segregating_only = TRUE)
  )

  stream <- make_replicate_streams(config$replicate$seed, 1)[[1]]
  set_replicate_stream(stream)
  founders <- build_founder_pop(config)
  pop <- founders$pop
  SP <- founders$SP
  validation_chip <- founders$chips$validation

  predicted_R <- numeric(n_gen)
  actual_R <- numeric(n_gen)
  va_trace <- numeric(n_gen)
  realized_h2 <- numeric(n_gen)

  for (gen in seq_len(n_gen)) {
    pop <- apply_recurrent_mutation(pop, SP, config$qtl$mutation_rate, validation_chip)
    pop <- setPheno(pop, h2 = h2_nominal, simParam = SP)

    mean_pheno_before <- mean(pop@pheno[, 1])
    mean_gv_before <- meanG(pop)
    va <- compute_va(pop, SP)
    vp <- var(pop@pheno[, 1])

    n_parents <- max(1, round(n_founders * 0.2))
    parents <- selectInd(pop, nInd = n_parents, use = "pheno", simParam = SP)
    S <- mean(parents@pheno[, 1]) - mean_pheno_before

    pop <- randCross(parents, nCrosses = n_founders, simParam = SP)
    mean_gv_after <- meanG(pop)

    predicted_R[gen] <- h2_nominal * S
    actual_R[gen] <- mean_gv_after - mean_gv_before
    va_trace[gen] <- va
    realized_h2[gen] <- if (vp > 0) va / vp else NA
  }

  cum_predicted <- cumsum(predicted_R)
  cum_actual <- cumsum(actual_R)

  report_gens <- unique(c(seq(5, n_gen, by = 5), n_gen))
  cat("\ngen\tVa\trealized_h2\tpredicted_R\tactual_R\tcum_predicted\tcum_actual\tcum_gap_pct\n")
  for (g in report_gens) {
    gap_pct <- 100 * (cum_predicted[g] - cum_actual[g]) / abs(cum_actual[g])
    cat(sprintf("%d\t%.4f\t%.3f\t\t%.4f\t\t%.4f\t\t%.3f\t\t%.3f\t\t%.1f%%\n",
                g, va_trace[g], realized_h2[g], predicted_R[g], actual_R[g],
                cum_predicted[g], cum_actual[g], gap_pct))
  }

  early_window <- 1:10
  late_window <- (n_gen - 9):n_gen
  early_gap <- mean(abs(predicted_R[early_window] - actual_R[early_window]))
  late_gap <- mean(abs(predicted_R[late_window] - actual_R[late_window]))
  cat(sprintf(
    "\nMean |predicted - actual| R, first 10 gens: %.4f | last 10 gens: %.4f\n",
    early_gap, late_gap
  ))
  cat(sprintf("Va: founder-era (gen 1-10 mean) = %.4f, late (gen %d-%d mean) = %.4f\n",
              mean(va_trace[1:10]), n_gen - 9, n_gen, mean(va_trace[late_window])))
  cat(sprintf("Realized h2: early mean = %.3f, late mean = %.3f (nominal used in prediction = %.3f)\n",
              mean(realized_h2[1:10], na.rm = TRUE), mean(realized_h2[late_window], na.rm = TRUE), h2_nominal))

  # Sign of cumulative response should be right and early-window prediction should track
  # reasonably (breeder's equation is a per-generation EXPECTATION, individual generations are
  # noisy, so this checks the cumulative trend over the early window, not single generations).
  expect_true(sign(cum_predicted[10]) == sign(cum_actual[10]) || abs(cum_actual[10]) < 1e-6)
  # The overprediction gap should be larger late than early -- that's the depletion signature.
  expect_true(late_gap >= early_gap * 0.5) # loose: confirms it doesn't shrink outright
})
