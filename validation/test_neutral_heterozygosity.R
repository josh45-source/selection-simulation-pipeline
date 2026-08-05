# Validates neutral heterozygosity decay against theory: under pure drift, expected heterozygosity
# (gene diversity, 2pq) decays as H_t = H_0 * (1 - 1/(2*Ne))^t per generation.
# ref: falconer1996qtg, lynch1998qtl
#
# Uses ONLY the validation_markers panel (mu = 0, protected from mutate() by the snapshot/restore
# in R/mutation.R -- see its file header) -- NOT the mutating neutral_markers panel, which would be
# pulled toward mutation-drift equilibrium rather than following the pure-drift decay curve this
# theory assumes.
#
# Runs two conditions across multiple replicates (never a single trajectory -- see design report):
#   - no_selection: selection_fraction = 1.0 (all individuals become parents), the closest this
#     harness gets to the idealized Wright-Fisher assumptions the theory is derived under. Realized
#     Ne here should track nominal Ne (n_founders) reasonably closely.
#   - low_selection: selection_fraction = 0.2, matching the "low" scenario elsewhere in this repo.
#     Only ~20% of individuals contribute each generation, so realized Ne should come in
#     noticeably BELOW nominal Ne -- that gap is the point of running this condition, not an error.

library(testthat)

source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")

run_condition <- function(config, selection_fraction, replicate_id, rng_stream) {
  set_replicate_stream(rng_stream)
  founders <- build_founder_pop(config)
  pop <- founders$pop
  SP <- founders$SP
  validation_chip <- founders$chips$validation

  n_gen <- config$generations$n_generations
  n_founders <- config$population$n_founders
  h2 <- config$qtl$heritability
  mu <- config$qtl$mutation_rate

  h_trace <- numeric(n_gen + 1L)
  h_trace[1] <- heterozygosity(pullSnpGeno(pop, snpChip = validation_chip, simParam = SP))

  for (gen in seq_len(n_gen)) {
    pop <- apply_recurrent_mutation(pop, SP, mu, validation_chip)
    pop <- setPheno(pop, h2 = h2, simParam = SP)
    n_parents <- max(1, round(n_founders * selection_fraction))
    parents <- selectInd(pop, nInd = n_parents, use = "pheno", simParam = SP)
    pop <- randCross(parents, nCrosses = n_founders, simParam = SP)
    h_trace[gen + 1L] <- heterozygosity(pullSnpGeno(pop, snpChip = validation_chip, simParam = SP))
  }
  h_trace
}

implied_ne <- function(h_trace) {
  t <- seq_along(h_trace) - 1L
  keep <- h_trace > 0 # log-linear fit undefined once fully fixed
  fit <- lm(log(h_trace[keep]) ~ t[keep])
  slope <- unname(coef(fit)[2])
  decay_factor <- exp(slope) # should approximate (1 - 1/(2*Ne))
  ne <- 1 / (2 * (1 - decay_factor))
  list(slope = slope, decay_factor = decay_factor, ne = ne)
}

test_that("neutral heterozygosity decay matches theory, and the selection-induced Ne gap is visible", {
  n_founders <- 100
  n_replicates <- 5
  n_gen <- 100

  config_base <- list(
    scenario_id = "hetdecay_test",
    replicate = list(n_replicates = n_replicates, seed = 424242),
    population = list(n_founders = n_founders, species = "generic_diploid"),
    genome = list(n_chr = 3, segsites_per_chr = 200),
    qtl = list(n_active_per_chr = 10, n_reservoir_per_chr = 40, heritability = 0.3,
               mean = 0, var_add = 1, mutation_rate = 2.3e-4),
    neutral_markers = list(n_active_per_chr = 10, n_reservoir_per_chr = 40, mutation_rate = 2.3e-4),
    validation_markers = list(n_per_chr = 100, mutation_rate = 0),
    generations = list(n_generations = n_gen, checkpoint_every = n_gen),
    sampling = list(locus_cadence = 10, segregating_only = TRUE)
  )

  streams <- make_replicate_streams(config_base$replicate$seed, n_replicates)

  results <- list()
  for (cond in c("no_selection", "low_selection")) {
    frac <- if (cond == "no_selection") 1.0 else 0.2
    ne_estimates <- numeric(n_replicates)
    h_start <- numeric(n_replicates)
    h_end <- numeric(n_replicates)
    for (r in seq_len(n_replicates)) {
      h_trace <- run_condition(config_base, frac, r, streams[[r]])
      est <- implied_ne(h_trace)
      ne_estimates[r] <- est$ne
      h_start[r] <- h_trace[1]
      h_end[r] <- h_trace[length(h_trace)]
    }
    results[[cond]] <- list(ne_estimates = ne_estimates, h_start = h_start, h_end = h_end)

    cat(sprintf(
      "\n[%s] H_0: %.4f (mean)  H_%d: %.4f (mean)\n",
      cond, mean(h_start), n_gen, mean(h_end)
    ))
    cat(sprintf(
      "[%s] Implied Ne per replicate: %s\n",
      cond, paste(round(ne_estimates, 1), collapse = ", ")
    ))
    cat(sprintf(
      "[%s] Implied Ne: mean=%.1f, range=[%.1f, %.1f]  (nominal Ne = %d)\n",
      cond, mean(ne_estimates), min(ne_estimates), max(ne_estimates), n_founders
    ))
  }

  gap_pct <- 100 * (mean(results$low_selection$ne_estimates) - mean(results$no_selection$ne_estimates)) /
    mean(results$no_selection$ne_estimates)
  cat(sprintf(
    "\nRealized Ne gap, low_selection vs no_selection: %.1f%% (nominal Ne = %d in both)\n",
    gap_pct, n_founders
  ))

  # no_selection should track nominal Ne within a generous band -- this is the theory-matches
  # sanity check. Tolerance is loose (50%) because a single log-linear fit over 100 generations
  # from real stochastic trajectories is noisy; the point is order-of-magnitude agreement, not
  # precision.
  expect_true(mean(results$no_selection$ne_estimates) > 0.5 * n_founders)
  expect_true(mean(results$no_selection$ne_estimates) < 2 * n_founders)

  # low_selection's realized Ne should be MEANINGFULLY below nominal -- that's the whole point of
  # running this condition. Not a tight bound, just confirming the gap goes the expected direction
  # and is not negligible.
  expect_true(mean(results$low_selection$ne_estimates) < 0.9 * mean(results$no_selection$ne_estimates))
})
