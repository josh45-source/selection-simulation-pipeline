# Confirms replicates don't accidentally share founders or correlated RNG streams -- the flip
# side of the reproducibility bugs fixed in the checkpoint/resume work (where the goal was making
# a SINGLE replicate reproducible; here the goal is confirming DIFFERENT replicates are NOT
# reproductions of each other).
#
# Caveat that matters for reading this test: AlphaSimR:::MaCS() (runMacs's founder generator) is
# not deterministic even given identical seeds (see dev/verify_macs_direct.R) -- so "founder
# genotypes differ between replicates" is expected regardless of whether streams are properly
# independent, and isn't by itself proof of non-overlapping L'Ecuyer streams. What IS a clean,
# deterministic check is that make_replicate_streams() actually hands out distinct seed vectors
# per replicate, and that a full run doesn't pathologically share state across replicates (e.g. a
# bug in run.R's loop reusing one replicate's post-run RNG state for the next). Both are checked
# below; the founder/trajectory comparisons are reported as corroborating evidence, not proof.

library(testthat)

source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")
source("R/checkpoint.R")
source("R/sim_loop.R")

test_that("make_replicate_streams gives distinct, non-identical seeds per replicate", {
  streams <- make_replicate_streams(909090, 5)
  for (i in 2:5) {
    expect_false(identical(streams[[1]], streams[[i]]))
  }
  # pairwise, not just each-vs-first
  for (i in 1:4) for (j in (i + 1):5) {
    expect_false(identical(streams[[i]], streams[[j]]))
  }
})

test_that("two replicates in the same run produce different founders and different trajectories", {
  config <- list(
    scenario_id = "indep_test",
    replicate = list(n_replicates = 2, seed = 909090),
    population = list(n_founders = 20, species = "generic_diploid"),
    genome = list(n_chr = 2, segsites_per_chr = 210),
    qtl = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, heritability = 0.3,
               mean = 0, var_add = 1, mutation_rate = 2.3e-4),
    neutral_markers = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, mutation_rate = 2.3e-4),
    validation_markers = list(n_per_chr = 10, mutation_rate = 0),
    generations = list(n_generations = 15, checkpoint_every = 15),
    sampling = list(locus_cadence = 5, segregating_only = TRUE)
  )

  streams <- make_replicate_streams(config$replicate$seed, 2)

  rep1 <- run_replicate(config, selection_fraction = 0.2, replicate_id = 1, rng_stream = streams[[1]])
  rep2 <- run_replicate(config, selection_fraction = 0.2, replicate_id = 2, rng_stream = streams[[2]])

  geno1 <- pullQtlGeno(rep1$pop, simParam = rep1$SP)
  geno2 <- pullQtlGeno(rep2$pop, simParam = rep2$SP)
  cat("Final-generation QTL genotype mismatches between rep1/rep2:",
      sum(geno1 != geno2), "of", length(geno1), "\n")
  expect_false(identical(geno1, geno2))

  cat("rep1 va trajectory:", round(rep1$summary$va, 5), "\n")
  cat("rep2 va trajectory:", round(rep2$summary$va, 5), "\n")
  expect_false(identical(rep1$summary$va, rep2$summary$va))

  # Corroborating check, MaCS-noise-free by construction: build founders ONCE, checkpoint them,
  # then run twice from that SAME checkpoint. This does NOT call build_founder_pop() twice with
  # the same stream (that would fail due to MaCS's own non-determinism regardless of replicate
  # independence -- see the file-level caveat). It isolates whether, given the identical starting
  # population, the harness reproduces the same trajectory -- confirming the rep1-vs-rep2
  # difference above is attributable to genuinely different input, not harness-level noise.
  tmp_dir <- tempfile("indep_test_")
  dir.create(tmp_dir)
  founders_path <- file.path(tmp_dir, "founders.rds")
  set_replicate_stream(streams[[1]])
  founders1 <- build_founder_pop(config)
  save_checkpoint(founders_path, founders1$pop, founders1$SP, founders1$chips, gen = 0L,
                   summary_rows = vector("list", config$generations$n_generations), locus_rows = list())

  rerun_a <- run_replicate(config, selection_fraction = 0.2, replicate_id = 1,
                            rng_stream = NULL, resume_path = founders_path)
  rerun_b <- run_replicate(config, selection_fraction = 0.2, replicate_id = 1,
                            rng_stream = NULL, resume_path = founders_path)
  geno_a <- pullQtlGeno(rerun_a$pop, simParam = rerun_a$SP)
  geno_b <- pullQtlGeno(rerun_b$pop, simParam = rerun_b$SP)
  cat("Same-founders-checkpoint rerun mismatches:", sum(geno_a != geno_b), "of", length(geno_a), "\n")
  expect_identical(geno_a, geno_b)
  expect_identical(rerun_a$summary, rerun_b$summary)
  unlink(tmp_dir, recursive = TRUE)
})
