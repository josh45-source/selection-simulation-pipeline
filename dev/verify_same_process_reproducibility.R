source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")
source("R/checkpoint.R")
source("R/sim_loop.R")

config <- list(
  scenario_id = "repro_test",
  replicate = list(n_replicates = 1, seed = 909090),
  population = list(n_founders = 20, species = "generic_diploid"),
  genome = list(n_chr = 2, segsites_per_chr = 210),
  qtl = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, heritability = 0.3,
             mean = 0, var_add = 1, mutation_rate = 2.3e-4),
  neutral_markers = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, mutation_rate = 2.3e-4),
  validation_markers = list(n_per_chr = 10, mutation_rate = 0),
  generations = list(n_generations = 10, checkpoint_every = 10),
  sampling = list(locus_cadence = 5, segregating_only = TRUE)
)

stream <- make_replicate_streams(config$replicate$seed, 1)[[1]]

run1 <- run_replicate(config, selection_fraction = 0.2, replicate_id = 1, rng_stream = stream)
run2 <- run_replicate(config, selection_fraction = 0.2, replicate_id = 1, rng_stream = stream)

geno1 <- pullQtlGeno(run1$pop, simParam = run1$SP)
geno2 <- pullQtlGeno(run2$pop, simParam = run2$SP)

cat("identical(geno1, geno2):", identical(geno1, geno2), "\n")
cat("Mismatched cells:", sum(geno1 != geno2), "of", length(geno1), "\n")
cat("identical(summary1, summary2):", identical(run1$summary, run2$summary), "\n")
