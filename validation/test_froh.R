# Validates F_ROH (genomic inbreeding from runs of homozygosity, via detectRUNS) against pedigree
# F (expected inbreeding from relationship structure, via kinship2). They should correlate but are
# NOT the same quantity -- F_ROH captures realized IBD from actual recombination and Mendelian
# sampling variance, pedigree F is only the population-average expectation for that relationship.
# ref: curik2014roh, keller2011roh (ROH as the genomic inbreeding estimator over pedigree-only F);
# kardos2015genomic (why genomic/ROH-based F is preferred over pedigree F)
#
# Uses the validation_markers panel (mu = 0) so ROH detection isn't confounded by mutation-induced
# homozygosity artifacts (see R/mutation.R) -- the same reason this panel exists for the neutral
# heterozygosity-decay test.

library(testthat)

source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")
source("R/froh.R")
source("R/pedigree.R")

test_that("F_ROH correlates with but is not identical to pedigree F, and both rise under inbreeding", {
  n_founders <- 40
  n_gen <- 40
  n_parents <- 8 # intense selection on a small population, to generate real inbreeding quickly

  config <- list(
    scenario_id = "froh_test",
    replicate = list(n_replicates = 1, seed = 246810),
    population = list(n_founders = n_founders, species = "generic_diploid"),
    genome = list(n_chr = 3, segsites_per_chr = 200),
    qtl = list(n_active_per_chr = 10, n_reservoir_per_chr = 40, heritability = 0.3,
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

  ped_all <- data.frame(id = as.numeric(pop@id), mother = as.numeric(pop@mother),
                         father = as.numeric(pop@father))

  froh_mean_trace <- numeric(n_gen)
  fped_mean_trace <- numeric(n_gen)

  for (gen in seq_len(n_gen)) {
    pop <- apply_recurrent_mutation(pop, SP, config$qtl$mutation_rate, validation_chip)
    pop <- setPheno(pop, h2 = 0.3, simParam = SP)
    parents <- selectInd(pop, nInd = n_parents, use = "pheno", simParam = SP)
    pop <- randCross(parents, nCrosses = n_founders, simParam = SP)

    ped_all <- rbind(ped_all, data.frame(
      id = as.numeric(pop@id), mother = as.numeric(pop@mother), father = as.numeric(pop@father)
    ))

    froh_df <- compute_froh(pop, SP, validation_chip, minSNP = 15, minLengthBps = 50000)
    fped_df <- pedigree_inbreeding(ped_all)
    fped_now <- fped_df[fped_df$id %in% as.numeric(pop@id), ]

    froh_mean_trace[gen] <- mean(froh_df$froh)
    fped_mean_trace[gen] <- mean(fped_now$f_pedigree)
  }

  report_gens <- unique(c(seq(5, n_gen, by = 5), n_gen))
  cat("\ngen\tmean_F_ROH\tmean_F_pedigree\n")
  for (g in report_gens) {
    cat(sprintf("%d\t%.4f\t\t%.4f\n", g, froh_mean_trace[g], fped_mean_trace[g]))
  }

  # Final-generation individual-level comparison
  froh_final <- compute_froh(pop, SP, validation_chip, minSNP = 15, minLengthBps = 50000)
  fped_final_df <- pedigree_inbreeding(ped_all)
  fped_final <- fped_final_df[fped_final_df$id %in% as.numeric(pop@id), ]
  merged <- merge(froh_final, fped_final, by = "id")

  cat("\nFinal generation, per-individual F_ROH vs F_pedigree:\n")
  cat("F_ROH:      ", sprintf("%.3f", merged$froh), "\n")
  cat("F_pedigree: ", sprintf("%.3f", merged$f_pedigree), "\n")
  correlation <- cor(merged$froh, merged$f_pedigree)
  cat(sprintf("\nCorrelation(F_ROH, F_pedigree) at final generation: %.3f\n", correlation))
  cat(sprintf("Mean F_ROH: %.4f, Mean F_pedigree: %.4f, Mean |F_ROH - F_pedigree|: %.4f\n",
              mean(merged$froh), mean(merged$f_pedigree), mean(abs(merged$froh - merged$f_pedigree))))

  # Both should rise substantially from generation 5 to the final generation under this
  # small-population + intense-selection scheme.
  expect_true(froh_mean_trace[n_gen] > froh_mean_trace[5])
  expect_true(fped_mean_trace[n_gen] > fped_mean_trace[5])

  # They should correlate (same underlying inbreeding process) but not be identical (F_ROH
  # reflects realized recombination/Mendelian sampling variance, pedigree F does not).
  expect_true(correlation > 0.3)
  expect_true(mean(abs(merged$froh - merged$f_pedigree)) > 0.001)
})
