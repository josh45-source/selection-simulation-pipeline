# Serial generation loop for one replicate. Generation t+1 depends on t (serial within a
# replicate); replicates/scenarios are independent and safe to parallelize by the caller (run.R).
#
# Checkpointing and full RNG-state save/restore are layered on after this slice is proven
# end-to-end (see validation/test_resume_equals_uninterrupted.R once added).

library(AlphaSimR)

run_replicate <- function(config, selection_fraction, replicate_id, rng_stream) {
  set_replicate_stream(rng_stream)

  founders <- build_founder_pop(config)
  pop <- founders$pop
  SP <- founders$SP
  validation_chip <- founders$chips$validation

  n_gen      <- config$generations$n_generations
  n_founders <- config$population$n_founders
  h2         <- config$qtl$heritability
  mu         <- config$qtl$mutation_rate
  cadence    <- config$sampling$locus_cadence
  seg_only   <- isTRUE(config$sampling$segregating_only)

  summary_rows <- vector("list", n_gen)
  locus_rows <- list()

  for (gen in seq_len(n_gen)) {
    pop <- apply_recurrent_mutation(pop, SP, mu, validation_chip)
    pop <- setPheno(pop, h2 = h2, simParam = SP)

    n_parents <- max(1, round(n_founders * selection_fraction))
    parents <- selectInd(pop, nInd = n_parents, use = "pheno", simParam = SP)
    pop <- randCross(parents, nCrosses = n_founders, simParam = SP)

    summary_rows[[gen]] <- data.frame(
      generation = gen,
      va = compute_va(pop),
      mean_gv = meanG(pop)
    )

    if (gen %% cadence == 0) {
      geno <- pullQtlGeno(pop, simParam = SP)
      freq <- allele_freq(geno)
      site_ids <- colnames(geno)
      if (seg_only) {
        keep <- !is_fixed(freq)
        freq <- freq[keep]
        site_ids <- site_ids[keep]
      }
      if (length(freq) > 0) {
        locus_rows[[length(locus_rows) + 1]] <- data.frame(
          generation = gen, site_id = site_ids, freq = freq
        )
      }
    }
  }

  per_locus <- if (length(locus_rows) > 0) {
    do.call(rbind, locus_rows)
  } else {
    data.frame(generation = integer(0), site_id = character(0), freq = numeric(0))
  }

  list(
    summary = do.call(rbind, summary_rows),
    per_locus = per_locus,
    replicate_id = replicate_id
  )
}
