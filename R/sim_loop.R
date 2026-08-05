# Serial generation loop for one replicate. Generation t+1 depends on t (serial within a
# replicate); replicates/scenarios are independent and safe to parallelize by the caller (run.R).
#
# checkpoint_path: if given, saveRDS()s pop + SP + accumulated output rows + RNG state every
# config$generations$checkpoint_every generations (see R/checkpoint.R).
# resume_path: if given, restores from that checkpoint (pop, SP, RNG state, and rows accumulated
# so far) and continues from checkpoint_gen + 1, instead of building a fresh founder population.
# When resuming, rng_stream is ignored -- the checkpoint's restored RNG state IS the stream state.

library(AlphaSimR)

run_replicate <- function(config, selection_fraction, replicate_id, rng_stream,
                           checkpoint_path = NULL, resume_path = NULL) {
  n_gen      <- config$generations$n_generations
  n_founders <- config$population$n_founders
  h2         <- config$qtl$heritability
  mu         <- config$qtl$mutation_rate
  cadence    <- config$sampling$locus_cadence
  seg_only   <- isTRUE(config$sampling$segregating_only)
  checkpoint_every <- config$generations$checkpoint_every

  if (!is.null(resume_path)) {
    ckpt <- load_checkpoint(resume_path) # restores .Random.seed as a side effect
    pop <- ckpt$pop
    SP <- ckpt$SP
    validation_chip <- ckpt$chips$validation
    start_gen <- ckpt$gen + 1L
    summary_rows <- ckpt$summary_rows
    locus_rows <- ckpt$locus_rows
  } else {
    set_replicate_stream(rng_stream)
    founders <- build_founder_pop(config)
    pop <- founders$pop
    SP <- founders$SP
    validation_chip <- founders$chips$validation
    start_gen <- 1L
    summary_rows <- vector("list", n_gen)
    locus_rows <- list()
  }

  if (start_gen <= n_gen) {
    for (gen in start_gen:n_gen) {
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

      if (!is.null(checkpoint_path) && gen %% checkpoint_every == 0) {
        chips <- list(markers = 1L, validation = validation_chip)
        save_checkpoint(checkpoint_path, pop, SP, chips, gen, summary_rows, locus_rows)
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
    replicate_id = replicate_id,
    pop = pop,
    SP = SP
  )
}
