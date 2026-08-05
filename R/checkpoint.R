# Checkpoint/resume: periodically persist the population object, its SimParam (needed to
# interpret the population's raw genotypes -- QTL effects, SNP chip site sets, genetic map -- so
# checkpointing pop alone isn't enough), the chip-index map, the accumulated per-generation output
# rows (otherwise an interruption loses everything recorded before the last checkpoint), and the
# full RNG *state* (.Random.seed, not a seed integer) together. Resuming must restore that RNG
# state rather than re-seed, so a resumed run is bit-for-bit identical to an uninterrupted one --
# see validation/test_resume_equals_uninterrupted.R.

save_checkpoint <- function(path, pop, SP, chips, gen, summary_rows, locus_rows) {
  saveRDS(
    list(
      pop = pop,
      SP = SP,
      chips = chips,
      gen = gen,
      rng_state = save_rng_state(),
      summary_rows = summary_rows,
      locus_rows = locus_rows
    ),
    path
  )
}

load_checkpoint <- function(path) {
  ckpt <- readRDS(path)
  restore_rng_state(ckpt$rng_state)
  ckpt
}
