# Independent, non-overlapping RNG streams per replicate via L'Ecuyer-CMRG
# (parallel::nextRNGStream), not just different integer seeds -- correlated replicates would
# quietly invalidate cross-replicate comparisons in the summary metrics table.
# ref: lecuyer2002streams, lecuyer1999cmrg
#
# Checkpoint/resume must restore RNG *state*, not re-seed (see validation/test_resume_equals_uninterrupted.R,
# added once checkpointing is wired into sim_loop.R).

make_replicate_streams <- function(base_seed, n_replicates) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)
  streams <- vector("list", n_replicates)
  seed <- .Random.seed
  for (i in seq_len(n_replicates)) {
    streams[[i]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }
  streams
}

set_replicate_stream <- function(stream_seed) {
  assign(".Random.seed", stream_seed, envir = .GlobalEnv)
}

save_rng_state <- function() {
  .Random.seed
}

restore_rng_state <- function(state) {
  assign(".Random.seed", state, envir = .GlobalEnv)
}
