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
  # RNGkind() must be set to L'Ecuyer-CMRG BEFORE assigning the saved .Random.seed vector. A fresh
  # R session (e.g. the resume subprocess) starts on the default Mersenne-Twister generator;
  # directly overwriting .Random.seed with a foreign-kind state vector doesn't reliably switch the
  # session's RNG dispatch to match it, even though the bytes look right. Confirmed empirically --
  # omitting this line caused resumed runs to silently diverge from the reference run despite the
  # checkpoint file containing the "correct" state (see validation/test_resume_equals_uninterrupted.R).
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", state, envir = .GlobalEnv)
}
