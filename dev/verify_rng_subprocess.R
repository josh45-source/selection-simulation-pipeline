args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

if (mode == "save") {
  state_path <- args[2]
  ref_path <- args[3]
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  invisible(runif(5)) # advance state a bit, like founder gen + a few operations would
  saveRDS(.Random.seed, state_path)
  ref_draws <- runif(3)
  saveRDS(ref_draws, ref_path)
  cat("Reference next-3 draws:", ref_draws, "\n")
} else if (mode == "restore") {
  state_path <- args[2]
  state <- readRDS(state_path)
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", state, envir = .GlobalEnv)
  draws <- runif(3)
  cat("Restored-in-subprocess next-3 draws:", draws, "\n")
}
