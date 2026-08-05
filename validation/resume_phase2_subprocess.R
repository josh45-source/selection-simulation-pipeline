# Helper invoked as a genuinely separate Rscript process by
# validation/test_resume_equals_uninterrupted.R, so the resumed run can only succeed using what's
# actually in the checkpoint file -- not leftover in-memory state from the process that wrote it.
# args: <config_rds_path> <checkpoint_path> <result_out_path>

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1]
checkpoint_path <- args[2]
result_out_path <- args[3]

setwd(Sys.getenv("PROJECT_ROOT", unset = "."))
source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")
source("R/checkpoint.R")
source("R/io.R")
source("R/sim_loop.R")

config <- readRDS(config_path)

result <- run_replicate(
  config,
  selection_fraction = 0.2,
  replicate_id = 1,
  rng_stream = NULL, # ignored on resume -- checkpoint's RNG state is authoritative
  resume_path = checkpoint_path
)

saveRDS(
  list(
    summary = result$summary,
    per_locus = result$per_locus,
    qtl_geno = pullQtlGeno(result$pop, simParam = result$SP),
    marker_geno = pullSnpGeno(result$pop, snpChip = 1, simParam = result$SP),
    validation_geno = pullSnpGeno(result$pop, snpChip = 2, simParam = result$SP)
  ),
  result_out_path
)
