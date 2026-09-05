# Proves checkpoint/resume restores full RNG *state* (not just re-seeding) alongside the
# population object, so a resumed run is bit-for-bit identical to continuing without interruption
# FROM THE SAME STARTING POPULATION. This is the trust anchor for using checkpointing on real
# long runs at all.
#
# Important scope note: this test does NOT (and cannot currently) prove that founder generation
# itself is reproducible from a scenario config + seed. AlphaSimR's compiled coalescent function
# (AlphaSimR:::MaCS, called by runMacs) was confirmed empirically (dev/verify_macs_direct.R) to
# produce DIFFERENT output across calls given IDENTICAL literal seeds -- a non-determinism baked
# into AlphaSimR's compiled code, not something fixable from this project. Both branches below
# therefore start from a single, already-generated founder population loaded from one checkpoint
# file, rather than each calling build_founder_pop() independently -- that's the actual guarantee
# checkpoint/resume needs to provide for long runs, and it's what this test isolates.
#
# The resumed half runs in a genuinely separate Rscript process (validation/resume_phase2_subprocess.R),
# not just a second function call in this process -- otherwise leftover in-memory state (the
# global .Random.seed, loaded objects) could mask a checkpoint that doesn't actually capture
# everything it needs to.
#
# RUN THIS THE DOCUMENTED WAY -- against the image's own baked-in code and library:
#
#   docker run --rm --user "$(id -u):$(id -g)" selection:latest validation/test_resume_equals_uninterrupted.R
#
# Do NOT run it from a mounted host path. This test is the one script in the suite that is
# sensitive to how it is invoked, because of the system2("Rscript", ...) call below: that
# subprocess is launched WITHOUT --vanilla, so it reads whatever .Rprofile sits in its working
# directory. Inside the image that is the intended one. With the repository mounted over the
# working directory instead, the subprocess picks up the repo's .Rprofile, activates renv against
# the HOST renv/library rather than the image's, and dies before writing result.rds -- surfacing
# here as a confusing triple failure (non-zero exit status, missing result.rds, then a readRDS
# connection error) that looks like a checkpoint bug and is not one. The host library is also not
# guaranteed to match renv.lock. If you need to test modified sources, either rebuild the image or
# overlay only the source directories, leaving the image's /project/renv in place:
#
#   docker run --rm -v "$PWD/R:/project/R:ro" -v "$PWD/validation:/project/validation:ro" \
#     selection:latest validation/test_resume_equals_uninterrupted.R

library(testthat)

source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")
source("R/checkpoint.R")
source("R/io.R")
source("R/sim_loop.R")

test_that("resume from a checkpoint reproduces an uninterrupted run bit-for-bit", {
  # Small scenario: just needs to exercise founder gen, mutation, selection, crossing, and a
  # checkpoint boundary -- not be scientifically meaningful on its own.
  config <- list(
    scenario_id = "resume_test",
    replicate = list(n_replicates = 1, seed = 909090),
    population = list(n_founders = 20, species = "generic_diploid"),
    genome = list(n_chr = 2, segsites_per_chr = 210),
    qtl = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, heritability = 0.3,
               mean = 0, var_add = 1, mutation_rate = 2.3e-4),
    neutral_markers = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, mutation_rate = 2.3e-4),
    validation_markers = list(n_per_chr = 10, mutation_rate = 0),
    generations = list(n_generations = 20, checkpoint_every = 10),
    sampling = list(locus_cadence = 5, segregating_only = TRUE)
  )
  expected_segsites <- (config$qtl$n_active_per_chr + config$qtl$n_reservoir_per_chr) +
    (config$neutral_markers$n_active_per_chr + config$neutral_markers$n_reservoir_per_chr) +
    config$validation_markers$n_per_chr
  stopifnot(config$genome$segsites_per_chr == expected_segsites)

  stream <- make_replicate_streams(config$replicate$seed, 1)[[1]]

  tmp_dir <- tempfile("resume_test_")
  dir.create(tmp_dir)

  # Build founders ONCE. Both branches below resume from this same checkpoint rather than each
  # calling build_founder_pop() independently -- see the file-level comment on why.
  set_replicate_stream(stream)
  founders <- build_founder_pop(config)
  gen0_path <- file.path(tmp_dir, "gen0.rds")
  save_checkpoint(gen0_path, founders$pop, founders$SP, founders$chips, gen = 0L,
                   summary_rows = vector("list", config$generations$n_generations),
                   locus_rows = list())

  # Reference: resume from gen 0, run straight through to generation 20 in one call.
  ref <- run_replicate(config, selection_fraction = 0.2, replicate_id = 1,
                        rng_stream = NULL, resume_path = gen0_path)

  # Phase 1: SEPARATELY resume from the same gen 0 checkpoint, run only to generation 10,
  # checkpoint there.
  checkpoint_path <- file.path(tmp_dir, "checkpoint.rds")
  config_phase1 <- config
  config_phase1$generations$n_generations <- 10L
  invisible(run_replicate(config_phase1, selection_fraction = 0.2, replicate_id = 1,
                           rng_stream = NULL, resume_path = gen0_path,
                           checkpoint_path = checkpoint_path))
  expect_true(file.exists(checkpoint_path))

  # Phase 2: resume in a genuinely separate Rscript process, continuing gens 11-20.
  config_path <- file.path(tmp_dir, "config.rds")
  result_path <- file.path(tmp_dir, "result.rds")
  saveRDS(config, config_path)
  project_root <- normalizePath(".")

  out <- system2(
    "Rscript",
    args = c("validation/resume_phase2_subprocess.R", config_path, checkpoint_path, result_path),
    env = paste0("PROJECT_ROOT=", project_root),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(paste(out, collapse = "\n"), "\n")
  }
  expect_true(is.null(status) || status == 0)
  expect_true(file.exists(result_path))

  resumed <- readRDS(result_path)

  ref_qtl <- pullQtlGeno(ref$pop, simParam = ref$SP)
  ref_markers <- pullSnpGeno(ref$pop, snpChip = 1, simParam = ref$SP)
  ref_validation <- pullSnpGeno(ref$pop, snpChip = 2, simParam = ref$SP)

  expect_identical(ref_qtl, resumed$qtl_geno)
  expect_identical(ref_markers, resumed$marker_geno)
  expect_identical(ref_validation, resumed$validation_geno)

  # Full recorded output series must match exactly too, not just the final genotypes.
  expect_identical(ref$summary, resumed$summary)
  expect_identical(ref$per_locus, resumed$per_locus)

  unlink(tmp_dir, recursive = TRUE)
})
