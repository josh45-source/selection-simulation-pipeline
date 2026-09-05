# Thin interface over AlphaSimR for founder population + genetic architecture setup.
# Kept separate from sim_loop.R/mutation.R so a future MoBPS backend (better native support for
# recurrent mutation and OCS) can be substituted here without touching the rest of the pipeline.

library(AlphaSimR)

# Builds the founder population and SimParam from a parsed scenario config (see
# config/scenario.example.yaml). Three site classes are designated at init:
#   - qtl: drives the trait; active + reservoir sites both get effects assigned by addTraitA,
#     since AlphaSimR only assigns effects to sites designated at init (see scenario yaml comments
#     on the active/reservoir split).
#   - neutral_markers (SNP chip 1): mutates at the same rate as qtl (see R/mutation.R for why).
#   - validation_markers (SNP chip 2): held at mu = 0 by R/mutation.R; used only for the
#     neutral heterozygosity-decay validation test and as the F_ROH panel, so it must stay clean
#     of mutation-induced homozygosity artifacts.
#
# These classes are NOT mutually exclusive, and an earlier version of this comment wrongly said
# they were. QTL are excluded from SNP eligibility, so neither chip overlaps the QTL -- that part
# holds. But addSnpChip() "randomly assigns eligible SNPs to a SNP chip" independently on each
# call, and nothing here calls SimParam$restrSegSites(overlap=), so the two chips are each drawn
# from the same pool of non-QTL sites and do overlap heavily. Measured on scenario.slice.yaml:
# 134 of chip 2's 150 sites are also on chip 1. That is close to forced by the site budget --
# QTL take 500 of 1050 sites per chromosome and chip 1 then takes 500 of the 550 that remain.
#
# The segsites_per_chr check below therefore budgets for a partition that AlphaSimR does not
# actually produce. It is kept because it still guarantees enough segregating sites exist for
# every class to be drawn, but it should not be read as enforcing disjointness.
#
# This is documented rather than fixed because it affects no reported metric. Chip 1 is never
# used for a scientific result -- it appears only as a bit-comparison target in
# validation/test_resume_equals_uninterrupted.R and validation/resume_phase2_subprocess.R.
# Everything reported comes from the QTL (additive variance, response, fixation) or from chip 2
# (heterozygosity decay, F_ROH), and chip 2 is correctly neutral and correctly frozen. Fixing the
# allocation via restrSegSites(excludeSnp=) or addSnpChipByName() would re-sample every site and
# change every published number, for no gain in any result that is actually reported.
build_founder_pop <- function(config) {
  n_founders <- config$population$n_founders
  n_chr <- config$genome$n_chr

  qtl_total        <- config$qtl$n_active_per_chr + config$qtl$n_reservoir_per_chr
  marker_total     <- config$neutral_markers$n_active_per_chr + config$neutral_markers$n_reservoir_per_chr
  validation_total <- config$validation_markers$n_per_chr

  segsites_per_chr <- config$genome$segsites_per_chr
  expected_total <- qtl_total + marker_total + validation_total
  if (segsites_per_chr != expected_total) {
    stop(sprintf(
      "genome$segsites_per_chr (%d) must equal qtl + neutral_markers + validation_markers totals per chr (%d).",
      segsites_per_chr, expected_total
    ))
  }

  founderPop <- runMacs(
    nInd = n_founders,
    nChr = n_chr,
    segSites = segsites_per_chr,
    species = "GENERIC",
    nThreads = 1L # separate from SimParam$nThreads below -- runMacs has its own thread count,
                  # defaulting to auto-detected (multi-threaded), which breaks reproducibility
                  # from a seed just as much as SimParam's did (see R/engine.R's other comment).
  )

  SP <- SimParam$new(founderPop)
  # Force single-threaded execution. SimParam$nThreads defaults to all available cores via OpenMP,
  # and AlphaSimR's own documentation examples all hide a `SP$nThreads = 1L` line -- an implicit
  # acknowledgment that multi-threaded execution isn't bit-reproducible from a seed (thread
  # scheduling order isn't deterministic). Confirmed empirically: with the default thread count,
  # resuming from a checkpoint diverged substantially from an uninterrupted run despite correct RNG
  # state restoration (see validation/test_resume_equals_uninterrupted.R); setting this to 1L
  # before any AlphaSimR operations use SP fixed it. This is required for reproducibility in
  # general, not just checkpoint/resume.
  SP$nThreads <- 1L
  SP$addTraitA(qtl_total, mean = config$qtl$mean, var = config$qtl$var_add)
  SP$addSnpChip(marker_total)      # chip 1: mutating neutral markers
  SP$addSnpChip(validation_total) # chip 2: validation markers, kept at mu = 0

  pop <- newPop(founderPop, simParam = SP)

  list(
    pop = pop,
    SP = SP,
    chips = list(markers = 1L, validation = 2L)
  )
}
