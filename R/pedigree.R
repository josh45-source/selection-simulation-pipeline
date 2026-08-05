# Pedigree-based inbreeding coefficient via the classic recursive tabular method (kinship matrix),
# for comparison against F_ROH. Pedigree F is the EXPECTED inbreeding coefficient from the
# relationship structure; realized (genomic/ROH-based) F varies around it due to Mendelian
# sampling variance in actual recombination -- they should correlate but are not the same quantity.
# ref: kardos2015genomic
#
# Implemented directly (not via kinship2) because AlphaSimR is monoecious -- the same individual
# can be recorded as "mother" in one cross and "father" in another, which kinship2::pedigree()
# rejects outright (it enforces a strict two-sex model and validates dadid/momid role consistency
# per individual across the whole pedigree). Confirmed empirically this is common, not an edge
# case (dev/check_role_conflict.R: 79 conflicting individuals in a 10-generation test). The
# tabular method below only uses mother/father ID links and has no sex model, so it doesn't hit
# this mismatch at all.

# ped_df: data.frame with columns id, mother, father accumulated across all generations
# (founders' mother/father should be 0 or NA), ordered with parents before offspring (true for
# AlphaSimR's sequential id assignment).
pedigree_inbreeding <- function(ped_df) {
  ped_df$mother[is.na(ped_df$mother) | ped_df$mother == 0] <- NA
  ped_df$father[is.na(ped_df$father) | ped_df$father == 0] <- NA
  ids <- ped_df$id
  n <- length(ids)
  id_to_idx <- setNames(seq_len(n), as.character(ids))

  mother_idx <- unname(id_to_idx[as.character(ped_df$mother)])
  father_idx <- unname(id_to_idx[as.character(ped_df$father)])

  K <- matrix(0, n, n)
  for (i in seq_len(n)) {
    mi <- mother_idx[i]
    fi <- father_idx[i]
    if (i > 1) {
      row <- rep(0, i - 1)
      if (!is.na(mi)) row <- row + 0.5 * K[mi, 1:(i - 1)]
      if (!is.na(fi)) row <- row + 0.5 * K[fi, 1:(i - 1)]
      K[i, 1:(i - 1)] <- row
      K[1:(i - 1), i] <- row
    }
    K[i, i] <- if (!is.na(mi) && !is.na(fi)) 0.5 * (1 + K[mi, fi]) else 0.5
  }
  f <- 2 * diag(K) - 1
  data.frame(id = ids, f_pedigree = f)
}
