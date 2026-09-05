# What this script proves, and what it cannot see.
#
# It confirms the two things the pipeline depends on: mutate() does place mutations on chip 1, and
# the snapshot/restore brings chip 2 back bit-identical so the validation panel is genuinely held
# at mu = 0. Both still hold.
#
# It CANNOT detect partial freezing of chip 1, and the PASS below should not be read as proving
# its absence. The test only asks whether chip 1 changed *somewhere* (sum(...) > 0), so it passes
# whether the restore rolls back none of chip 1's mutations or most of them. AlphaSimR draws the
# two chips independently from the same pool of non-QTL sites, so they overlap, and restoring
# chip 2 does silently revert mutations at the sites they share -- about 8% of chip 1's mutations
# at the slice config's rate. Detecting that requires comparing chip 1 immediately after mutate()
# against chip 1 after the restore, which this script does not do. See R/engine.R and R/mutation.R
# for the measurement and for why it is documented rather than fixed.

library(AlphaSimR)

set.seed(42)
founderPop <- runMacs(nInd = 100, nChr = 3, segSites = 300, species = "GENERIC")
SP <- SimParam$new(founderPop)
SP$addTraitA(100)          # QTL, not touched by this test
SP$addSnpChip(100)         # chip 1: "mutating markers"
SP$addSnpChip(100)         # chip 2: "validation markers" -- must stay untouched

pop <- newPop(founderPop, simParam = SP)

mutating_before   <- pullSnpHaplo(pop, snpChip = 1, haplo = "all", simParam = SP)
validation_before <- pullSnpHaplo(pop, snpChip = 2, haplo = "all", simParam = SP)

pop2 <- mutate(pop, mutRate = 0.05, simParam = SP)   # high rate, deliberately genome-wide
pop2 <- setMarkerHaplo(pop2, haplo = validation_before, simParam = SP)

mutating_after   <- pullSnpHaplo(pop2, snpChip = 1, haplo = "all", simParam = SP)
validation_after <- pullSnpHaplo(pop2, snpChip = 2, haplo = "all", simParam = SP)

cat("Mutating-chip cells changed by mutate() (should be > 0):",
    sum(mutating_before != mutating_after), "\n")
cat("Validation-chip cells changed after revert (must be exactly 0):",
    sum(validation_before != validation_after), "\n")

if (sum(mutating_before != mutating_after) > 0 && sum(validation_before != validation_after) == 0) {
  cat("CONCLUSION: snapshot/restore via pullSnpHaplo + setMarkerHaplo correctly protects one chip while leaving mutation active elsewhere -- PASS.\n")
} else {
  cat("CONCLUSION: revert mechanism did NOT behave as expected -- FAIL, needs a different approach.\n")
}
