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
