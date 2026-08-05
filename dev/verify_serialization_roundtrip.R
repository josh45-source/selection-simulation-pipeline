library(AlphaSimR)

set.seed(1)
founderPop <- runMacs(nInd = 20, nChr = 2, segSites = 100, species = "GENERIC")
SP <- SimParam$new(founderPop)
SP$addTraitA(100, mean = 0, var = 1)
pop <- newPop(founderPop, simParam = SP)

geno_before <- pullQtlGeno(pop, simParam = SP)

tmp <- tempfile(fileext = ".rds")
saveRDS(list(pop = pop, SP = SP), tmp)
rt <- readRDS(tmp)

geno_after <- pullQtlGeno(rt$pop, simParam = rt$SP)

cat("identical(geno_before, geno_after):", identical(geno_before, geno_after), "\n")
cat("Mismatched cells:", sum(geno_before != geno_after), "of", length(geno_before), "\n")

# Also check that continuing the simulation after round-trip behaves identically to continuing
# without round-tripping, GIVEN THE SAME starting .Random.seed for the next step.
state_before_step <- .Random.seed
pop2_direct <- setPheno(pop, h2 = 0.3, simParam = SP)

.Random.seed <<- state_before_step
pop2_roundtrip <- setPheno(rt$pop, h2 = 0.3, simParam = rt$SP)

cat("identical(pheno after setPheno, direct vs roundtripped pop/SP):",
    identical(pop2_direct@pheno, pop2_roundtrip@pheno), "\n")
