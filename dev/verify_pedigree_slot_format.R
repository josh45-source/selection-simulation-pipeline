library(AlphaSimR)

founderPop <- runMacs(nInd = 10, nChr = 1, segSites = 50, species = "GENERIC", nThreads = 1L)
SP <- SimParam$new(founderPop)
SP$nThreads <- 1L
SP$addTraitA(50)
pop <- newPop(founderPop, simParam = SP)

cat("Founder id class:", class(pop@id), "\n")
cat("Founder id:", pop@id, "\n")
cat("Founder mother:", pop@mother, "\n")
cat("Founder father:", pop@father, "\n")

pop2 <- randCross(pop, nCrosses = 5, simParam = SP)
cat("\nGen1 id:", pop2@id, "\n")
cat("Gen1 mother:", pop2@mother, "\n")
cat("Gen1 father:", pop2@father, "\n")
