library(AlphaSimR)

founderPop <- runMacs(nInd = 40, nChr = 1, segSites = 50, species = "GENERIC", nThreads = 1L)
SP <- SimParam$new(founderPop)
SP$nThreads <- 1L
SP$addTraitA(50)
pop <- newPop(founderPop, simParam = SP)

ped <- data.frame(id = as.numeric(pop@id), mother = as.numeric(pop@mother), father = as.numeric(pop@father))
for (gen in 1:10) {
  pop <- setPheno(pop, h2 = 0.3, simParam = SP)
  parents <- selectInd(pop, nInd = 10, use = "pheno", simParam = SP)
  pop <- randCross(parents, nCrosses = 40, simParam = SP)
  ped <- rbind(ped, data.frame(id = as.numeric(pop@id), mother = as.numeric(pop@mother), father = as.numeric(pop@father)))
}

used_as_father <- unique(ped$father[ped$father != 0])
used_as_mother <- unique(ped$mother[ped$mother != 0])
conflict_ids <- intersect(used_as_father, used_as_mother)
cat("Individuals used as BOTH father (some row) and mother (another row):", length(conflict_ids), "\n")
if (length(conflict_ids) > 0) cat("Example ids:", head(conflict_ids, 10), "\n")
