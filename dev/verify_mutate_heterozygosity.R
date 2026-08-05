library(AlphaSimR)

set.seed(1)
founderPop <- quickHaplo(nInd = 200, nChr = 1, segSites = 50)
SP <- SimParam$new(founderPop)
SP$addTraitA(50)
pop <- newPop(founderPop, simParam = SP)

geno_before <- pullSegSiteGeno(pop, simParam = SP)
cat("Genotype value counts before mutate():\n")
print(table(geno_before))

pop2 <- mutate(pop, mutRate = 0.05, simParam = SP)
geno_after <- pullSegSiteGeno(pop2, simParam = SP)

n_new_het <- sum((geno_before == 0 & geno_after == 1) | (geno_before == 2 & geno_after == 1))
n_hom_flipped <- sum(geno_before == 0 & geno_after == 2) + sum(geno_before == 2 & geno_after == 0)
n_het_to_hom <- sum(geno_before == 1 & geno_after != 1)
cat("NEW heterozygotes created (homozygous -> het, i.e. single-copy mutation):", n_new_het, "\n")
cat("Direct 0->2 or 2->0 (both-copy-at-once) flips after mutate():", n_hom_flipped, "\n")
cat("Existing hets that became homozygous (back-mutation on one copy):", n_het_to_hom, "\n")
cat("Total genotype cells that changed at all:", sum(geno_before != geno_after), "\n")

if (n_new_het > 0) {
  cat("CONCLUSION: mutate() produces single-haplotype (heterozygous-capable) mutations.\n")
} else {
  cat("CONCLUSION: mutate() did NOT produce any heterozygotes -- forces homozygosity like editGenome, or no mutations occurred.\n")
}
