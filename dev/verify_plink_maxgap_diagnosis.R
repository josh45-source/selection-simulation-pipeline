library(AlphaSimR)
source("R/engine.R")
source("R/mutation.R")
source("R/froh.R")

config <- list(
  population = list(n_founders = 20, species = "generic_diploid"),
  genome = list(n_chr = 1, segsites_per_chr = 60),
  qtl = list(n_active_per_chr = 5, n_reservoir_per_chr = 15, heritability = 0.3,
             mean = 0, var_add = 1, mutation_rate = 2.3e-4),
  neutral_markers = list(n_active_per_chr = 5, n_reservoir_per_chr = 15, mutation_rate = 2.3e-4),
  validation_markers = list(n_per_chr = 20, mutation_rate = 0)
)
set.seed(1)
founders <- build_founder_pop(config)
pop <- founders$pop
SP <- founders$SP
chip <- founders$chips$validation

# Inbreed hard: 3 parents only, 20 generations
for (gen in 1:20) {
  pop <- apply_recurrent_mutation(pop, SP, config$qtl$mutation_rate, chip)
  pop <- setPheno(pop, h2 = 0.3, simParam = SP)
  parents <- selectInd(pop, nInd = 3, use = "pheno", simParam = SP)
  pop <- randCross(parents, nCrosses = 20, simParam = SP)
}

geno <- pullSnpGeno(pop, snpChip = chip, simParam = SP)
cat("Per-individual heterozygosity fraction after 20 gens of intense inbreeding:\n")
print(round(rowMeans(geno == 1), 3))

map <- getSnpMap(snpChip = chip, simParam = SP)
tmp_dir <- tempfile("inspect_")
dir.create(tmp_dir)
ped_path <- file.path(tmp_dir, "g.ped")
map_path <- file.path(tmp_dir, "g.map")
write_plink_files(pop, SP, chip, ped_path, map_path)

cat("\nFirst individual's genotype row (raw dosage):\n")
print(geno[1, ])
cat("\nFirst individual's .ped line:\n")
cat(readLines(ped_path, n = 1), "\n")

cat("\n--- consecutiveRUNS.run(), minSNP=5, minLengthBps=1000 ---\n")
runs <- consecutiveRUNS.run(genotypeFile = ped_path, mapFile = map_path,
                             ROHet = FALSE, minSNP = 5, minLengthBps = 1000)
cat("nrow(runs):", nrow(runs), "\n")
print(head(runs, 10))
