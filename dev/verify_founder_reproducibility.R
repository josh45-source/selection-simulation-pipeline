source("R/engine.R")
source("R/rng.R")

config <- list(
  population = list(n_founders = 20, species = "generic_diploid"),
  genome = list(n_chr = 2, segsites_per_chr = 210),
  qtl = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, heritability = 0.3,
             mean = 0, var_add = 1, mutation_rate = 2.3e-4),
  neutral_markers = list(n_active_per_chr = 20, n_reservoir_per_chr = 80, mutation_rate = 2.3e-4),
  validation_markers = list(n_per_chr = 10, mutation_rate = 0)
)

RNGkind("L'Ecuyer-CMRG")
set.seed(909090)
stream <- .Random.seed

set_replicate_stream(stream)
f1 <- build_founder_pop(config)
g1 <- pullQtlGeno(f1$pop, simParam = f1$SP)

set_replicate_stream(stream)
f2 <- build_founder_pop(config)
g2 <- pullQtlGeno(f2$pop, simParam = f2$SP)

cat("identical(founder geno 1, founder geno 2):", identical(g1, g2), "\n")
cat("Mismatched cells:", sum(g1 != g2), "of", length(g1), "\n")
