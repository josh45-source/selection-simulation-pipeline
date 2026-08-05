source("R/rng.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(909090)
stream <- .Random.seed

set_replicate_stream(stream)
seed1 <- sapply(1:2, function(x) as.character(sample.int(1e+08, 1)))

set_replicate_stream(stream)
seed2 <- sapply(1:2, function(x) as.character(sample.int(1e+08, 1)))

cat("seed1:", seed1, "\n")
cat("seed2:", seed2, "\n")
cat("identical:", identical(seed1, seed2), "\n")
