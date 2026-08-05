library(AlphaSimR)

command <- "40 1E8 -t 1E-5 -r 4E-6  -eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0 -s "
segSites <- c(210L, 210L)
seed <- c("11111111", "22222222")

out1 <- AlphaSimR:::MaCS(command, segSites, FALSE, 2L, 1L, seed)
out2 <- AlphaSimR:::MaCS(command, segSites, FALSE, 2L, 1L, seed)

cat("identical(out1$geno, out2$geno):", identical(out1$geno, out2$geno), "\n")
if (!identical(out1$geno, out2$geno)) {
  g1 <- unlist(out1$geno)
  g2 <- unlist(out2$geno)
  cat("Mismatched values:", sum(g1 != g2), "of", length(g1), "\n")
}
cat("identical(out1$genMap, out2$genMap):", identical(out1$genMap, out2$genMap), "\n")
