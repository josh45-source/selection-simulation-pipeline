# Re-run of dev/verify_macs_direct.R against a *specified* AlphaSimR build, to check whether
# PR #265 actually fixed MaCS's seed non-determinism rather than assuming it did from the
# GitHub thread alone. Deliberately run with NO project directory mounted and NO renv
# activation, to rule out contamination from this project's own pinned 1.5.3 library.
#
# The branch trap this script exists to avoid: PR #265 was merged to AlphaSimR's `devel`
# branch (2026-04-03), NOT to its default branch. `remotes::install_github("gaynorr/AlphaSimR")`
# installs the default branch (`master`), which is the 2.1.0 release and does NOT contain the
# fix. Installing the fix requires ref = "devel" explicitly. Verified results:
#
#   AlphaSimR 1.5.3   (CRAN, the version renv.lock pins)          -> NOT reproducible
#   AlphaSimR 2.1.0   (CRAN 2025-11-08 / master, no PR #265)      -> NOT reproducible
#   AlphaSimR 2.1.0.9004 (devel @ e667a37, contains PR #265)      -> reproducible
#
# PR #265 also changed MaCS()'s seed argument from character to numeric (runMacs() now draws
# seeds with sample.int(), so R's set.seed() reaches them). This script picks the seed type
# from the loaded version so it runs against both old and new builds.
#
# Install the fixed build into an isolated library:
#   remotes::install_github("gaynorr/AlphaSimR", ref = "devel", lib = "<isolated-lib>")

cat("Library paths in use:", paste(.libPaths(), collapse = " | "), "\n")
cat("AlphaSimR loaded from:", find.package("AlphaSimR"), "\n")
desc <- packageDescription("AlphaSimR")
cat("Version field:", desc$Version, "\n")
cat("RemoteRef (branch installed from, if remotes tracked it):",
    if (!is.null(desc$RemoteRef)) desc$RemoteRef else "NOT RECORDED", "\n")
cat("RemoteSha (git commit installed from, if remotes tracked it):",
    if (!is.null(desc$RemoteSha)) desc$RemoteSha else "NOT RECORDED", "\n")
cat("Packaged/Built:", if (!is.null(desc$Packaged)) desc$Packaged else desc$Built, "\n\n")

library(AlphaSimR)

# Post-#265 builds take a numeric seed; 1.5.3 and 2.1.0 take character.
fixed_seed_api <- utils::compareVersion(as.character(packageVersion("AlphaSimR")), "2.1.0") > 0
seed <- if (fixed_seed_api) c(11111111, 22222222) else c("11111111", "22222222")
cat("Seed argument passed as:", class(seed), "\n\n")

command <- "40 1E8 -t 1E-5 -r 4E-6  -eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0 -s "
segSites <- c(210L, 210L)

cat("--- direct AlphaSimR:::MaCS(), identical literal seeds ---\n")
out1 <- AlphaSimR:::MaCS(command, segSites, FALSE, 2L, 1L, seed)
out2 <- AlphaSimR:::MaCS(command, segSites, FALSE, 2L, 1L, seed)
cat("identical(out1$geno, out2$geno):", identical(out1$geno, out2$geno), "\n")
if (!identical(out1$geno, out2$geno)) {
  g1 <- unlist(out1$geno)
  g2 <- unlist(out2$geno)
  cat("Mismatched values:", sum(g1 != g2), "of", length(g1), "\n")
}
cat("identical(out1$genMap, out2$genMap):", identical(out1$genMap, out2$genMap), "\n")

cat("\n--- issue #228's own repro (set.seed + runMacs, nThreads=1) ---\n")
set.seed(123); result1 <- runMacs(nInd = 100L, nChr = 1L, segSites = 50L, nThreads = 1L)
set.seed(123); result2 <- runMacs(nInd = 100L, nChr = 1L, segSites = 50L, nThreads = 1L)
cat("identical(result1, result2):", identical(result1, result2), "\n")

cat("\n--- multi-chromosome / multi-threaded under set.seed ---\n")
set.seed(456); m1 <- runMacs(nInd = 50L, nChr = 2L, segSites = 100L, nThreads = 2L)
set.seed(456); m2 <- runMacs(nInd = 50L, nChr = 2L, segSites = 100L, nThreads = 2L)
cat("identical(m1, m2):", identical(m1, m2), "\n")

# Negative control: a build that returned a constant would pass every test above.
cat("\n--- negative control: different seeds must still differ ---\n")
set.seed(1); n1 <- runMacs(nInd = 50L, nChr = 1L, segSites = 50L, nThreads = 1L)
set.seed(2); n2 <- runMacs(nInd = 50L, nChr = 1L, segSites = 50L, nThreads = 1L)
cat("identical(n1, n2) [expect FALSE]:", identical(n1, n2), "\n")
