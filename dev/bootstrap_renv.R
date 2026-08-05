# One-time bootstrap: initializes renv in the project, installs the pinned package set via renv
# (so it resolves through CRAN_REPO, the pinned PPM jammy snapshot), and snapshots renv.lock.
# Run once inside a container with the project mounted; renv.lock + renv/activate.R + .Rprofile +
# renv/settings.json are then committed, and the Dockerfile switches to renv::restore().

repo <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = repo), timeout = 600)

install.packages("renv")
renv::init(bare = TRUE)

pkgs <- c(
  "AlphaSimR", "arrow", "duckdb", "yaml", "dplyr", "tibble", "purrr",
  "detectRUNS", # covers F_ROH directly; optiSel omitted -- unresolvable from the pinned repo
                # snapshot for R 4.3.3 (not its own R floor, which is only >=3.5.0; more likely a
                # missing binary/transitive dep at this R-version+platform combo). Revisit if OCS
                # features specifically need it later.
  "ggplot2", "testthat"
)

for (attempt in 1:8) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) break
  if (attempt > 1) {
    wait <- min(60, 10 * attempt)
    cat(sprintf("Retry %d/8 for: %s (waiting %ds)\n", attempt, paste(missing, collapse = ", "), wait))
    Sys.sleep(wait)
  }
  try(renv::install(missing))
}
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Failed to install via renv: ", paste(missing, collapse = ", "))

renv::snapshot(prompt = FALSE)
cat("renv.lock written.\n")
