#!/usr/bin/env Rscript
# Entrypoint: reads a scenario config, runs all replicates x selection intensities, writes
# Parquet, runs one DuckDB query over it, and produces one plot. This is the smallest
# end-to-end slice: sim -> Parquet -> DuckDB query -> plot.

library(yaml)
library(arrow)
library(duckdb)
library(DBI)
library(ggplot2)

source("R/engine.R")
source("R/mutation.R")
source("R/metrics.R")
source("R/rng.R")
source("R/io.R")
source("R/sim_loop.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "config/scenario.example.yaml"
config <- yaml::read_yaml(config_path)

n_replicates <- config$replicate$n_replicates
base_seed    <- config$replicate$seed
streams      <- make_replicate_streams(base_seed, n_replicates)

scenario_id <- config$scenario_id
data_dir    <- config$output$data_dir
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

for (rep_id in seq_len(n_replicates)) {
  for (intensity_name in names(config$selection$intensities)) {
    frac <- config$selection$intensities[[intensity_name]]
    cat(sprintf("Replicate %d, selection=%s (fraction=%.2f)\n", rep_id, intensity_name, frac))

    result <- run_replicate(config, frac, rep_id, streams[[rep_id]])
    scen_tag <- paste0(scenario_id, "_", intensity_name)

    write_partitioned_parquet(result$summary, data_dir, "summary_metrics",
                               scenario_id = scen_tag, replicate_id = rep_id)
    write_partitioned_parquet(result$per_locus, data_dir, "per_locus_freq",
                               scenario_id = scen_tag, replicate_id = rep_id)
  }
}

con <- dbConnect(duckdb::duckdb())
query <- sprintf(
  "SELECT scenario_id, generation, AVG(va) AS mean_va
   FROM read_parquet('%s/summary_metrics/**/*.parquet', hive_partitioning=1)
   GROUP BY scenario_id, generation
   ORDER BY scenario_id, generation",
  data_dir
)
va_by_gen <- dbGetQuery(con, query)
dbDisconnect(con, shutdown = TRUE)

p <- ggplot(va_by_gen, aes(x = generation, y = mean_va, color = scenario_id)) +
  geom_line() +
  labs(title = "Mean additive genetic variance across replicates",
       x = "Generation", y = "Mean Va") +
  theme_minimal()

plot_path <- file.path(data_dir, "va_trajectory.png")
ggsave(plot_path, p, width = 8, height = 5)
cat("Wrote plot to", plot_path, "\n")
