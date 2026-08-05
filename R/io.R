# Parquet writers, partitioned by scenario_id/replicate_id per the data model in the design doc:
# summary metrics (narrow, one row per generation x replicate), per-locus allele frequencies (the
# firehose -- segregating sites only, sampled per cadence), and provenance (static per run).

library(arrow)
library(dplyr) # write_dataset() uses dplyr internally for character-vector `partitioning=`;
                # renv's implicit snapshot only picks up packages referenced in source, so this
                # needs to be an explicit library() call to land in renv.lock.

write_partitioned_parquet <- function(df, data_dir, dataset_name, scenario_id, replicate_id) {
  df$scenario_id <- scenario_id
  df$replicate_id <- replicate_id
  path <- file.path(data_dir, dataset_name)
  write_dataset(
    df,
    path = path,
    partitioning = c("scenario_id", "replicate_id"),
    format = "parquet",
    existing_data_behavior = "overwrite"
  )
}
