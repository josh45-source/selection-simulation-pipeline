library(detectRUNS)

genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")

map_dt <- data.table::fread(mapFile, header = FALSE)
colnames(map_dt) <- c("Chrom", "SNP", "cM", "bps")
cat("Chromosomes present in official example map:", unique(map_dt$Chrom), "\n")
cat("Count per chromosome:\n")
print(table(map_dt$Chrom))

# Subset both ped and map to ONLY chromosome 2 (mimicking my single-chromosome test case)
keep_snps <- map_dt$SNP[map_dt$Chrom == 2]
cat("\nN SNPs on chrom 2:", length(keep_snps), "\n")

# Build a filtered map file
map_sub <- map_dt[map_dt$Chrom == 2, ]
map_sub_path <- tempfile(fileext = ".map")
write.table(map_sub, map_sub_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Build a filtered ped: keep first 6 cols + only the allele-pairs for chrom2 SNPs, in map order
all_lines <- readLines(genotypeFile)
first <- strsplit(all_lines[1], " ")[[1]]
n_markers_total <- (length(first) - 6) / 2
cat("Total markers in full ped:", n_markers_total, "\n")
snp_idx <- which(map_dt$Chrom == 2) # column positions correspond to map row order

ped_sub_path <- tempfile(fileext = ".ped")
con_out <- file(ped_sub_path, "w")
for (ln in all_lines) {
  fields <- strsplit(ln, " ")[[1]]
  header6 <- fields[1:6]
  alleles <- fields[7:length(fields)]
  keep_cols <- as.vector(rbind(2 * snp_idx - 1, 2 * snp_idx))
  new_line <- paste(c(header6, alleles[keep_cols]), collapse = " ")
  writeLines(new_line, con_out)
}
close(con_out)

cat("\n--- Running consecutiveRUNS.run() on single-chromosome subset of the OFFICIAL data ---\n")
runs_sub <- consecutiveRUNS.run(ped_sub_path, map_sub_path, minSNP = 15, ROHet = FALSE,
                                 maxOppRun = 0, maxMissRun = 0, maxGap = 10^6,
                                 minLengthBps = 100000)
cat("nrow(runs) on single-chromosome subset:", nrow(runs_sub), "\n")
