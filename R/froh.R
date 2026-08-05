# F_ROH (genomic inbreeding from runs of homozygosity) via detectRUNS.
# ref: curik2014roh, keller2011roh -- ROH as the genomic inbreeding estimator, over pedigree-only F
#
# detectRUNS only accepts PLINK .ped/.map files (no in-memory data frame interface), so genotypes
# are converted from AlphaSimR's dosage (0/1/2) matrices. Allele labels ("1"/"2") are arbitrary --
# detectRUNS only needs per-locus zygosity, not true allele identity or haplotype phase, since it
# detects consecutive-homozygous-site runs along the physical map, not phased haplotype blocks.
# Physical position (bp) is derived from AlphaSimR's genetic map position (Morgans) using the
# GENERIC species model's 1 Morgan : 100 Mb ratio (genLen=1, physical length 1E8bp -- see
# AlphaSimR::runMacs's GENERIC species parameters).

library(detectRUNS)

write_plink_files <- function(pop, SP, chip, ped_path, map_path) {
  geno <- pullSnpGeno(pop, snpChip = chip, simParam = SP)
  map <- getSnpMap(snpChip = chip, simParam = SP)
  stopifnot(identical(colnames(geno), map$id))

  map_out <- data.frame(
    chrom = map$chr,
    id = map$id,
    cM = map$pos * 100, # AlphaSimR pos is in Morgans; PLINK .map wants cM. Non-zero/monotonic --
                         # a constant 0 here was part of what caused detectRUNS to silently find
                         # no runs regardless of genotype content (see dev/inspect_plink_files.R).
    bp = round(map$pos * 1e8)
  )
  # Physical order must be strictly increasing within a chromosome for run detection to be
  # meaningful; AlphaSimR site naming (chr_index) already reflects map order.
  write.table(map_out, map_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  # Nucleotide-letter alleles, not numeric 1/2 -- matching detectRUNS' own example data format.
  allele1 <- ifelse(geno == 0, "A", "G")
  allele2 <- ifelse(geno == 2, "G", "A")
  allele_pairs <- matrix(NA_character_, nrow = nrow(geno), ncol = 2 * ncol(geno))
  allele_pairs[, seq(1, 2 * ncol(geno), by = 2)] <- allele1
  allele_pairs[, seq(2, 2 * ncol(geno), by = 2)] <- allele2

  ped_out <- cbind(
    famid = 1,
    id = rownames(geno),
    patid = 0,
    matid = 0,
    sex = 0,
    pheno = -9,
    as.data.frame(allele_pairs)
  )
  write.table(ped_out, ped_path, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
}

# Returns a data.frame: id, froh (genome-wide F_ROH).
#
# maxGap ("max distance between consecutive SNPs to still count as part of a run") defaults to
# 10^6 (1Mb) in detectRUNS, calibrated for real SNP-chip density. Our simulated genome (GENERIC
# species, 1 Morgan : 100Mb) with a few hundred markers per chromosome has average inter-marker
# spacing of several Mb -- with the default maxGap, nearly every consecutive marker pair exceeds
# it, so the run-building algorithm fragments at almost every site and finds nothing regardless of
# actual homozygosity (confirmed empirically: even a 100%-homozygous synthetic genotype vector
# produced zero runs at the default maxGap -- see dev/test_pedconvert.R). maxGap is scaled here to
# comfortably exceed the map's actual max inter-marker gap instead.
compute_froh <- function(pop, SP, chip, minSNP = 15, minLengthBps = 100000, ROHet = FALSE) {
  tmp_dir <- tempfile("plink_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))
  ped_path <- file.path(tmp_dir, "geno.ped")
  map_path <- file.path(tmp_dir, "geno.map")
  write_plink_files(pop, SP, chip, ped_path, map_path)

  map <- getSnpMap(snpChip = chip, simParam = SP)
  bp <- round(map$pos * 1e8)
  max_observed_gap <- max(tapply(bp, map$chr, function(x) max(diff(sort(x)), 0)))
  maxGap <- max(10^6, max_observed_gap * 2)

  runs <- consecutiveRUNS.run(
    genotypeFile = ped_path, mapFile = map_path,
    ROHet = ROHet, minSNP = minSNP, minLengthBps = minLengthBps, maxGap = maxGap
  )
  if (nrow(runs) == 0) {
    ids <- rownames(pullSnpGeno(pop, snpChip = chip, simParam = SP))
    return(data.frame(id = ids, froh = 0))
  }
  froh_df <- Froh_inbreeding(runs = runs, mapFile = map_path, genome_wide = TRUE)
  data.frame(id = froh_df$id, froh = froh_df$Froh_genome)
}
