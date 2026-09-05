# Recurrent mutation via AlphaSimR's mutate() [ref: gaynor2021alphasimr], which applies true
# per-haplotype-copy mutation (heterozygous-capable), NOT forced homozygosity like editGenome().
# editGenome() was rejected for this purpose: it "edits selected loci ... to a homozygous state,"
# which would fabricate homozygosity at every mutated site -- directly contaminating F_ROH and the
# heterozygosity-decay validation. Confirmed empirically (dev/verify_mutate_heterozygosity.R):
# mutate() produces genuine homozygous -> heterozygous transitions.
#
# Recurrent mutation is required over long horizons because a fixed founder haplotype set
# exhausts standing variation and selection response flatlines as a model artifact.
# ref: robertson1960limits, hill1986recurrent
#
# mutate(pop, mutRate, simParam) applies ONE rate genome-wide per call -- there is no argument to
# restrict which chromosome/site subset is affected. To keep validation_markers at mu = 0 (required
# so the neutral heterozygosity-decay test isn't confounded by mutation-drift equilibrium), we
# snapshot that chip's haplotypes before the call and restore them after.
#
# The restore does exactly what the validation panel needs: chip 2 comes back bit-identical, so it
# is genuinely held at mu = 0. It does NOT, however, leave the rest of the genome untouched, which
# an earlier version of this comment claimed. AlphaSimR's two SNP chips are drawn independently
# from the same pool of non-QTL sites and overlap heavily (see R/engine.R), so restoring chip 2
# also rolls back any mutation that landed on a site chip 2 shares with chip 1. Measured: about 8%
# of the mutations mutate() places on chip 1 are silently reverted this way (3 of 37 in one run at
# the slice config's rate). The QTL are unaffected -- they overlap neither chip -- so the trait,
# and every metric derived from it, is untouched.
#
# This is accepted rather than fixed because chip 1 is not used for any reported metric; see the
# note in R/engine.R for why fixing the allocation would cost more than it buys.

library(AlphaSimR)

apply_recurrent_mutation <- function(pop, SP, mutRate, validation_chip) {
  if (mutRate <= 0) return(pop)
  validation_before <- pullSnpHaplo(pop, snpChip = validation_chip, haplo = "all", simParam = SP)
  # Namespace-qualified: dplyr (loaded by R/io.R for write_dataset()) also exports mutate(), and
  # masks AlphaSimR::mutate() on the search path depending on library() load order.
  pop <- AlphaSimR::mutate(pop, mutRate = mutRate, simParam = SP)
  pop <- setMarkerHaplo(pop, haplo = validation_before, simParam = SP)
  pop
}
