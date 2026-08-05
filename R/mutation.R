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
# snapshot that chip's haplotypes before the call and restore them after. Confirmed empirically
# (dev/verify_mutation_revert.R) that this reverts only the targeted chip and leaves mutations
# elsewhere in the genome intact.

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
