Load-bearing (evidence behind real design decisions, keep): `verify_mutate_heterozygosity.R`,
`verify_va_stability.R`, `verify_mutation_revert.R`, `verify_pedigree_inbreeding_known_answer.R`
(proves the tabular-method F calculator against a hand-calculable full-sib-mating case).

MaCS reproducibility investigation (one-off root-cause analysis, superseded by the README/
sim_loop.R writeup, kept only as evidence): `verify_founder_reproducibility.R`,
`verify_macs_direct.R`, `verify_rng_subprocess.R`, `verify_same_process_reproducibility.R`,
`verify_seed_derivation.R`, `verify_serialization_roundtrip.R`.

detectRUNS/pedigree investigation (one-off, superseded by R/froh.R and R/pedigree.R's own
comments, kept only as evidence): `verify_plink_maxgap_diagnosis.R` (found the maxGap-too-small
bug), `verify_detectruns_internals.R` (isolated it to consecutiveRunsCpp via detectRUNS' own
bundled example data), `verify_monoecious_role_conflict.R` (found kinship2 unusable for AlphaSimR
pedigrees), `verify_pedigree_slot_format.R` (AlphaSimR id/mother/father format).

Utilities: `bootstrap_renv.R` (regenerates `renv.lock`), `Dockerfile.dev` (fast-iteration image,
not the production build).
