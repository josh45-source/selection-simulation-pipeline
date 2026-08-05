Load-bearing (evidence behind real design decisions, keep): `verify_mutate_heterozygosity.R`,
`verify_va_stability.R`, `verify_mutation_revert.R`. MaCS root-cause investigation (one-off,
superseded by the README/sim_loop.R writeup, kept only as evidence): `verify_founder_reproducibility.R`,
`verify_macs_direct.R`, `verify_rng_subprocess.R`, `verify_same_process_reproducibility.R`,
`verify_seed_derivation.R`, `verify_serialization_roundtrip.R`. Utilities: `bootstrap_renv.R`
(regenerates `renv.lock`), `Dockerfile.dev` (fast-iteration image, not the production build).
