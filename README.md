# Selection simulation MVP

Status: early scaffold. Sim -> Parquet -> DuckDB -> plot slice is working end-to-end; mutation
model is empirically verified (see `dev/`); checkpoint/resume + RNG streams are working and
verified bit-for-bit (`validation/test_resume_equals_uninterrupted.R`); F_ROH, the rest of the
validation suite, and the Quarto dashboard are still to come. Full setup/build/run instructions
land here once the dashboard step is done -- this is a stub covering just what's needed now: how
to run a container without host-side permission surgery, and a known reproducibility limitation.

## Known limitation: founder generation is not reproducible from a seed

`AlphaSimR:::MaCS()` (the compiled coalescent function `runMacs()` calls to build founder
haplotypes) does not produce identical output across calls given identical literal seeds --
confirmed empirically, bypassing R's own RNG entirely (see `dev/verify_macs_direct.R`). This is a
limitation in AlphaSimR's compiled code, not something fixable from this project, and it's
undocumented -- AlphaSimR's changelog shows the maintainers forced single-threading in their own
doc examples (v1.5.1) but never describes fixing MaCS's seed handling itself.

Practical effect: a scenario config + seed is enough to reproduce everything *after* founders
exist (mutation, selection, crossing, checkpoint/resume -- all proven bit-identical), but not the
founder population itself. `run.R` therefore persists each replicate's founder population once, to
`data/founders/<scenario>_rep<n>.rds`, immediately after it's built (reusing the checkpoint
mechanism at generation 0 -- see `R/sim_loop.R`'s `founders_path` argument). Reproducing a specific
run's founders requires that file, not just its seed. A full Provenance dataset (git SHA, container
digest, full parameter set) is still to be built; whatever it references for "how to reproduce this
run" must point at this file.

## Running a scenario

The image runs as a non-root user, so the bind-mounted `data/` output directory must be owned by
whatever UID actually runs the container -- otherwise Docker auto-creates it as root on first
mount and the container can't write to it. `--user "$(id -u):$(id -g)"` runs the container as your
own host UID/GID instead of the image's baked-in default, so this works regardless of what UID you
are on your machine (not just uid 1000):

```bash
docker build -t selection:latest .

mkdir -p data
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd)/data:/project/data" \
  selection:latest run.R config/scenario.slice.yaml   # fast smoke test
  # or: config/scenario.example.yaml for the full-scale scenario
```
