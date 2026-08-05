# Selection simulation MVP

Status: early scaffold. Sim -> Parquet -> DuckDB -> plot slice is working end-to-end; mutation
model is empirically verified (see `dev/`); checkpoint/resume, F_ROH, validation suite, and the
Quarto dashboard are still to come. Full setup/build/run instructions land here once the dashboard
step is done -- this is a stub covering just the one thing needed now: how to run a container
against this repo without host-side permission surgery.

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
