> **This is the public code of PlanetFinder.** The day-to-day analyses, candidate
> reviews and status notes referred to below (`investigations/`, `records/`, `STATUS.md`,
> `targets/`) live in a private working repo. Candidates are published, with their
> evidence, in [PlanetCandidates](https://github.com/moosesteve21/PlanetCandidates) once
> they have been submitted.

# PlanetFinder

An unattended blind transit search of TESS light curves, aimed at small, lightly
observed dwarf stars, with automatic vetting to reject false positives. It also
flags interesting non-planets, such as eclipsing binaries and spotted stars, for
later investigation.

## Project overview

**Who and why.** A one-person citizen-science project (Steve Price; Zooniverse:
moosehunter21), inspired by Pavel Rabtsevich's TESS candidate around TIC 4206066, found with
Claude Code (preprint doi:10.5281/zenodo.22967456; that star is now the regression benchmark). It runs on a home Ubuntu server (8 cores / 16 threads, 30 GB) with an
AI coding assistant (Claude Code). The aim is to find transiting planets that the
official TESS pipelines and the TOI process have missed, around bright (T 8–11), small
(K and M) dwarfs that TESS has observed in only 1–6 sectors. Variable stars found along
the way (eclipsing binaries, spotted stars, wrong catalogue periods) are a second output.

**Where results go.**
- **Planet candidates** are posted to the Planet Hunters TESS Talk "Planets!" board for
  community and science-team review. ExoFOP now accepts community TOIs (CTOIs) only with
  a peer-reviewed publication, so Talk is the route.
- **Variable stars** go to the AAVSO VSX catalogue, in one session a week, from an
  automatically built shortlist.
- **The public repo** [PlanetCandidates](https://github.com/moosesteve21/PlanetCandidates)
  is the permanent record of what was submitted: plots, discovery records, Talk posts,
  VSX form data and reports (see "Submitting" below).
- This repo holds the code and the analyses (`investigations/`). The database lives on
  the server, with its irreplaceable tables (notes, watchlist) copied to `records/`.

**How it's run.**
- The sky is worked through in numbered **batches** (b001, b002, ...) of the
  highest-scoring stars not yet searched: 300 stars for daytime runs (about 1¼ hours),
  500 overnight.
- **The owner starts every batch** and reviews its results before starting the next.
  Batches are never chained or started automatically.
- **Preparing the next batch is a manual step,** never done by the pipeline. The owner
  (or the assistant, when asked) selects it with `tess_target_finder.py` and prefetches
  it, usually while the current batch is still running, so it's ready to start.
- When a batch ends, a Discord message reports the outcome (`server/notify.sh`).
- Commits are made on the server; the owner pushes to GitHub (the server holds no
  GitHub credentials). GitHub is also the off-site backup.

**After each batch.** The pipeline ends by loading the database, backing it up and
posting the outcome to Discord. Then the owner and assistant:
1. Check the Discord message and `server/batch_status.sh` (complete, stopped or failed).
2. Review any clean candidates: an independent pipeline, per-sector depths, pixel
   position, density and FPP. Then decide: post, watchlist, or reject. Record the
   verdict as a note with `tessdb.py note`. For a candidate worth keeping (posted later,
   or on the watchlist), run `investigations/fingerprint.py` to log a dated fingerprint
   in the public PlanetCandidates `discovery_log.md`.
3. Look over the batch's VSX shortlist, saving submissions for the weekly session.
4. Commit `records/` and update `STATUS.md`.
5. Prepare the next batch, if one isn't ready, and the owner starts it.

**Submitting.** A Talk post or a VSX submission is prepared in full (README, plots,
post or form data, report) in the PlanetCandidates repo. That folder is committed as
the last step before the owner posts or submits, so the repo records exactly what went
out. Every public submission (Talk post, VSX report, email) links both the candidate's
PlanetCandidates folder and the public code repo, https://github.com/moosesteve21/PlanetFinder.

**What makes a candidate worth posting.** Surviving the automatic checks ("clean") is
only the start. Before posting, a candidate should also have:
- the dip in every sector, at a consistent depth;
- confirmation by a pipeline the search didn't use (such as TARS);
- a transit shape implying the catalogued star's density;
- odd and even transits that match, with no secondary eclipse;
- evidence that the dip is on the target (centroids, or a per-pixel depth map when there
  is no 2-min data);
- a TRICERATOPS false-positive probability that isn't dominated by nearby-star scenarios.

Weaker or doubtful signals go on the **watchlist**, to be re-tested as new sectors arrive.

**Roadmap.** Finish the current target list (stars with 1–6 sectors). Then refresh the
catalogue for new sectors, and extend the search to stars with 7–10 sectors, then 11–30
(after capping the per-sector checks), then TESS's polar continuous-viewing zones.

**Where to look next.** `STATUS.md` (current state and results), `EXAMPLE_RUN.txt` (one
full batch, step by step), `investigations/README.md` (every close look at a star), and
`server/README.md` (boot-time sessions, notifications).

## The pipeline at a glance

1. **Choose targets:** `tess_target_finder.py` builds a whole-sky catalogue from the TESS
   Input Catalog and writes each batch: the highest-scoring unsearched stars.
   - **Prefetch:** `prefetch.py` downloads a batch's light curves ahead of its search.
2. **Search:** `run_parallel.sh` runs `find_planets.py` (BLS search and automatic
   vetting) on a shared work queue, with several workers.
3. **Catalogue checks:** `enrich_candidates.py` adds ExoFOP, SIMBAD, Gaia and
   known-variable data.
4. **Pixel checks:** `pixel_check.py` tests whether each surviving dip comes from the target.
5. **Reduction check:** `reduction_check.py` tests whether the dip appears in independent
   reductions of the same data.
   - **VSX shortlist:** `vsx_shortlist.py` lists the variable-star signals for VSX.
6. **Transit fit:** `transit_fit.py` (batman + emcee) checks clean candidates' implied
   stellar density and duration.
   - **Depth map:** `depth_map_check.py` finds which star the dip comes from, star by star.
7. **False-positive probability:** `fpp_check.py` runs TRICERATOPS on consistent fits
   that the depth map doesn't put on another star.

`tessdb.py` keeps every run in one SQLite database, browsable with Datasette. After
choosing (and ideally prefetching) a batch, one command runs steps 2–7 (with the VSX shortlist), loads the
results into the database and backs it up, unattended:

```bash
./pipeline.sh targets/b0xx.txt results/b0xx 8 "batch b0xx"
```

**Two stages.** `pipeline.sh` runs the search (CPU-bound). When it's complete, Discord
gets "b0xx search done … next batch can start", and the rest (steps 3–7, ingest, backup;
mostly waiting on MAST) moves to `followup.sh` in a `tess` tmux window named `fu-b0xx`,
logging to the same `pipeline.log` and posting "b0xx complete" at the end. So the next
batch's search can start as soon as the previous search is done. Follow-ups queue and run
one at a time, in order, at low CPU priority (`nice`), so a search running alongside keeps
its speed. `./followup.sh results/b0xx` redoes a batch's follow-up by hand;
`FOLLOWUP=inline ./pipeline.sh …` runs everything in one go, as before. (In b025 the
follow-up took about 2 h, waiting on a slow MAST: the next search no longer waits for it.)

`watchlist.py` re-tests chosen candidates as new sectors are released.

`pipeline.sh` logs each step to `results/<name>/pipeline.log` and stops at the first
failure. The follow-up steps run only when the search is complete (see step 2), and it
ends by backing up the database. Rerunning the same command resumes. To pause a long
search, `touch results/<name>/STOP`: the pipeline exits with status 6 once the workers
have finished their current stars, and rerunning continues.

**Prefetching (recommended for long batches):** fetch a batch's MAST data ahead of
its search, whenever MAST is healthy (for example right after creating the batch):

```bash
.venv/bin/python prefetch.py targets/<name>.txt results/<name>     # rerun until complete
```

For each star it runs the search's own MAST steps (product search and choice,
downloads, TIC radius and mass) and records them in `results/<name>/prefetch.jsonl`.
When the batch runs, `run_parallel.py` points the workers at the manifest and they read
prefetched stars entirely from the cache, so a MAST outage during the run doesn't stall
it. Stars not in the manifest, or whose files have gone, are searched the usual way. A search started while the prefetch is still running waits for it to finish (`.prefetch.lock`; the log says "waiting for the prefetch"), and a second prefetch on the same folder is refused. Tested
with all network access blocked: a prefetched star gave the same result as its original
search. A manifest more than two weeks old may miss new sectors; use `--refresh`. The
follow-up checks still need the network, but only for the few signals.

**Rolling updates:** during the search, every 100 stars finished (`ROLLING_INGEST`, set
to 0 to turn off) the results so far are loaded into the database, so a long batch can be
browsed in Datasette as it runs. The run's note shows `[in progress: x/y searched]`.
A run paused with `STOP` (or interrupted) is loaded on the way out too, marked
`[stopped: x/y searched]`. These loads hold the search's signals, flags and plots; the catalogue, pixel, reduction,
fit and FPP checks are loaded at the end, with the final load replacing the interim ones.
Plot links work for a run started after `tessdb.py serve`, because the whole `results/`
folder is mounted.

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Developed with Python 3.14, astropy 8 and numpy 2.5 (pinned in `requirements.txt`).
`run_parallel.sh` uses `.venv/bin/python` automatically.

**Vetting environment** (for steps 6 and 7): TRICERATOPS needs an old PyTransit (2.2)
and numpy (1.23), so the transit fit and false-positive check run in a separate
Python 3.11 environment, called by the main scripts as a subprocess:

```bash
uv venv --python 3.11 .venv-vet
uv pip install --python .venv-vet/bin/python -r requirements-vet.txt
```

Without `.venv-vet`, `pipeline.sh` skips those two steps. The current PyTransit (2.9)
can't replace it yet: it needs celerite, which has no Python 3.14 build, and the
server has no compiler. When TRICERATOPS drops its PyTransit 2.2 pin, a single Python
3.12 vetting environment could hold both.

## Workflow

### 1. Choose targets

**Recommended: batches from the whole-sky catalogue.** The best targets are a property of
each star rather than of a patch of sky, so the finder can rank every suitable star
across the sky and hand out the best ones not yet searched:

```bash
.venv/bin/python tess_target_finder.py --build-catalog          # once (about 30 min; resumable)
.venv/bin/python tess_target_finder.py --batch 90 --name b001   # each new batch
tmux attach -t tess                                             # then, at its prompt:
./pipeline.sh targets/b001.txt results/b001 8 "batch b001"
```

Run batches at the prompt of the `tess` tmux session (created at boot; if it's missing,
`server/start_sessions.sh` recreates just that one), not as the session's own command:
`tmux new -s tess './pipeline.sh …'` closes the session when the run ends, taking its
final output with it.

- **`--build-catalog`** builds the catalogue in `db/tess.db`:
  - it queries the TIC in 15° sky boxes, 4 at a time, applying the stellar cuts below,
    and records finished boxes so an interrupted build resumes;
  - every star gets its sector count (tess-point), a known-system flag (TOI, confirmed
    planet, TESS eclipsing binary), its Gaia RUWE, and a **score**.
- **The score** is the expected SNR of a 2 R⊕ planet on a 10-day orbit. It comes from the
  star's T magnitude (an approximate TESS noise model), radius, mass and number of
  sectors. For comparison, injection tests give 50% completeness at an expected SNR
  of about 12.
- **`--batch N --name NAME`** picks the N highest-scoring eligible stars:
  - eligible means 1–6 sectors, not a known system, and RUWE ≤ 1.4;
  - it skips stars already searched in any run in `db/tess.db`, and stars handed out in an
    earlier batch;
  - it checks the picked stars for SPOC 2-minute data, writes
    `targets/NAME.csv` / `.txt` (with criteria header, best first) and records the
    batch in the `batches` table.
- **Preparing ahead:** batches can be made in advance, since each one's stars are held
  for it. The known-system list is fixed when the catalogue is built, so a batch made
  early picks the same stars as one made later; light curves are downloaded when the
  search runs, so new sectors are still used.
- **`--release NAME`** returns an unused batch's stars to the pool and deletes its
  target files, for example after changing the selection criteria. It refuses once the
  batch has been run or started (`results/NAME` holds search output, or any of its stars
  has been searched). A batch that has only been prefetched can be released; its
  `prefetch.jsonl` is deleted too.
- **Refreshing:** rebuilding the catalogue (for example after new sectors) updates
  magnitudes and sectors and keeps batch assignments. Delete `catalog_tiles` rows to
  re-query boxes.

**Alternative: one patch of sky (cone).** Running the script with no options does the
original cone search, around `FIELD_RA`, `FIELD_DEC` within `FIELD_RADIUS`.

The cone search queries the TESS Input Catalog for dwarf stars in one patch of sky.
It keeps stars that are lightly observed (1–6 sectors) and have little light from
neighbouring stars. It removes stars that already have a TOI or a confirmed planet,
and known TESS eclipsing binaries. It also removes likely binaries using Gaia's
RUWE measure, which is high when a star wobbles on the sky.

For a new field, edit the settings at the top of the script:

- `FIELD_NAME`, `FIELD_RA`, `FIELD_DEC` and `FIELD_RADIUS`;
- the stellar cuts, such as `TMAG_RANGE`;
- the output names, `OUTPUT_CSV` and `OUTPUT_TXT`.

The script writes two files:

- **The CSV** (for example `targets/aries.csv`): full details, ranked.
- **The target list** (for example `targets/aries.txt`), read by the search:
  - its header records the selection criteria, the run date and how many stars
    survived each step;
  - stars are ordered with the least contaminated first;
  - a `2min` comment marks stars with SPOC 2-minute data, which the official pipeline
    has already searched.

Change the output names for each field, or the previous field's files are overwritten.

**Choosing a field:**

- **Near the ecliptic** gives few sectors and fresh data.
- **Far from the Galactic plane** (|b| ≳ 30°) gives little crowding. The first field, on
  the plane, was dominated by blends.
- **Outside K2 campaign footprints**, whose bright dwarfs K2 already searched.
- **T ≤ 11** is needed to reach planets smaller than about 3 R⊕ (see Sensitivity).
  At that limit a radius of 8–9° gives about 70–100 targets.

| Field | Centre | Target list |
|---|---|---|
| First FFI field (Galactic plane) | RA 95°, Dec +15°, r 4°, T 8–13 | `targets/ffi.txt` (70) |
| Aries | RA 40°, Dec +15°, r 8.5°, T 8–11 | `targets/aries.txt` (91) |

### 2. Search

```bash
./run_parallel.sh targets/<name>.txt results/<name> 8      # at the tess session's prompt
```

Arguments are: the target list, an output folder, the number of workers
(default 8), and optionally the script to run (default `find_planets.py`).
Run it inside `tmux` so it survives a closed terminal or editor. Detach with
`Ctrl-b d` and reattach with `tmux attach -t tess`.

- **Shared queue:** the targets go into `results/<name>/queue.db` (SQLite) in file
  order, which is score order for batches. Each worker claims the best-ranked star
  that is due, searches it, records the outcome and claims the next, so no worker
  sits idle while stars are left. Worker *k* writes only to `results/<name>/w<k>/`
  and uses one CPU thread. `run_parallel.sh` calls `run_parallel.py`, which
  supervises the workers.
- **Worker count:** about 8 works well. The MAST archive limits how fast it serves
  downloads, so more workers mostly cause more download errors. The count can differ
  between runs of the same folder, for example fewer workers while you use the server.
- **Automatic retries:** a star that ends in an error goes back into the queue and is
  retried after `RETRY_WAIT` seconds (default 300), up to `MAX_ATTEMPTS` attempts
  (default 3). A worker that dies has its star returned to the queue (that counts as
  an attempt) and is replaced.
- **Pausing:** `touch results/<name>/STOP`. Workers finish their current star and exit,
  and the run stops with status 6. Rerunning the same command removes `STOP` and
  continues. Ctrl-C or `kill` on the runner does the same, and returns the stars being
  searched to the queue without charging them an attempt.
- **Merged output:** at the end, `candidates_all.csv` (each star's rows from the
  attempt that finished it), `done_all.txt` (one line per star: star, status, time,
  worker), `failed_targets.txt` and `attention.txt` are written from the queue.
- **Exit status:**
  - `0`: complete.
  - `3`: incomplete: stars used all their attempts; they are listed in
    `failed_targets.txt`, and rerunning the same command gives them fresh attempts.
  - `4`: another run is using the folder; a lock prevents clashes.
  - `6`: stopped (`STOP` file or interrupted) with stars left; rerun to continue.

  Stars whose search stopped partway (`partial_…`) count as complete, but are listed in
  `attention.txt` because their candidate lists may be short.
- **Rerunning** the same command resumes: finished stars are skipped, claims left by a
  dead worker are returned to the queue, and stars added to the target file are queued
  at the end. The original commit and settings stay in `run_info.json`, and each resume
  is added under `resumes`. Folders from the older round-robin runner (no `queue.db`)
  are resumed by seeding the queue from their `w*/done.txt` files.
- **Separate runs:** use a new output folder for each batch of targets.
- **Settings:** `STAGGER` sets the seconds between worker starts (default 20).
- **Testing changes:** `tests/test_run_parallel.sh` checks the queue, retry, stop and
  resume logic in a few seconds, using a fake search (`tests/fake_search.py`).

`find_planets.py` run on a list by itself exits with 0 when every target has a final
status and 3 if any ended in an error. A failed download counts as an error, so it's
retried. Under `run_parallel` it runs in worker mode (`--queue`).

Light curves are cached in `~/.cache/lightkurve`, so reruns are much faster.
Expect a few minutes per star per worker on the first run.

**MAST timeouts:** the per-star scripts give up on a MAST request after 90 s of
silence (`mastutil.py`) instead of astroquery's default 600 s, and their own retries,
and the run's queue, try again. In b008, two waves of MAST trouble cost 24 stars about
600 s per attempt. `tess_target_finder.py` keeps the default, because its whole-sky
catalogue queries can legitimately take minutes.

**Disk space:** about 90% of stars log no signal. Their downloaded light curves
(about 6 MB a star) are deleted from the cache once their search ends; stars with any
signal keep theirs for the follow-up checks and later review. Searching a deleted star
again just downloads it again. `--keep-lightcurves` (on `find_planets.py`) keeps
everything. `tessdb.py clean-cache` does the same tidy-up for stars searched earlier.

### 3. Catalogue checks

```bash
.venv/bin/python enrich_candidates.py results/<name>/candidates_all.csv
```

For each star with a candidate, the script queries:

- **ExoFOP-TESS:** TOIs and community TOIs, follow-up observations and known companions.
- **SIMBAD:** identifier, object type and spectral type.
- **Known variables within 3′:** VSX, Gaia DR3 eclipsing binaries and classifications,
  and ZTF periodic variables.
- **TIC neighbours within 5′.**
- **Gaia DR3 binarity measures for the target:** RUWE, a binary-orbit solution,
  radial-velocity variability, BP/RP flux excess and double-peaked images.

Every candidate is then compared with these known variables. A nearby variable whose
period matches P, 2P, P/2, or a 3:2-family alias (3P/2, 2P/3, 3P, P/3) is almost
certainly the real source of the signal. Ground-based surveys observe about once per night,
so their catalogued period can be a one-day alias of the true one (frequencies differing by
1.0027 cycles per day). Those count as matches too, labelled `1-day alias of P`.

The script also lists neighbouring stars bright enough to produce the observed depth
if they had a deep eclipse (`blend_suspects`). Queries are cached in `enrich_cache/`.
The output is `enriched.csv`. Rerunning later picks up new ExoFOP data, but delete the
star's file in `enrich_cache/` first, or the old cached answer is reused.

**VSX shortlist.** `vsx_shortlist.py results/<name>` (a pipeline step after the pixel and
reduction checks) lists the run's signals that look like variable stars rather than transits, for
the weekly VSX session. The catalogue step (`enrich_candidates.py`) checks VSX, Gaia DR3 (eclipsing binaries and
classes), ZTF and the TESS Ten Thousand Catalog. It looks for sinusoidal variability, a secondary eclipse, a
phase-curve residual at P ≤ 3 d, or an odd/even difference over enough transits. TOIs,
CTOIs and signals matching a variable on another star are skipped. Each listed star is
in one of four groups:
- `new`: no VSX entry at the target. Other catalogues that have it, such as Gaia DR3
  classes or ZTF, are listed alongside.
- `published`: not in VSX (or in VSX without this period), but the **TESS Ten Thousand
  Catalog** (Kostov et al. 2025, 2025ApJS..279...50K; VizieR J/ApJS/279/50, tables 3 and 4)
  has it as an eclipsing binary at this period. It still belongs in VSX, but as a published
  star with that reference, not as our discovery. Most of its 7,936 new EBs aren't in VSX yet.
- `period`: in VSX, but with no period, or a period TESS doesn't confirm. A match only as
  a one-day alias counts as unconfirmed, because TESS observes continuously, so it's the
  catalogue period that's the alias.
- `blend`: the pixel check puts the dip off target, and a star there is bright enough
  to cause it (an eclipse under 50% deep). It's probably that variable neighbour, named in
  `blend_suspects`, which may itself be new to VSX. If no star there is bright enough,
  the signal stays in its own group with a note.

Signals the reduction check calls `reduction_dependent` are left out, because they
are most likely systematics rather than a variable.

Stars already noted as `vsx`, `vsx-candidate` or `submitted` show that note. The
output is `vsx_shortlist.csv`, loaded into the database as `vsx_shortlist_view`. Run
it by hand on older runs; it reads only local files and loads them into the database.

### 4. Pixel checks

```bash
.venv/bin/python pixel_check.py results/<name>/enriched.csv
.venv/bin/python pixel_check.py results/<name>/enriched.csv --target "TIC 123456789"
```

The script downloads full-frame-image cutouts for each sector and runs these tests:

- **Aperture test:** compares the dip depth in the target's own pixel with a 3×3 box.
  If the depth grows in the larger box, the signal comes from a neighbour.
- **Difference image:** subtracts the average in-dip image from the average out-of-dip
  image. The centre of the dimming is converted to a sky position and combined across
  sectors.
- **Source star (off-target dips only):** downloads the source star's own light curve
  and reports two things:
  - its strongest periodicity, with amplitude, and how that period relates to the
    candidate's;
  - how much it dims at the candidate's dip times. If the source's own variation is
    unrelated to the candidate's period, it is removed first. The significance is
    judged against the same measurement at 60 random times.

  A deep, significant dip means the source is an eclipsing binary at that period.
  No dip means the candidate lined up with some other variation of the source by
  chance, for example the minima of a spotted star.

The verdict is either *consistent with the target* or *OFF-TARGET*. An off-target
verdict lists the nearest catalogued stars and any known variables at the dimming
position, each with its distance from that position. The search stays well short of
the target, so the target itself is never named. With fewer than 3 usable sectors
the verdict is marked tentative.

How to read the plots, one per candidate and sector:

- **Left panel (median image):** each pixel's typical brightness. Colours are scaled
  to the brightest pixel, so a brighter neighbour can outshine the target.
- **Right panel (difference image):** how much each pixel dims during the dips.
  A real dip shows as a compact, star-shaped bright patch on its source. Random
  speckle, together with a very small colour scale, means the dip isn't detected
  in that sector.
- **Markers:**

| Marker | Meaning |
|---|---|
| Red × | The target |
| Magenta + | Centre of the dimming; only drawn if some pixel is more than 3σ above the noise |
| White circle with a number | A catalogued star up to 1 magnitude fainter than the target, or brighter; the number is its TESS magnitude |
| Orange diamond | A known variable (VSX, Gaia DR3 eclipsing binary or ZTF), labelled with its type and, if known, its period |

Variable types use the VSX codes, for example: `EA`/`EB`/`EW` eclipsing binaries
(detached / semi-detached / contact), `RS` RS CVn-type spotted active binary,
`BY` spotted dwarf, `ROT` rotational variable, `VAR` variable of unspecified type.

By default the script checks only candidates still worth testing:

- `clean` rows;
- `double_period_variable` rows with enough transits;
- and in both cases, only rows without a catalogue period match.

Use `--target` (and optionally `--candidate`) to check any star.

Results go to `pixel_check/`:

- `pixel_check.csv`: one row per candidate and sector.
- `pixel_summary.csv`: one verdict per candidate, with the source-star analysis.
- A plot for each candidate and sector.

**Each run overwrites both CSV files**, so use `--out` to keep runs separate. The cutout service often asks clients to slow down, so new sectors can
take 1–3 minutes each. Downloads are retried automatically.

## Sensitivity (injection–recovery)

```bash
.venv/bin/python inject_recover.py targets/ffi.txt --n 500 --workers 8 --out results/injection
.venv/bin/python inject_recover.py --summary results/injection
.venv/bin/python inject_recover.py targets/aries.txt --rp 0.7 6 --out results/injection_aries
.venv/bin/python inject_recover.py targets/aries.txt --invert --out results/falsealarm_aries
```

- **Ranges:** `--rp MIN MAX` and `--period MIN MAX` override the default radius and
  period ranges.
- **False-alarm test:** `--invert` searches each star's light curves flipped upside down.
  Real transits become brightenings the search ignores, while noise and spot variability
  still make dips. Every signal it finds is therefore a false alarm, which gives the field's
  false-alarm rate. Results go to `falsealarms.csv`.

Each job injects one synthetic transit into a field star's real light curves:

- **Planet:** radius 1–8 R⊕ and period 0.5–30 d (both drawn evenly on a log scale), random
  impact parameter and epoch.
- **Shape:** a trapezoid sized from the star's TIC radius and mass.

The unchanged search then runs on the result. The summary gives the recovered fraction
by expected SNR, and by planet radius and period. Jobs come from a fixed random seed and
results are appended as they finish, so rerunning the same command resumes a stopped
run. Expect 2–5 minutes per job per worker.

**First field (70 K dwarfs, T ≈ 11–13, 2–6 sectors each), 600 injections, 11 October 2026:**

- **By expected SNR:** 50% of transits are recovered at an SNR of 12.4 and 90% at 18.5.
  None below SNR 7 are recovered.
- **By planet size:**

| Radius (R⊕) | P 0.5–1 d | P 1–3 d | P 3–10 d | P 10–30 d |
|---|---|---|---|---|
| 1–1.5 | 4% | 9% | 0% | 4% |
| 1.5–2 | 20% | 13% | 6% | 5% |
| 2–3 | 29% | 30% | 21% | 11% |
| 3–4 | 71% | 72% | 46% | 38% |
| 4–6 | 82% | 93% | 83% | 59% |
| 6–8 | 100% | 93% | 94% | 84% |

So for stars like these, "no candidates" rules out most planets larger than about 4 R⊕
inside 10 days, but says little about anything below about 3 R⊕.

**Aries field (91 K dwarfs, T ≤ 11, mostly 4–5 sectors), 300 injections of 0.7–6 R⊕:** a
similar SNR threshold (50% at an expected SNR of about 10–15), but the brighter stars reach
much smaller planets:

| Radius (R⊕) | P 1–3 d | P 3–10 d | P 10–30 d |
|---|---|---|---|
| 1–1.5 | 0% | 7% | 0% |
| 1.5–2 | 56% | 50% | 44% |
| 2–3 | 85% | 92% | 83% |
| 3–4 | 100% | 95% | 60% |

**False alarms (Aries, `--invert`):** 91 flipped stars gave 9 false signals, **none `clean`**,
so clean false alarms occur at fewer than about 3% per star (95% limit).

- **Flagging:** 14% of recovered injections were flagged. That was mostly `few_transits`
  at long periods (fewer than 5 transits), occasionally `odd_even_mismatch` (about 3%).
  This test led to `long_duration` being used only as a fallback, after it flagged genuine
  ultra-short-period transits.

### 5. Reduction check (survivors only)

```bash
.venv/bin/python reduction_check.py results/<name>/enriched.csv
```

A real transit is in the star's light, so it appears in every reduction of the same
observations. A processing artefact, for example from a systematics correction near a
momentum dump, appears in only one. For each surviving candidate (selected as in step 4)
the script:

- **loads every reduction:** SPOC 2-minute and TESS-SPOC (PDCSAP and SAP flux), plus
  QLP's raw and detrended flux;
- **measures each predicted transit** in each one;
- **compares stacked depths** with the reduction the search used, on the events both cover.

A reduction **supports** the signal at > 3σ and at least half the depth, or **contradicts**
it when it's significantly shallower than that. The signal is `reduction_dependent` when
independent reductions contradict it and none supports it.

**Testing:** all 10 control planets and the 500 ppm benchmark were `consistent`, with 5 of 5
reductions supporting. An early candidate from the Aries field was `reduction_dependent`, with 5 of 5
contradicting: its strongest dips appear only in SPOC 2-minute PDCSAP flux. Results go to
`reduction_check/`, and `tessdb.py` loads them.

### 6. Transit fit and density check (clean candidates)

```bash
.venv/bin/python transit_fit.py results/<name>/enriched.csv
```

A transit's shape and duration fix a/R*, and with the period the star's density
(ρ* = 3π (a/R*)³ / G P²). If the dip is on the catalogued star, that density matches
the TIC's. A blend with a different star, or an eclipsing binary, usually doesn't.
For each clean candidate, `transit_fit.py`:

1. normalises up to 60 transits by a local quadratic;
2. looks up the star's TIC mass, radius and Teff;
3. runs a batman + emcee fit (`vetting/transit_mcmc.py`, in `.venv-vet`) with a
   circular orbit, quadratic limb darkening with Teff-based priors, and free epoch,
   period, Rp/R*, b, log ρ* and noise scale.

It compares the fitted and catalogued densities in log10, with both uncertainties:
within 3σ is `consistent`; beyond is `density_mismatch` (also added to the
candidate's flags). A fitted T14 over 1.3× the longest transit the star allows
(central, circular, `max_t14_h`) is `duration_too_long` instead, whatever the density
(flag `fit_duration_too_long`): such fits can pass the density test through wide
errors (a b023 blend: 7.7 h against 4.4 h, density 1.9σ; a b017 fit that wandered to
9.9 h). Real candidates so far reach 1.16×. Only `consistent` fits go on to the FPP.
Results go to `transit_fit/transit_fit.csv`: Rp, b, a/R*, T14, ρ* and z.

- **Timing variations:** planets with transit-timing variations (e.g. near a
  resonance) smear the stacked dip, which the fit reads as a grazing transit of a
  low-density star. Each sector's timing offset is measured (`ttv_p`,
  `ttv_max_min`), and a `density_mismatch` whose timing scatters significantly
  carries a note saying so. Aligning sectors before fitting was tried and
  rejected: it fixed TOI-270 c and d, but for weak signals it aligned noise
  (TOI-700 e "corrected" by 116 min).
- **Calibration:** 11 of 13 known planets are `consistent` (pi Men c, LHS 3844 b,
  TOI-700 b–e, TOI-270 b, WASP-18 b, the benchmark's two signals; all within 1σ).
  TOI-270 c and d, which have known TTVs, are flagged with the timing note. TIC
  111784356, whose dip lasts longer than its star allows, is `density_mismatch`
  (3.1σ). The reviewed candidates TIC 302453503, 14566446 and HD 229793 are
  `consistent`.

### 6b. Per-source depth map (clean candidates)

```bash
.venv/bin/python depth_map_check.py results/<name>/enriched.csv
```

`depth_map_check.py` (method in `depth_map.py`) measures the in-transit flux deficit in
every pixel of TESScut cutouts, then asks which star's image explains that map best:
the target, or any neighbour within 2.5′ up to 5.5 mag fainter. The chi-square of each
source is summed over sectors. It uses the transit fit's ephemeris when the fit is
consistent (over a long baseline the search's period can drift by tens of minutes).

- **Verdicts** (`depth_map/depth_map.csv`): `on_target` (the target within Δχ² 9 of the
  best; others that fit as well are listed), `off_target` (another star better by over
  25 and best in at least half the sectors; flag `depth_map_off_target`, and the FPP
  isn't run), `inconclusive`, `not_detected` (no source beats "no dip" by Δχ² 16: too
  shallow to locate, not evidence against the target).
- **Leave one sector out:** if dropping any single sector flips the verdict between on
  and off target, it's `inconclusive` and the detail names the sector (one sector of
  scattered light can dominate the sum; leave it out by hand with `DEPTHMAP_SKIP`).
- **Why:** `pixel_check.py`'s centroid is pulled about in crowded fields. In b023 it
  called a blend "consistent" (centre 28″ off), and TRICERATOPS gave it NFPP 0; the
  depth map put the dip on a brighter star 82″ away (Δχ² 1,260, best in 4 of 5 sectors).
- **By hand:** `investigations/depth_map.py TIC P T0 T14_HOURS OUTDIR` runs the same
  map with any ephemeris (`DEPTHMAP_SKIP=63` leaves a sector out).
- About 1–3 minutes per candidate (cutouts are cached).

### 7. False-positive probability (consistent fits)

```bash
.venv/bin/python fpp_check.py results/<name>/enriched.csv     # [--draws 1000000] [--runs 10]
```

For candidates whose fit is `consistent`, `fpp_check.py` runs TRICERATOPS
(Giacalone et al. 2021; `vetting/fpp_triceratops.py`, in `.venv-vet`).

- **Inputs:** the same normalised transits, folded and binned, with SPOC and
  TESS-SPOC apertures. QLP publishes none, so QLP-only data can't be used. There's
  no contrast curve.
- **Practice** (after Eschen & Kunimoto 2024):
  - TRICERATOPS is a Monte Carlo estimate, so it runs 10 times (`--runs`) and FPP and NFPP
    are the means; `detail` gives the spread.
  - The fold is binned to 10 minutes, or to the longest cadence if longer (mixed 2-min
    and 30-min data are binned to 30 min), and the model is integrated over the bin.
  - Every nearby star is kept as a possible source. Only after a per-source pixel depth
    map (`investigations/depth_map.py`) has ruled the others out, rerun with
    `--target` and `--neb-radius`. Not automatic: `pixel_check.py`'s centroid is pulled
    about in crowded fields (a b023 blend: centre 28″ ± 16″, true source 82″ away).
  - A low NFPP doesn't clear a candidate: TRICERATOPS's aperture model can dismiss a star a
    few pixels away whose light does reach the pixels. That blend got NFPP 0 with every star
    allowed; only the depth map caught it.
  - Above 8 R⊕ it isn't run (`too_large`): its priors are for planets.
- **Output:** FPP, NFPP and the main scenarios, with a label: `validation_level`
  (FPP < 0.015, NFPP < 0.001), `likely_planet` (FPP < 0.5, NFPP < 0.001),
  `likely_nearby_fp` (NFPP > 0.1) or `unclear`. Written to `fpp_check/fpp.csv`.
- **Network:** it queries Gaia and MAST, and failures are retried 3 times.
- **It's a record, not a verdict:** without high-resolution imaging, the
  unresolved-companion scenarios stay open. TIC 302453503: FPP 0.04–0.06, NFPP
  ≤ 0.0003, `likely_planet`.

Both steps run in `pipeline.sh` when `.venv-vet` exists. `tessdb.py` loads them into
`transit_fits` and `fpp_checks`, and `candidates_view` shows them in its `density`
and `fpp` columns.

## Watchlist: re-testing candidates as new sectors arrive

A candidate found in some sectors predicts a dip of known depth and duration at a known
time in every later sector. `watchlist.py` keeps those candidates and tests the
prediction when new data appears. This is more sensitive than searching again, because
the period doesn't have to be found.

```bash
.venv/bin/python watchlist.py add "TIC 123456789" --run b011 --candidate 1 --reason "doubtful"
.venv/bin/python watchlist.py check          # every watched star; or name some
.venv/bin/python watchlist.py list
.venv/bin/python watchlist.py remove "TIC 123456789"
```

- **add** takes the period, epoch, duration and depth from the run's transit fit, or from
  the search if there was no fit. The sectors the candidate was found in become its
  baseline.
- **check** asks MAST which sectors exist and downloads any not yet measured, using the
  search's products (SPOC, then TESS-SPOC, then QLP). It then measures each transit's
  depth at the predicted time, with the shape fixed and a local quadratic baseline. The
  error is the larger of the formal error and the scatter of the same measurement at
  times away from the transits, which accounts for red noise.
- **Baseline sectors** are measured too, to show the method sees the signal, but only
  later sectors count towards the verdict:
  - `supports`: new depth > 0 at 3σ and consistent with the prediction.
  - `contradicts`: new depth at least 3σ below the prediction.
  - `inconclusive`: neither yet; the report estimates how many more transits are needed.
  - `no new data`: no transit in a later sector yet.
- **Report:** `db/watchlist/<date>.md`, with event plots, plus the sectors tess-point
  expects to observe the star next. Results are in the `watchlist`, `watch_sectors`
  and `watch_events` tables; `rebuild` keeps them.

Run `check` after new sectors are released. Light curves reach MAST a month or two
after each sector, and @MAST_News announces ingests. Sectors 97 and 98 were four
orbits long, and "sector 1751" was a short special pointing (interstellar comet
3I/ATLAS), so sector numbers aren't a timeline.

## Results database and browsing

`tessdb.py` keeps one SQLite file, `db/tess.db`, holding every run. It isn't in git, so back
it up as described below. Runs are loaded after they finish; the per-worker CSVs stay as
the raw record.

```bash
.venv/bin/python tessdb.py ingest results/aries        # after enrich and pixel checks
.venv/bin/python tessdb.py note "TIC 123456789" "Possible new eclipsing binary..." --tag vsx-candidate
.venv/bin/python tessdb.py list                        # runs, signals, clean counts
.venv/bin/python tessdb.py backup                      # dated copy in db/backups/ + records/*.csv
.venv/bin/python tessdb.py restore-records             # load records/*.csv into a (new) database
.venv/bin/python tessdb.py remove test                  # delete a run (named after its folder)
.venv/bin/python tessdb.py rebuild                     # after a schema change: re-ingest all, keep notes
.venv/bin/python tessdb.py clean-cache [--dry-run]     # delete cached light curves of no-signal stars
```

- **What `ingest` loads:**
  - a folder's searches, candidates and flags (flags are split out so they can be queried);
  - its `enriched.csv`, `pixel_check*/` and `injections.csv`, where present;
  - its `run_info.json`, which `run_parallel.sh` writes at the start of every run with the
    git commit, any uncommitted changes, the target file and the search settings.
- **Re-ingesting** a run replaces its rows. For folders made before `run_info.json`
  existed, pass `--commit`, `--targets`, `--started` and `--note`.
- **Duplicates fail loudly.** Each check may appear once per run: one pixel check per
  star, candidate and sector, and one pixel verdict, reduction check, transit fit and
  FPP result per candidate (also one row per injection job). If a run's folders hold the
  same check twice (for example an extra `pixel_check_*/` that repeats a candidate),
  ingest stops and names the duplicate row. The transaction is rolled back, so the
  database keeps what it had. Remove the duplicate and ingest again.
- **Records in git.** Most of the database can be recreated by re-ingesting the run
  folders, or by rerunning the batches. Notes, the watchlist and the batch list can't:
  they record decisions. `backup`, run at the end of every pipeline, also writes those
  tables to `records/*.csv` (notes, watchlist, watch_sectors, watch_events, batches),
  which are tracked in git, so a commit and push puts them on GitHub. After losing the
  database: `rebuild` or re-ingest the runs, then `tessdb.py restore-records`. Notes
  already present are skipped. The dated database copies in `db/backups/` stay on this
  disk (about 7 MB each, 3–4 MB gzipped), so copy one elsewhere now and then.
- **Testing:** set `TESSDB=scratch/test.db` so test runs don't go into `db/tess.db`.
  To remove a run that does, use `tessdb.py remove <run>`.
- **Notes** record decisions and follow-ups against a star, such as VSX submissions or
  blend sources. They appear beside the star's candidates in `candidates_view`.

**Browsing in a browser on your PC:**

```bash
tmux new -s db '.venv/bin/python tessdb.py serve'      # Datasette on 127.0.0.1:8001
```

- **Opening it:** VS Code Remote-SSH forwards port 8001 automatically (see its Ports
  panel); otherwise add `-L 8001:localhost:8001` to your ssh command. Then open
  `http://localhost:8001`.
- **Using it:** start with `candidates_view`. It lists every signal with its run, the
  worker that found it, flags, notes and a clickable `plot_url` that opens the diagnostic
  plot. **Restart the server after ingesting a new run**, because plot folders are mounted
  when it starts. Tables can be filtered and faceted (by status or flag), exported as
  CSV, or queried with SQL. `stars_searched` shows which stars have been searched and when.
- **Security:** only the results folders are served, for plots. The repository root,
  `.git` and the venv are not, and the server listens only on localhost.

For a quick look at one folder without the database, plain Python is enough:
`python3 -m http.server 8002 --bind 127.0.0.1 --directory results/aries`.

## Output files

| File | Contents |
|---|---|
| `results/<name>/candidates_all.csv` | One row per detected signal: period, epoch, depth, SNR, per-sector depths and timing, status, flags |
| `results/<name>/done_all.txt` | Each star's outcome; errors are retried on a rerun |
| `results/<name>/w<k>/run.log` | Detailed search log for worker *k* |
| `results/<name>/w<k>/plots/<star>/` | Per-candidate plot: periodogram, folded light curve, odd vs even transits |
| `results/<name>/enriched.csv` | Candidates plus catalogue information and catalogue flags |
| `results/<name>/pixel_check/` | Pixel-check tables (per sector and per candidate) and difference-image plots |

## Status and flags

Each detected signal gets one of three statuses:

- `clean`: no flags raised, apart from informational ones (`v_shaped`); worth a closer look.
- `flagged`: one or more vetting flags.
- `same_period`: related to a signal already found at the same period: a
  `secondary_eclipse`, a `phase_curve_residual` or a `harmonic`.

### Search flags (`flags` column)

| Flag | Meaning |
|---|---|
| `alias_of_<P>` | Period is a simple ratio of a signal already found |
| `few_transits` | Fewer than 5 usable transits; dip counts mean little |
| `transits_at_data_edges` | Most transits sit against a gap in the data, where detrending is weakest |
| `odd_even_mismatch` | Alternate dips differ in depth |
| `double_period_variable` | Comes with an odd/even mismatch: the true period is probably 2P (`period_x2_d`). Usually an eclipsing binary, or two spot groups on a star rotating in 2P |
| `single_event_dominated` | One event carries most of the signal |
| `sinusoidal_variability` | A sine wave fits better than a box-shaped dip |
| `long_duration` | Duration is more than 10% of the period. Used only when the star's TIC radius and mass are missing; otherwise `duration_too_long_for_star` applies |
| `non_positive_depth` | Not a dip |
| `too_large_for_planet` | Implied radius above 2 Jupiter radii |
| `duration_too_long_for_star` | Longer than any transit the catalogued star allows. Uses the fitted duration (`dur_fit_h`) with a 1.3× margin; the coarse search duration with 1.5× if the fit fails. The shape fit starts from 11 durations and ingress shapes and keeps the best, because a single start at the search duration can stop in a short local minimum (TIC 111784356) |
| `depth_varies_by_sector` | Depth differs between sectors; typical of a blend, since each sector's aperture lets in different amounts of a neighbour's light. Each sector is measured at its own dip time (a linear drift fitted to the sector timings, since the search's period can be slightly off over a long baseline), and its depth is corrected for exposure smearing, using the dip's own shape from the shortest-cadence data: a 30-min exposure makes a short transit look shallower. Sectors keeping less than half the dip are shown as `(smeared xF, not compared)`; corrected ones as `(/F)`. TIC 12999193 (TOI-7689.01, 0.6 h) was wrongly flagged before this |
| `dip_missing_in_sector` | The dip should be visible in some sector but isn't |
| `timing_varies_by_sector` | Dips wander in time; eclipses and transits keep time, spot minima don't |
| `significant_secondary` | A dip at phase 0.5 above 5σ and deeper than 10% of the main dip or 1000 ppm: a stellar companion's second eclipse (depth in `secondary_ppm`) |
| `near_momentum_dumps` | More than half the transits (at least 3) lie within half a duration plus 3 h of a momentum dump, a likely pointing artefact. Counts per sector are in `dump_near` |
| `asymmetric_dip` | The two sides of the dip (just outside, or the two halves inside) differ by more than 3σ and 25% of the depth (`asym_sigma`) |
| `ambiguous_period` | The periodogram's best peak at an unrelated period (not a ratio m/n of the candidate's, m, n ≤ 5, nor of another signal found in the star) is within 1.0 SDE (`sde_margin` = margin@rival period). Noise dips can be strung together at several periods almost equally well (TIC 14566446: 10.79 vs 10.51 at 19.10 d). Calibrated on 65 recovered injections (1–3 R⊕; smallest margin 1.06) and the known planets: none flagged |
| `not_unique` | Model-shift test, after Kepler's Robovetter: folded at the period, another dip (or a brightening) is about as significant as the candidate, allowing for the number of independent phases (`model_shift` = ms2;ms3, flagged when either is below 0). Calibrated on 11 known planets and 50 recovered injections: none flagged |
| `fit_duration_too_long` | Added by `transit_fit.py` (step 6): the fitted T14 is over 1.3× the longest transit the TIC star allows. Such a fit can pass the density test through wide errors (a b017 fit wandered to 9.9 h; a b023 blend fitted 7.7 h against 4.4 h). Real candidates so far reach 1.16× |
| `density_mismatch` | Added by `transit_fit.py` (step 6): the transit's stellar density differs from the TIC star's by more than 3σ. A blend or eclipsing binary, or (with the timing note) timing variations smearing the fit |
| `odd_even_timing` | Odd and even transits are centred at different times: more than 3σ and 15% of the duration (`oe_dt_min`) |
| `v_shaped` | **Informational: doesn't stop `clean`.** The trapezoid fit's ingress takes more than 40% of the half-duration (`ingress_frac`; 1 is a pure V), and a flatter-bottomed fit is worse by Δχ² > 9 (scaled by the fit's reduced χ²). Grazing: an eclipsing binary, or a planet at high impact parameter (TOI-7689.01 is one) |

### Catalogue flags (`enrich_flags` column)

| Flag | Meaning |
|---|---|
| `nearby_variable_period_match` | A known variable nearby shares the period (or a simple multiple of it); it is almost certainly the source |
| `known_planet_period_match` | The target is a catalogued planet host at this period |
| `target_known_variable_period_match` | The target itself is a catalogued variable at this period |
| `target_catalogued_variable` | A catalogued variable lies at the target's own position, with or without a period (Gaia's `SOLAR_LIKE` class is listed in `target_variable` but doesn't raise the flag) |
| `has_toi` / `has_ctoi` | ExoFOP lists a TOI or community TOI |
| `catalogue_check_incomplete` | A catalogue query (ExoFOP, SIMBAD, Gaia, MAST, a VizieR catalogue) still failed after 3 attempts with timeouts. The failure is never cached as "nothing found": rerun `enrich_candidates.py` and only the missing checks are fetched again. The step's log lists the stars and checks affected. Until then, a missing `has_toi` or variable match means nothing |
| `has_spoc_tce` | SPOC's own pipeline found this period (or 2× or ½×) as a threshold-crossing event, from its data-validation reports on MAST (the widest-range and latest reports). Known to SPOC even if never promoted to a TOI (one b020 candidate was such an unpromoted TCE). Details in `spoc_tces` |
| `simbad_<type>` | SIMBAD classifies the target as something notable, such as `SB*` |
| `gaia_ruwe_high` | RUWE above 1.4: extra wobble on the sky, often an unresolved companion |
| `gaia_non_single_star` | Gaia fitted a binary (or companion) orbit |
| `gaia_rv_variable` | Radial velocity varies (Katz et al. 2023: p ≤ 0.01, renormalised gof > 4, ≥ 10 epochs) |
| `gaia_bp_rp_excess` | Corrected BP/RP flux excess beyond 3σ (Riello et al. 2021): possible blend within a few arcseconds |
| `gaia_multi_peak` | More than 10% of Gaia scans show a double image |

The Gaia flags show that the star has a companion or changing velocity. They don't by
themselves mean the dip is a false positive: WASP-18 raises `gaia_non_single_star` and
`gaia_rv_variable` because its 10-Jupiter-mass planet makes the star wobble. The values
behind them are in the `gaia_ruwe`, `gaia_rv_chisq_pvalue` and `gaia_c_star_sigmas` columns.

## Tips

- **Controls:** check the pipeline on `targets/controls.txt`, a list of stars with known planets
  (π Men, LHS 3844, TOI-700, TOI-270, WASP-18). All 10 planets should come out `clean`.
- **Regression set (after any major change):** `targets/regression.txt` is the control
  stars plus stars whose outcome is known: the benchmark (two signals), the candidate
  TIC 302453503 (must stay clean), and three reviewed false or doubtful signals (must
  not turn clean). Run it, then check against `tests/regression_expect.csv`:

  ```bash
  ./run_parallel.sh targets/regression.txt scratch/regression 8
  .venv/bin/python tests/check_regression.py scratch/regression
  ```

  This prints one line per expected signal and exits 1 if any fail (about 20 minutes).
  It is kept separate from `targets/controls.txt`, which lists only confirmed planets.
- **Period-grid warning:** "period grid capped" means the data span is long enough that
  the trial periods are coarser than intended. The shortest transits lose a little
  sensitivity. Raising `MAX_GRID` fixes it at the cost of run time.
- **Lessons from the first field:** every candidate that reached the catalogue and pixel
  checks turned out to come from a neighbour 30–80″ away. Most were blends with
  catalogued eclipsing binaries, caught by step 3. One was a chance alignment with the
  1.38-day spot cycle of an uncatalogued variable (TIC 434331752), caught by step 4.
  TESS pixels are 21″ across, so always check neighbours well beyond the target itself.
- **QLP quality flags:** the search drops every QLP cadence with a non-zero quality flag
  (`QLP_STRICT`). QLP's flag bits aren't covered by lightkurve's default mask, and in
  some sectors up to half the cadences are flagged. Those points carry most of the
  outliers and once produced fake 4%-deep transits.
- **Edge trimming:** the search always trims at least 0.5 days from each end of every
  data segment. It trims up to 3 days while the edge still looks disturbed, because
  scattered-light ramps in QLP data can last a day or more.

## Files

| File | Purpose |
|---|---|
| `tess_target_finder.py` | Step 1: target selection |
| `run_parallel.sh` | Step 2: parallel runner (calls `run_parallel.py`) |
| `run_parallel.py` | Supervises the workers: shared queue, retries, crash recovery, pausing, merged output |
| `workqueue.py` | The run's queue (`queue.db`), used by the supervisor and the workers |
| `mastutil.py` | Shorter MAST timeouts (90 s per request) for the per-star scripts |
| `prefetch.py` | Fetches a batch's MAST data ahead of its search (`prefetch.jsonl`), so the search runs without MAST |
| `transit_fit.py` | Step 6: physical transit fit and density check (runs `vetting/transit_mcmc.py` in `.venv-vet`) |
| `depth_map.py` | Per-source pixel depth map from TESScut cutouts (used by the two below) |
| `depth_map_check.py` | Step 6b: depth map verdict for clean candidates |
| `followup.sh` | Steps 3–7 plus ingest and backup for one batch, started by `pipeline.sh` in a `tess` window when the search is done; one at a time |
| `fpp_check.py` | Step 7: TRICERATOPS false-positive probability (runs `vetting/fpp_triceratops.py` in `.venv-vet`) |
| `vetting/` | Scripts that run in the vetting environment (`requirements-vet.txt`) |
| `find_planets.py` | Step 2: transit search and vetting |
| `enrich_candidates.py` | Step 3: catalogue checks |
| `watchlist.py` | Re-tests watched candidates at their predicted transit times as new sectors are released |
| `vsx_shortlist.py` | Step 3b: variable-star signals for the weekly VSX session (`vsx_shortlist.csv`) |
| `pixel_check.py` | Step 4: pixel-level checks |
| `inject_recover.py` | Sensitivity test: injects synthetic transits and measures how many are recovered |
| `reduction_check.py` | Step 5: is the dip in every independent reduction? |
| `pipeline.sh` | Runs steps 2–5 and the database load, unattended |
| `tessdb.py` | Results database: ingest runs, notes, backups, Datasette browsing |
| `tests/` | `test_run_parallel.sh`: fast checks of the run logic with a fake search; `check_regression.py` + `regression_expect.csv`: expected outcomes for `targets/regression.txt` |
| `investigations/` | Close looks at individual stars: scripts, plots and notes (index in its README) |
| `targets/` | Target lists (`.txt`, read by the search) and the finder's ranked CSVs: `controls.txt` (known-planet control stars), `ffi.txt` (first field), `aries.txt`, `b001.txt`, ... |
| `records/` | CSV copies of the tables only the database holds (notes, watchlist, batches), written by `tessdb.py backup`; commit them after each batch |
| `server/` | Boot-time setup: the tmux sessions (claude, db, jupyter, tess) recreated by a systemd user service, and `batch_status.sh` for batches a reboot interrupted |
| `publish/` | Private-only: `sync_public.py` copies the allowlisted code (`allowlist.txt`) into the public PlanetFinder repo, refusing if a file names an unpublished star |
| `docs/` | `literature_review.md` (how comparable pipelines vet; improvements to make) and `references.bib` (citations for a write-up) |
| `STATUS.md` | Current state: progress, candidates, VSX queue, open items. Updated after each batch review |

Not in git (made by running the pipeline):

| Folder | Contents |
|---|---|
| `results/<name>/` | One folder per run (search, checks, logs, plots); its database run name is `<name>` |
| `db/` | `tess.db`, dated backups in `db/backups/`, Datasette settings in `db/datasette/` |
| `scratch/` | Throwaway tests and logs |

Earlier versions of the search script are in the git history at tag `pre-cleanup`: the original
`find_planets.py` (v1), `find_planets_v2.py`, and `find_planets_v3.py`, which is today's `find_planets.py`.
For example: `git show pre-cleanup:find_planets_v2.py`.
