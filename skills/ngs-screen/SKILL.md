---
name: ngs-screen
description: >-
  Process pooled sensor-screen sequencing data: paired-end FASTQs to per-guide
  counts, CRISPResso2 editing rates, lane- and replicate-merged editing
  tables, and a guide count matrix. Covers base-editing (15 nt barcode / 42 nt
  sensor) and prime-editing extended-barcode (16 nt matching / 55 nt sensor)
  layouts. Use for: process my screen sequencing data, make a config file for
  a run, run CRISPResso on the sensors, rerun the array tasks that failed,
  merge NGS lanes or replicates, compute editing percentages, build a count
  matrix, debug zero counts or low read assignment. Not for LFC/FDR hit
  calling (screen-stats) or arrayed single-locus amplicons (genotyping).
---

# ngs-screen — FASTQ → counts → editing tables → count matrix

## Where this skill starts and stops

```
BMC datahub ──► run folder ──► STEP4 ──► STEP5 ──► STEP6 ──► STEP7 ──► STEP8 ──┐
                 counts/        per-guide  CRISPResso  aggregated  editing    count   │
                 classification/ FASTQs    per guide   CSV/sample  tables    matrix   │
                                                                                      ▼
                                                              screen-stats: LFC, z, p, FDR
```

**It ends at two files.** The MLE editing tables (`MLE/<condition>.csv`) and the
count matrix (`*_COUNTS.txt`). LFC, z-scores, empirical p-values, FDR and hit
calling are the **screen-stats** skill, which reads exactly those two. Do not
compute an LFC here.

Neighbours: **library-design** designs the library this skill validates;
**genotyping** is for arrayed single-locus amplicons (same tools, one locus, no
barcode pool); **sci-plot** styles any figure.

**MAGeCK is deliberately not run.** The count matrix is MAGeCK-shaped and the
lab's `analysis.md` has the invocation, but MAGeCK was dropped for these screens:
it collapses guides into gene-level scores and consistently overcalled hits,
including in noisy tissues. The live path is screen-stats' non-parametric
control-guide null. If someone asks for MAGeCK, say this first.

## This is not a laptop task

STEPs 4–6 need a cluster. STEP4 is a pure-Python per-read loop that does a
full-library boolean DataFrame scan **per read**, plus an O(library) Hamming
fallback per unmatched read. Budget tens of minutes to a few hours per sample on
a multi-million-read FASTQ. It is embarrassingly parallel, so as a SLURM array
the wall clock is one sample, not 188. Memory stays under a few GB; reads are
streamed. STEP5 is thousands of independent CRISPResso invocations (minutes to an
hour per sample). STEP6 is I/O-bound, a minute or two per sample.

STEPs 7–8 are local and cheap — seconds to a couple of minutes for a 48-lane
screen. The disk hog is one directory per guide per sample; put it on the lab
volume, not home.

## Inputs

**1. A run folder**, one subdirectory per NGS sample, each with a pair of
uncompressed FASTQs named `*_1_sequence.fastq` (R1 = sensor read) and
`*_2_sequence.fastq` (R2 = protospacer read). Mounted locally this is
`/Volumes/sanchezrivera/<user>/<run>`; on Luria it is `$LABROOT/<run>`, with
`$LABROOT` defaulting to `/net/bmc-lab2/data/lab/sanchezrivera/$USER`.

Stage it first — this is step one of every real run:

```bash
ssh <you>@luria.mit.edu
srun --pty bash
srun rsync -av /net/bmc-pub17/data/bmc/public/datahub/<datafolder> $LABROOT/
gunzip $LABROOT/<run>/*/*.fastq.gz      # see the GZ gotcha below
```

**2. A library CSV** with these exact column names:

| column | BE | PE | notes |
|---|---|---|---|
| `gRNA_id` / `pegRNA_id` | `gRNA_id` | `pegRNA_id` | guide id, becomes `Guide_ID` |
| `Protospacer` | 20 nt | 20 nt | leading G + 19; **must be 20**, see gotchas |
| `Hamming_BC` | 15 nt | 16 nt | sense strand; script revcomps it |
| `sensor_wt` | 42 nt | 55 nt | `--amplicon_seq`; sets `sensor_len` |
| `sensor_alt` | 42 nt | varies | `--expected_hdr_amplicon_seq` |
| `wt_qwc`, `edit_qwc` | — | required | per-pegRNA quantification window |
| `gene_name_m`, `classification`, `Editor` | for STEP8 | for STEP8 | gene names, control collapse, subpool filter |

`bc_len` and `sensor_len` are **not** parameters you tune — they must equal
`len(Hamming_BC)` and `len(sensor_wt)`. `setup_run.py` enforces this.

## The two layouts

| | BE | PE extended-barcode |
|---|---|---|
| `bc_len` (matching) | 15 | **16** |
| `real_bc_len` (sensor start) | 15 | **8** |
| `sensor_len` | 42 | 55 |
| CRISPResso | `--base_editor_output`, `--plot_window_size 5` | `--prime_editor_output`, `--plot_window_size 15` |
| quant window | hardcoded `6-35,6-35` | per-row `wt_qwc,edit_qwc` |

The PE layout exists because earlier PE screens showed unusually high
barcode–protospacer recombination. The matched barcode is extended 8 nt into the
sensor for specificity: R1 is `[8 nt real BC][55 nt sensor]`, the first **16 nt**
are matched against the library, but the sensor extracted into the per-guide
FASTQ starts at position **8**. Two different numbers, which is exactly where the
source had a bug.

## Procedure

### 1. Stage and validate (local, seconds)

```bash
python scripts/setup_run.py \
    --run-dir /Volumes/sanchezrivera/$USER/251030Hem \
    --library FINAL_focused_library.csv \
    --layout BE -o . --make-folders
```

Writes `CONFIG_<run>.txt`, the three sbatch scripts, and creates `counts/`,
`classification/`, `confusion_mats/`, `crispresso/`. It prints the library
report and the `--array=1-N` you need, and **refuses to write anything** if the
library has a problem that would corrupt the run. Read that report.

Then confirm the barcode strand on one sample before spending the array:

```bash
python scripts/qc_screen.py strand --library LIB.csv \
    --r1 <sample>/..._1_sequence.fastq --r2 <sample>/..._2_sequence.fastq
```

### 2. Copy to the cluster

Into `$LABROOT/<run>/`: the config, the library CSV, the three `.sh` files,
`sensor_extraction.py`, `crispresso_sensors.py`, and
`assets/crispresso_quant_blank.csv` (**STEP6 fails without it**).

### 3. Run the three arrays (cluster, hours)

```bash
cd $LABROOT/<run>
sbatch step4_sensor_extraction.sh     # sensor_env
sbatch step5_crispresso.sh            # crispresso_env, after STEP4
sbatch step6_aggregation.sh           # crispresso_env, after STEP5
```

Between steps, check what finished:

```bash
python scripts/qc_screen.py failed --run-dir . --config CONFIG_<run>.txt
python scripts/qc_screen.py report --run-dir .
```

`failed` prints a ready-to-paste `--array=3,17,58` for just the samples that
died. Edit that into the `.sh` and resubmit — do not rerun the whole array.

### 4. Copy back and finish locally (local, minutes)

Copy `counts/` and `crispresso/*_crispresso_aggregated.csv` to your machine.

```bash
# the sample manifest: raw STEP4 folder name -> friendly id with a lane suffix
python scripts/compile_editing.py samples CONFIG_<run>.txt > samples.tsv
$EDITOR samples.tsv          # fill in: input-1, input-2, d15-rep1-1, spleen1-2, ...

# STEP7a: aggregated CSVs -> per-sample editing tables, lanes merged
python scripts/compile_editing.py compile --samples samples.tsv \
    --crispresso-dir ./crispresso --out-dir ./compact

# STEP7b: merge biological replicates into conditions
printf 'bm,bm1,bm2,bm3,bm4,bm5\nspleen,spleen1,spleen2\nd5,d5\ninput,input\n' > groups.csv
python scripts/compile_editing.py mle --compact-dir ./compact \
    --groups groups.csv --out-dir ./MLE

# STEP8: the count matrix, one subpool at a time
python scripts/compile_editing.py counts --samples samples.tsv \
    --counts-dir ./counts --library LIB.csv \
    --out ABE_BC_COUNTS.txt --count-col bc_count --editor ABE
```

`MLE/<condition>.csv` and `ABE_BC_COUNTS.txt` are what screen-stats reads. Stop
here.

## Environment

Two micromamba envs on the cluster, specs in this repo's `envs/`:

- **sensor_env** (STEP4) — numpy, pandas, biopython, python-levenshtein
- **crispresso_env** (STEP5/6) — python 3.9/3.10, `crispresso2` **pinned**

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate sensor_env
```

micromamba replaced conda because the conda solver hung for hours.

STEP7/8 are local and need only pandas + numpy. **Do not run them in the pegg
env** — pegg's Python 3.9/3.10 + `scikit-learn==1.1.1` + `numpy<2` pin belongs to
library-design and buys you nothing here.

Pin CRISPResso2 and record the version in your methods: editing percentages are
version-dependent, so two people running the same reads months apart can get
different numbers with nothing in the output saying so.

## Gotchas

These came from reading the source. The first three are **real bugs** that this
skill's scripts fix; the fixes are marked `# FIX:` inline so they are visible
rather than silent.

**1. The BE extraction script cannot start.** `sensor_extraction_42nt_trimming_quality.py`
checks `len(sys.argv) != 10` and then reads `sys.argv[10]` — which does not exist
when argv has 10 entries. An unconditional `IndexError` before any work. (The PE
sibling is correct: `!= 12` then `argv[11]`, because `-o` takes a slot.) Fixed
here by using argparse, so the off-by-one cannot come back. The usage string it
printed on failure also omitted `quality_check_mode` entirely.

**2. The BE script never applied your `quality_check_mode`.** It called
`quality_checker(q1, q2)` while the signature is
`quality_checker(q1, q2, bc_len, sensor_len, proto_len=19, mode=...)` — two
required positional args missing, so `TypeError` on the first read. Worse: if
anyone "fixed" that by giving those params defaults, the mode read from argv
would still be ignored and **every run would silently use `region_average`**.
Fixed: all four passed. The PE sibling already did this right.

**3. Every PE sensor was extracted from the wrong window.** In
`sensor_extraction_peg_counts_ext_bc_16nt.py`, `to_IOSeq_rec(...)` slices
`sensor_seq[real_bc_len : real_bc_len+sensor_len]`, but all three call sites pass
`bc_len` (16) into that parameter instead of `real_bc_len` (8). Every PE sensor
read therefore started 8 nt downstream and ran 8 nt past the window, so the
per-guide FASTQs did not align to the 55 nt `sensor_wt` and CRISPResso reported
collapsed or zero alignment. Fixed and verified: the extracted read now equals
`sensor_wt` exactly. **Any PE sensor editing rate produced before this fix is
wrong.**

**4. The count matrix mislabelled every sample.** `counts.ipynb`'s `count_matrix`
did `new_sample_ids = np.unique([x.rsplit('-',1)[0] for x in samp_id])` —
alphabetically sorted — while `samps_to_merge` is in file order, then paired the
i-th lane pair with the i-th *alphabetical* name. On the real 46-lane sample set
(`lib, input, d5, d15-rep1..5, spleen1..5, bm1..5, men1..5`) **all 23 lane pairs
get the wrong label**: the library's counts land under `bm1`, input's under
`bm2`. Fixed with `dict.fromkeys`, which is what STEP7's `merge_NGS_replicates`
already used. **Rebuild any count matrix made with the notebook version.**

**5. `Protospacer` must be exactly 20 nt.** The scripts compare
`Protospacer[1:]` against `r2[:19]` using `Levenshtein.hamming`, which counts a
length difference as mismatches. With a 21 nt protospacer the best achievable
distance is 1, so `proto_identified_perfect` stays at **zero for the whole run**
and nothing matches. The lab's real libraries are 20 nt (leading G + 19).
`setup_run.py` and `sensor_extraction.py` both refuse to proceed otherwise.

**6. Compressed FASTQs read as zero matches, with no error.** `GZ=False` is a
module global in every extraction script and is not on the command line. Gunzip
first, or pass `--gz`. `setup_run.py` warns when it sees `.gz` in a run folder.

**7. Reruns double-count.** Per-guide FASTQs are opened in append mode `'a'`, and
the BE original began with a bare `os.mkdir()` that raises `FileExistsError`, so
a failed array task could neither resume nor safely rerun. Here a populated
folder is **refused** with an explanation; `--force` clears and redoes it. Never
rerun STEP4 into a folder that already has reads.

**8. The four output folders must pre-exist**, lowercase and exactly spelled:
`counts/`, `classification/`, `confusion_mats/`, `crispresso/`. The scripts write
straight to `'counts/' + ...` with no `makedirs`, so the task dies on the final
`to_csv` **after** all the compute. Use `--make-folders`.

**9. `Hamming_BC` is reverse-complemented before matching.** If your library
already stores the R1-strand barcode, every read lands in `no_match_bc` and
counts come back all zero. Check one sample with `qc_screen.py strand` before
launching — it reports the match rate for both strands and names the culprit.

**10. One missing CRISPResso directory used to lose a whole sample.** STEP6 did
`os.listdir(.../CRISPResso_on_<guide_id>)` with no existence check, so a single
guide without an output directory raised `FileNotFoundError` and killed the
aggregation — after the CRISPResso compute was already paid for. The blank-row
fallback only covered a *missing file in an existing directory*. Here both cases
take the all-zero blank rows and the substitution count is printed.

**11. The three count columns are not interchangeable.**
`matched_guide_count` (barcode *and* protospacer correct — stringent),
`bc_count` (barcode correct regardless of protospacer, includes recombinants —
highest, and the notebook default), `total_guide_count` (protospacer correct
regardless of barcode). Recombination rate sets how far apart they are. Pick one
and keep it fixed across a paper. For a duplicated protospacer
(`duplicate_sgRNA=True`), `total_guide_count` is inflated because the increment
hits every row sharing that protospacer.

**12. The BE quantification window is hardcoded `6-35,6-35`** — 5 nt flanks
around a 20 nt protospacer in a 42 nt sensor. Change the sensor design and it
silently quantifies the wrong bases. For PE the window lives per-row in the
library instead.

**13. Lane merging pairs consecutive rows blindly.** If any sample has one lane
instead of two, the pairing shifts and reads from two different biological
samples get merged with no error at all. `compile_editing.py` refuses unless
every sample has exactly two lanes named `<sample>-1`/`<sample>-2`; use
`--no-merge-lanes` and group explicitly if your run is genuinely uneven.

**14. Subpool hygiene — the one thing the notebook shouts about.** An ABE count
matrix must contain only ABE guides. Leaving CBE guides in does not error; it
changes every column's depth normalization and quietly destroys the LFC. Always
pass `--editor`.

**15. Only two exact spellings collapse to `controlsgRNA`:**
`non-targeting control` and `safe-targeting control`. Any other spelling leaves
each control as its own single-guide "gene" and screen-stats has no null to
build. `setup_run.py` flags near-misses.

**16. `byproduct_ambiguous` can go negative.** It is a residual
(`byproduct_all - sub_all - indel_all`) and the lifted code only `fillna(0)`s, it
does not clip. Filter on `Reads_aligned_all_amplicons` before interpreting
byproduct splits.

**17. `MLE_merge_auto` sums raw counts then recomputes percentages**, so a
replicate with 10× the reads carries 10× the weight. That is intentional — it is
a maximum-likelihood merge — but one deep replicate can dominate a condition. It
asserts `Guide_ID` order matches and fails loudly if subpool filters differed.
(The legacy `MLE_merge` handles exactly 3 replicates; it is not lifted here.)

**18. The shipped sbatch files contain `srun --pty bash` in the body** — a
copy-paste artifact from the interactive setup notes, which in a batch script
opens a shell the job cannot use. They also hardcode a run path, a library name
and `kexindon@mit.edu` or `XXX@mit.edu`. The generated scripts omit all of it.

**19. `--array=1-N` must equal the number of config rows.** Too small silently
skips samples; too large gives tasks whose awk lookup returns an empty
`R1_FILE` and fails confusingly. `setup_run.py` prints the right N.

**20. CRISPResso can fail on group write permissions** on the shared
filesystem: `chmod -R g+w $LABROOT/<run>/crispresso`.

## Verification

After STEP4, on the first sample, before launching the rest:

```bash
python scripts/qc_screen.py report --run-dir .
```

- **`good_%`** — low means a QC mode mismatch, not bad sequencing. Try
  `region_average` (the BE default was never actually applied, see gotcha 2).
- **`no_match_bc_%`** — should be low. Over ~50% is almost always the barcode
  strand; run `qc_screen.py strand`.
- **`perfect_%`** of `bc_identified` — near zero with a healthy `no_match_bc`
  means the protospacer length is wrong (gotcha 5).
- **`recomb_%`** — over ~20% and the three count columns diverge badly. Read
  `confusion_mats/` off-diagonal mass before trusting PE counts.
- **`zero_count`** and **`top10%_share`** — library skew, which propagates
  straight into LFC noise.

After STEP7, sanity-check the arithmetic on one guide: `WT + correct_edit +`
(bystander + byproducts) should reconcile against
`Reads_aligned_all_amplicons`, and all `_perc` columns are over
`Reads_aligned_all_amplicons`, **not** `Reads_in_input`.

After STEP8: the count-matrix column names must be the biological samples you
expect (gotcha 4 made them wrong for everyone), `controlsgRNA` must be present
with a sensible guide count, and no column should be all zeros. The script warns
on all three.

## Figures

This skill prints tables. Anything headed for a paper — editing-rate
distributions, recombination heatmaps, count-skew panels — hands off to the
**sci-plot** skill. No skill here defines its own figure style.

## Files

```
scripts/setup_run.py           STEP3: config + library check + folders + sbatch
scripts/sensor_extraction.py   STEP4 engine, BE + PE, 4 bugs fixed
scripts/crispresso_sensors.py  STEP5 run + STEP6 aggregate
scripts/compile_editing.py     STEP7 compile/mle + STEP8 count matrix
scripts/qc_screen.py           report | failed | strand
assets/crispresso_quant_blank.csv   ship this to the cluster for STEP6
```

Every engine is lifted verbatim from `fsrlab-sensor-analysis`, with the source
file and cell named at the top of each script so the two can be diffed and
cannot drift. The column lists in `compile_editing.py` are copied character for
character because screen-stats' `attach_editing_info` reads those exact names.
Prose condensed from `fsrlab-sensor-analysis/docs/source/{setup,extraction,analysis}.md`.
