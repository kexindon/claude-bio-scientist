---
name: cohort_analysis
description: Answers how often a mutation occurs in real patients, and what carrying it is associated with. Use for "how common is DNMT3A R882H", "how many AML patients have FLT3-ITD", a bare variant name ("IDH2 R140Q?", "TET2 hotspots"), "pull a B-ALL mutation table", "which genes co-occur with TP53", "is SRSF2 mutually exclusive with SF3B1", "fit a Cox model", "get the hazard ratio", "plot Kaplan-Meier", "redo Figure 2", 突变频率, 多少病人, 共突变, 互斥, 生存分析, 预后, 算一下 HR. Not for figure styling — that is sci-plot.
---

# cohort_analysis — patient cohorts: frequency, co-mutation, outcome

Two kinds of question, split by **where the data lives**:

| Question | Data source | Section |
| --- | --- | --- |
| How often is this mutated? Which cancer types? Which genes co-occur? | public cohorts under `$LAB_DATA_ROOT` | [Public cohorts](#public-cohorts) |
| Does carrying it change survival? HR, KM, log-rank | a **local** patient-level cohort file you pass by path | [Outcome analysis](#outcome-analysis) |

Nothing here styles a figure. Every plotting function returns a **plot-ready
frame** — what goes on which axis, which points are significant — and fonts,
palettes, panel size and PDF export come from the **sci-plot** skill. A quick
QC glance is fine to draw directly; anything going into a paper hands off.

## Where the data lives

Public cohorts are reached through `$LAB_DATA_ROOT` (see `docs/data-layout.md`);
no path is hardcoded. **AACR-GENIE is controlled-access under a data use
agreement** — read it in place, never copy it into a repo, and check with
whoever administers the lab's GENIE access before pointing a shared machine or
an automated agent at it. The patient-level outcome cohort
(`final_data_matrix.RData`) is row-level de-identified clinical data: it stays
outside the repo, you pass its path each time, and this skill never copies,
embeds or commits any of it.

```bash
export LAB_DATA_ROOT=~/Documents/GitHub/Database
```

## Public cohorts

| `--cohort` | what it is | notes |
| --- | --- | --- |
| `genie` | AACR-GENIE, pan-cancer, ~200k samples | `v17.0` only is exercised; GRCh37 |
| `heme` | `heme_msk_impact_2022` | 296 AML samples by `CANCER_TYPE_DETAILED` |
| `impact` | `msk_impact_2017` | registered, never exercised by the source notebooks |

### Variant and gene frequency

```bash
python ~/.claude/skills/cohort_analysis/scripts/cohort_analysis.py freq \
    --cohort heme --disease-detailed 'Acute Myeloid Leukemia' \
    --genes DNMT3A,TET2,ASXL1 --classifications LOF -o OUT
```

```
OUT/variant_frequency.csv   one row per variant: n_patient, n_cohort, freq, ci_low, ci_high
OUT/gene_frequency.csv      genes ranked, with freq and optional COSMIC role
OUT/provenance.json         cohort, version, build, which column was filtered, the denominator
```

MAF column names come out intact (`Hugo_Symbol`, `Chromosome`,
`Start_Position`, `End_Position`, `Variant_Type`, `Reference_Allele`,
`Tumor_Seq_Allele2`, `HGVSp_Short`), so the table feeds **h2m_variant** or a
library design step unchanged. This skill does **not** do the H2M rename —
that adapter belongs to `h2m_variant`, and two copies would drift.

Useful flags: `--classifications LOF` for the nonsense/frameshift/splice
preset, `--protein-change R882H` (the `p.` prefix is added and stripped both
ways), `--per sample` to reproduce the notebooks' sample-level counting,
`--fallback genie` to check a variant absent from the disease cohort recurs
pan-cancer, `--cosmic-role` for the oncogene/TSG annotation.

### Co-occurrence and mutual exclusivity

```bash
python scripts/cohort_analysis.py comut --cohort genie \
    --disease 'B-Lymphoblastic Leukemia/Lymphoma' \
    --focal TP53 --vs top50 -o OUT
```

Builds the patient × gene binary matrix, runs Fisher's exact for every focal ×
comparator pair, and reports the full `a/b/c/d` contingency table alongside
`odds_ratio`, `p`, `log_or` and `or_floored`. `--matrix` skips the cohort load
and takes a binary matrix you already have.

```
OUT/cooccurrence.csv                 a,b,c,d, odds_ratio, p, q_within_focal, q_global, or_floored
OUT/cooccurrence_volcano_frame.csv   x_log_or, y_neglog10, size_n_double, correction_used
OUT/patient_gene_binary.csv          the matrix, for reuse
```

**Two q columns, always both.** `q_within_focal` is BH within each focal gene
— one family per volcano panel, 20 tests — which is what `CH.ipynb` does.
`q_global` is BH across the whole grid, which is what the original
`Figure_2.R` (line 451) does. On the same data these differ, and the
within-focal family is systematically the more permissive. Pick one with
`--correction`, and print it into the figure; a reviewer comparing your q
values to the R code's will find the discrepancy otherwise.

## Outcome analysis

```bash
python scripts/cohort_analysis.py surv \
    --table /path/outside/repo/final_data_matrix.RData --subset de_novo \
    --focal-group ch --vs top20 --contrast four_arm --min-per-arm 10 -o OUT
```

Required columns: `Sample`, `Gene`, `Time_to_OS` (**days**), `Censor` (1 =
death observed). Optional and used: `Subset`, `variant_type` (for the FLT3
ITD/TKD split), `amino_acid_change`, and anything you name in `--covariates`.
`--focal-group ch` uses the 10-gene grouped CH list from `CH.ipynb`.

Every run prints a **cohort QC block first** — patients, events vs censored,
median follow-up, `Censor` levels, and a verdict on whether `Time_to_OS` looks
like days or years. Read it before trusting any HR.

```
OUT/pair_survival_summary.csv  HR, 95% CI, p, q_bonf, q_fdr, contrast, n per arm, fit_warnings
OUT/pairs_skipped.csv          every untested pair and WHY — missing ≠ non-significant
OUT/km_frame_<pair>.csv        per-arm times and events, for sci-plot
OUT/logrank_<pair>.csv         pairwise log-rank p and stars between arms
```

### Name the contrast, every time

`--contrast` picks the model, and the three answer **different questions and
give different HRs on identical data**:

| `--contrast` | who is compared |
| --- | --- |
| `double_vs_driver_only` *(default)* | among comparator-mutant patients only, does adding the focal mutation change outcome. Wild-type and focal-only patients are **dropped**. This is `CH.ipynb`. |
| `double_vs_rest` | double-mutant vs everyone else in the cohort. This is `Figure_2.R`'s `STRATA_all == 6 vs 0`. |
| `four_arm` | wild-type / focal-only / comparator-only / double as four arms — KM only, no single HR. The published Figure 2 curves are these. |

On a 600-patient test cohort the same pair gave HR 1.62 under
`double_vs_driver_only` and 2.00 under `double_vs_rest`. The contrast name is
written into every result row and into `provenance.json`; put it in every KM
legend too.

**"Redo the Figure 2 KM" means `--contrast four_arm`.** The notebook's default
2-arm contrast will silently hand back a different plot and a different HR
from the published one.

## Gotchas

**`CANCER_TYPE` and `CANCER_TYPE_DETAILED` are not interchangeable.** Verified
on GENIE v17: `CANCER_TYPE == 'B-Lymphoblastic Leukemia/Lymphoma'` gives 1459
samples, an exact `CANCER_TYPE_DETAILED` match gives 886, because the detailed
field splits B-ALL into ~10 subtypes (`...BCR-ABL1`, `...Hyperdiploidy`,
`...BCR-ABL1 Like`). The broad column is usually what you want for library
design. A README in `PhD-FSR-MH-Lab` claims the opposite of what its own code
does — anyone "fixing" the code to match it loses 573 samples. The skill prints
which column it used and records it in `provenance.json`.

**Address clinical columns by name, never by index.** GENIE's
`data_clinical_sample.txt` starts `PATIENT_ID, SAMPLE_ID, ...` with
`CANCER_TYPE` at index 6; MSK-HEME's starts `SAMPLE_ID, PATIENT_ID, ...` with
`CANCER_TYPE` at index 9 (both verified). Index-based code written against
GENIE returns the wrong rows on MSK-HEME, silently.

**A count without a denominator is not a frequency.** FLT3
`p.D593_F594insNS` is literally 0 patients in MSK-HEME AML (296 samples) and 68
in pan-cancer GENIE (~200k). `n_cohort`, `freq` and the 95% CI are columns, not
prose, so they survive into a figure. The pan-cancer fallback is a **recurrence
check, not a frequency source** — selection belongs on the disease cohort.

**Patient vs sample counting changes the number.** GENIE ids are
`GENIE-[CENTER]-[PATIENT]-[SAMPLE]` and one patient commonly has a primary and
a relapse sample. Default here is patient-level (`n_patient`); `--per sample`
reproduces the notebooks' sample-level `nunique`. The `counted_per` column says
which.

**`fisher_for_pair` in `CH.ipynb` fabricates an effect size — real bug.** When
the odds ratio is exactly 0 it substitutes the literal `0.06982398`, and NaN
becomes 1.0, purely so `log(OR)` plots without an infinity. On the real de novo
MetaAML cohort that fires for **15 of 200** focal × comparator cells, which then
land at `log(OR) = -2.66` on every volcano, bubble and heatmap as if that were
measured — when it actually means "zero double-mutant patients". The p-value is
unaffected, so such a cell reads as a strong, *significant* mutual exclusivity
at an invented effect size. **What to do:** the engine here uses `safe_log_or`
and flags those cells in `or_floored`; shape or annotate them differently, and
never quote their `log_or` as an effect size. `comut.fisher_for_pair` is kept
only to reproduce published figures byte-for-byte, and says so in its docstring.

**`Time_to_OS` is in days and nothing in the name says so.** Max observed in
MetaAML is 6780. A Cox model fits happily on days and returns an HR per day
that looks plausible and compares to nothing published. `Time_to_OS_years` is
built by `/365.0`; the QC block gives a verdict, and a max under ~30 means
someone already converted.

**Skipped pairs are missing, not null.** The per-arm minimum (`--min-per-arm`,
default 5) drops pairs entirely — 93 of 200 in the original run. A later
`.pivot` then yields NaN cells that read as null results. `pairs_skipped.csv`
names every one and why. Note `CH.ipynb` used 5 per arm while `Figure_2.R` used
`n_pts_both >= 10`, so the two disagree on which pairs exist at all.

**The Cox model is univariate unless you say otherwise.** No adjustment for
age, sex, cytogenetic risk or treatment, all of which sit unused in the table.
For AML outcome this is not a safe default: NPM1 and FLT3-ITD are both among
the usual comparators *and* the dominant prognostic confounders. Pass
`--covariates Age,Sex,Cytogenetics` and report the adjusted HR beside the
univariate one. Convergence warnings are captured per fit into
`fit_warnings` rather than suppressed — the notebook runs under a blanket
`warnings.filterwarnings('ignore')` twice over, so a barely-converged fit there
reports an HR and CI with no visible complaint.

**Three different significance thresholds in the source.** `survival_status`
uses raw `p <= 0.05`; the original's `label` and forest/KM selection used
`q_value_bonf <= 0.15` (unusually loose, and its own notes say it was used
"because that was explicitly requested"); its heatmap masked on raw `p > 0.05`.
All three are emitted here as separate columns. Never report "significant"
without naming the threshold, the correction and the contrast.

**`normalize_flt3` is not exhaustive.** It maps `variant_type` in
`{ITD, INDEL}` → `FLT3-ITD` and `{SNV, Deletion, other}` → `FLT3-TKD`, but
`Unknown` keeps the bare symbol `FLT3`. On the de novo subset that is 0 samples
so it is invisible today; on another `Subset` or a refreshed cohort you get
three competing FLT3 entities and the unsplit one can displace a real driver
out of a top-N list. The loader warns when any FLT3 row survives the split.

**Point mutations and indels only.** Fusions and rearrangements live in
`data_sv.txt`, a wholly separate table, and copy number in the CNA files;
neither is covered. So `BCR-ABL1` and `ETV6-RUNX1` queries return zero rows,
and "which genes co-mutate with TP53" silently means "by point mutation only" —
which makes a negative co-occurrence result for a CNA-driven gene like
`CDKN2A` or `IKZF1` uninterpretable from this table alone. FLT3-ITD partly
lives in the MAF as `In_Frame_Ins`, so an ITD query half-works, which is worse
than failing: you get a plausible undercount.

**The Cancer Gene Census is not provisioned yet.**
`$LAB_DATA_ROOT/PublicDatabase/COSMIC/` currently holds only
`Cosmic_MutantCensus_Tsv_v99_GRCh37` — an unrelated 500 MB table with
`UPPER_SNAKE` columns and no `Role in Cancer`. `--cosmic-role` therefore fails
with download instructions rather than a missing-path traceback: get
`Census_all.csv` from cancer.sanger.ac.uk (free academic login) and save it
there, or pass `--census`. The old notebook read it from a timestamped filename
inside an unrelated repo; that path is deliberately not carried over.

**`HGVSp_Short` matching is exact string matching.** `--protein-change R882H`
is expanded to try `p.R882H` too, but frameshifts carry their new stop
(`p.G646Wfs*12`), so a truncated `p.G646Wfs` matches nothing. Position parsing
takes the first integer, so an insertion like `p.Q577_M578insSPE` is placed at
residue 577 — its left anchor only.

**`--snp-only` exists and is not the default.** The source notebook restricted
GENIE to `Variant_Type == 'SNP'` before any co-occurrence test, dropping every
indel including the clinically important ones. The flag reproduces that and
warns.

**GENIE v17 is GRCh37.** `NCBI_Build` is uniformly GRCh37 and the local COSMIC
is `GRCh37`, so anything downstream must pass `ver=37` and pair with
`GCF_000001405.25_GRCh37.p13` / `gencode_v19_GRCh37.db`. Handing GRCh37
coordinates to a GRCh38 annotation produces plausible-looking wrong variants,
not an error.

**GENIE's MAF is ~3 GB.** `--genes` filters while reading in chunks, which is
what keeps it off the heap; without it the whole table is read as the notebooks
do, costing several GB of RAM and tens of seconds to minutes. Only `v17.0` is
exercised — the other releases are on disk and resolve, but warn.

## Environment

`freq` and `comut` need only pandas, numpy, scipy, statsmodels. `surv`
additionally needs `lifelines`, and `pyreadr` for an `.RData` cohort.

Use the Python 3.12 conda base (`/opt/miniconda3/bin/python` — verified:
lifelines 0.30.0, pyreadr 0.5.6, statsmodels 0.14.4), **not** the system
`/usr/bin/python3`, which has neither. Keep this out of the `pegg` envs: those
are pinned to Python 3.9/3.10 with `scikit-learn==1.1.1` and `numpy<2` for the
pickled PRIDICT2 models, and this skill is numpy-2/pandas-2 clean. Don't try to
make one env serve both.

The scripts are written to 3.9 syntax anyway, so they import in either: no
match statements, no nested-quote f-strings (`f'{os.environ["X"]}/y'` is a
`SyntaxError` before 3.11 — paths here use plain concatenation and
`os.path.join`).

## Verifying

Engine bodies are lifted verbatim from the notebooks named at the top of each
script, so the two cannot drift. Anchors from the existing
`MetaAML_results/CH/` outputs, worth asserting after a change:

- `DNMT3A`/`NPM1` co-occurrence: OR 5.3769112, p 5.3157e-60, n_double 345
- `DNMT3A + NPM1` survival: HR 1.2275637 (0.9855174–1.5290572), p 0.0672841,
  n_double 345, n_driver_only 362
- 107 of 200 pairs survive the n≥5 guard; 15 of 200 co-occurrence cells are
  `or_floored`; `len(CH_UNIQUE) == 10`

Because the engine deliberately **changes** one thing — `safe_log_or` replacing
the `0.06982398` literal — the right regression check is that `p` and
`n_double` reproduce exactly while `odds_ratio` differs **only** on the 15
`or_floored` cells.

```{note}
`CH.ipynb`'s markdown says "9 unique CH genes, so the volcano panel is a full
3x3 layout". The code's list has **10** and every grid is `subplots(3, 4)`.
Anything rebuilt from that prose makes the wrong grid. `survival.py` asserts
the count. Relatedly, cell 12 there has a real bug —
`for ax in axes[len(CH_UNIQUE):]` slices *rows* of a 2-D array and iterates
zero times, so two empty panels keep their spines in `co-occurrence.pdf`; use
`axes.ravel()[n:]`, which the neighbouring cells do correctly.
```

## Related skills

**sci-plot** owns every figure's appearance — load it for any plot built from
the frames here. **h2m_variant** owns human↔mouse conversion and the MAF→H2M
rename; hand it `variant_frequency.csv`. Guide and pegRNA design, and the
designability verdict, belong to the library-design and `pegrna_oligo` skills:
"most frequent" and "designable" are different queries, and this skill answers
the first honestly rather than pre-filtering for the second.
