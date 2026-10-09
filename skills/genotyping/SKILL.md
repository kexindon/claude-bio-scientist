---
name: genotyping
description: Call genome-editing outcomes from Sanger .ab1 traces — editing efficiency, bystander edits and homozygous/heterozygous/WT clone calls at a codon plus nucleotide and amino-acid sequence logos; also turns FlowJo gate-statistic .xls exports into population frequencies, absolute counts and fold expansion. Use for "did the edit work", "editing efficiency", "quantify editing from the ab1", "check for bystander editing", "what fraction is R248Q", "which clones are homozygous", "genotype the colonies", "my sequencing came back", "knock-in/HDR efficiency", "what percent LSK", "FACS", "fold expansion" (编辑效率, 测序结果, 纯合还是杂合, 流式). Not for FASTQ — amplicon NGS editing quantification is ngs_screen.
---

# genotyping — arrayed validation: one .ab1 per sample, one gate table per run

Two bench readouts that answer "did it work", on a per-sample basis:

- **Sanger** — one `.ab1` per sample, quantified at a codon.
- **Flow** — a FlowJo "export group statistics" `.xls`, turned into frequencies
  and absolute counts.

## When not to use this

| You have | Use |
|---|---|
| Pooled screen FASTQ, amplicon NGS, CRISPResso2 | **`ngs_screen`** |
| An exact allele frequency you will publish | **`ngs_screen`** — see *Peak heights are not allele fractions* |
| A design to make, not validate | `library_design`, `pegrna_oligo` |
| A paper figure, not a QC glance | **`sci-plot`** (this skill defers all styling) |
| Raw `.fcs`, no FlowJo workspace | Not supported. `pip install readfcs` or re-gate in FlowJo and export statistics |

**Arrayed, single-locus, one file per sample.** The moment the readout is a
pooled library or a FASTQ, it is `ngs_screen` — including when the question
sounds identical ("what's my editing efficiency at R248Q"). The two differ by
input, not by question.

## Step 0 — environment

`logomaker` is the constraint. It **is** importable from the system Python 3.9
(`/usr/bin/python3`, logomaker 0.8.7), which is what `python3` resolves to on
this machine and why the source notebook runs as written. It is **absent from
every conda env** — `base`, `myenv`, `pegg`, `pegg3`, `pegg_env`, `stat_env`.
So the notebook breaks the moment anyone activates an env or moves to the
server.

```bash
pip install 'logomaker>=0.8,<1.0'      # pin pre-1.0
```

Also needs `biopython>=1.81` (verified 1.83), `pandas`, `numpy`, `matplotlib`,
and for the flow path `xlrd==2.*` (verified 2.0.2).

This is **not** a `pegg` skill: no `scikit-learn==1.1.1`, no `numpy<2`. Do not
run it in `pegg_env`. Neither path needs `$LAB_DATA_ROOT`; only an optional
reference CDS for `locate_window.py` does, and that reads the env var with the
standard guard.

## Sanger: the procedure

### 1. Find the window — do not guess it

`start`/`end` are **1-based inclusive coordinates in the AB1's own called
sequence (PBAS1)** — what SnapGene or Benchling shows for that read. Not
genomic, not transcript, not amplicon coordinates.

Let the aligner derive them, together with frame and strand:

```bash
python scripts/locate_window.py --ab1 KP_Sort-p53_F.ab1 \
    --cds-file trp53_region.fa --cds-first-residue 230 \
    --residue 245 --flank-codons 2
```

It prints `start`, `end`, `reverse_complement`, `aa_anno_range`, the called
bases, and the coding-strand translation to eyeball — then the exact
`call_edit.py` command to run next. `--cds` is the **in-frame** coding sequence
of the amplified region, so frame comes out derived rather than assumed.

**Always read the printed bases and confirm the WT residues are what you
expect** (R at 248, R at 273, …) before trusting any number downstream.

### 2. Quantify and call

```bash
python scripts/call_edit.py --ab1 KP_Sort-p53_F.ab1 \
    --start 118 --end 138 --aa-first 249 --aa-last 255 --rc \
    --residue 252 --wt L --alt P \
    --control KP_Pre-p53_F.ab1 -o OUT
```

A whole plate, from a variant table — this replaces the notebook's 15
copy-pasted per-variant blocks:

```bash
python scripts/call_edit.py --table plate.csv -o OUT
```

`plate.csv` columns: `ab1, variant, start, end, reverse_complement, aa_first,
aa_last`, optionally `residue, wt_aa, alt_aa, control_ab1`.

```
OUT/
├── calls.csv                       one row per sample: genotype, props, QC verdict
├── <variant>_aa_proportions.csv    indexed by real residue number
├── <variant>_nt_proportions.csv
├── <variant>_nt_logo.pdf
└── <variant>_aa_logo.pdf
```

Genotype classes are `homozygous_edited` (≥0.85), `heterozygous` (0.25–0.75),
`wild_type` (≤0.10) and `ambiguous` — a clone at 0.80 is **not** filed as a
clean homozygote on purpose.

## Flow: the procedure

List first, always — never glob specimens:

```bash
python scripts/flowjo.py 08-May-2026.xls --list
```

```bash
python scripts/flowjo.py 08-May-2026.xls \
    --samples "Sample 1" "Sample 2" "Sample 3" \
    --initial-lsk 24000 --viability 0.452 --volume-ml 3 --conc-per-ml 12e6 \
    -o OUT
```

Writes `flow_summary.csv` plus the three QC panels `flow_bar_plot.pdf`,
`flow_stacked_linear.pdf`, `flow_stacked_log.pdf`. The hemocytometer numbers
are arguments, not literals: `final total live = viability × volume × conc`.

The gate tree in `summarize_lsk` is the LSK panel as exported by the
08-May-2026 workspace. A different panel calls `gate()` directly rather than
bending that function.

## Figures

**All styling belongs to `sci-plot`.** `flowjo.py` sets up through
`kdplot.setup()`/`save()` and decides only what goes on which axis; it defines
no fonts, sizes or palettes. The segment colors (`#762a83` LSK, `#4dac26` Lin−,
`tab:grey` Lin+) are carried from the source notebook as defaults so existing
panels still match, and take an override. The logo PDFs come from logomaker's
own `classic`/`chemistry` schemes.

These three panels are QC you glance at. **For a paper figure, hand off to
`sci-plot`** — it wants SD caps and replicate dots, which the source notebook's
`.mean()`-only code never had. `absolute_counts()` therefore returns
`n_samples`, `sd_lsk_frac` and `per_sample_LSK_counts` so the dots can be drawn.

## Gotchas

**Peak heights are not allele fractions.** This is the big one. Trace
proportions are peak-height ratios at the called peak, biased by dye mobility,
signal saturation, nearby indels and short read length. They are **only**
meaningful as a difference from a matched WT/unedited control trace read at the
same window and frame. Without `--control`, every row is tagged
`confidence=uncalibrated` and the run prints a warning — honor it. Never quote
an exact allele frequency from Sanger; that needs amplicon NGS (`ngs_screen`).

**An indel makes everything downstream uninterpretable.** In a non-clonal pool
an indel puts the two alleles out of register, so from the indel onward every
position is a superposition of two different bases. The trace does not look
noisy — it looks like a clean heterozygote, at every position. `indel_scan()`
flags a run of ≥8 consecutive positions whose top-base proportion is below
0.75, and `assert_interpretable()` **refuses to report** rather than printing a
confident wrong number. `--allow-indel` downgrades it to a flag for when you
only want the logo as a picture. It is a heuristic and deliberately paranoid: it
does not size or localize the indel.

**A shifted window produces a plausible wrong answer.** The source code does
`range(start-1, end)` with no alignment and no bounds check. Off-by-a-few gives
a clean-looking logo of the wrong codons; off-by-one in frame gives a complete
amino-acid logo that is entirely wrong, and the multiple-of-3 check cannot
detect it. This is why `locate_window.py` exists and why you confirm the WT
residues. (The out-of-range case now raises a readable error instead of the
notebook's opaque `IndexError` from `positions[idx]`.)

**Reverse-complement follows the strand, not the primer name.** For TP53 —
a minus-strand gene — `reverse_complement=True` goes with the **FW** primer
read. In the source notebook the flag is *mostly* the inverse of the primer
name — all 12 `_FW` invocations pass `True`, and 7 of the 8 `_RV` variants pass
`False` — because an FW-primer read is antisense to the CDS. The rule is **not**
absolute: `7_X329_splice_RV-6_RV.ab1` (cells 65/66 and 132/133, window 153–161)
is an `_RV` read passed `reverse_complement=True`. So the filename is a decent
guess and a bad guarantee. Do not reason from it; let `locate_window.py` decide
from sequence — it reports `strand evidence` as anchor counts on each strand.

**The AA collapse assumes codon positions are independent.** It multiplies the
three per-base marginals over all 64 codons. That is false for a real edited
pool — a cell carrying the edit usually carries its bystanders too — so with
more than one edited base in a codon the amino-acid proportions are an
approximation that can show haplotypes which do not exist.
`haplotype_warning()` lists those codons and `calls.csv` carries them in
`multi_mixed_codons`. Say so when reporting such a codon.

**Trace channel mapping is read from `FWO_1`, never assumed.** On these files
`FWO_1` is `b'GATC'`, so `DATA9=G, DATA10=A, DATA11=T, DATA12=C` — **not**
alphabetical. Hardcoding `DATA9..12 → A,C,G,T` silently transposes G/A and T/C.
All four channels are sampled, keyed by base letter, at `PLOC1[idx]` — the
trace index of the called peak, not the base index. `DATA9-12` are the analyzed
(baselined, mobility-corrected) traces; `DATA1-4` are raw and are not what this
reads.

**The source `sanger_weblogo` returns `None` for `logo_type='nc'`** after a
stray `print(nt_df)`, so `aa_df = sanger_weblogo(..., 'nc', ...)` silently binds
`None` in the notebook. The lifted copy here returns `nt_df` from both branches
and drops the print — the only two changes made to it.

**No `makedirs` in the original.** A missing output directory killed a long run
at the very last `savefig`. Fixed in the lift.

**FlowJo statistics are percent-of-PARENT, not percent-of-total.** Anything two
or more gates deep must be chained through every intermediate:
`pct_of_ancestor(9.83, 69.0) → 6.7827`. Skipping the chaining is the single
commonest flow readout error.

**A typo'd gate path returned `(nan, 0)` silently** in the source, turning a
typo into NaN columns and a plot of nothing. `gate()` now raises and prints the
paths sharing that leaf name. FlowJo's quadrant names have exact spacing —
`'Q2: APC-A+ , BV421-A+'` has a space *before* the comma.

**Compensation controls sit in the same export as the samples.**
`Specimen_001_APC.fcs`, `_BV421.fcs`, `_BV510.fcs`, `_PE.fcs`, `_PEcy7.fcs`
alongside `Specimen_001_Sample 1/2/3.fcs`. They have near-zero LSK and would
silently drag down every mean, which is why the source hardcoded three sample
names. Run `--list` and name the real samples.

**Deep gates have tiny denominators.** Sample 1's LSK CD150+ CD34− gate is
**4 cells of 260** LSK. Every `%` is carried with its `#Cells` and the CLI warns
when the deepest gate is under 20 cells. A percentage off single-digit counts is
noise.

**The log stacked bar is a presentational construction, not a log axis.** Each
segment's height is `log10(total) × (count/total)`, so heights are proportional
to composition while the bar's total height reads as `log10(count)`. Segment
boundaries do **not** correspond to their own counts on the 10^n ticks. Keep it
paired with the linear panel and caption it, or a reader will misread segment
sizes.

**`xlrd` 2.x reads only legacy BIFF `.xls`.** An `.xlsx` export raises; the
parser gives a clear message pointing at `openpyxl`.

**Cells 1 and 4 of `flow_plot.ipynb` both write `flow_stacked_log.pdf`** with
opposite stacking order, so whichever ran last won. Cell 4 (LSK at the bottom,
matching the linear panel) is the intended one and is what is lifted; cell 1 is
dropped.

## Verification

`scripts/sanger_logo.py` holds `sanger_weblogo` **verbatim** from
`endogenous_editing.ipynb` cell 2 and is the reference. `scripts/trace.py`
re-expresses the same arithmetic independently. Checked on
`KP_Sort-p53_F.ab1` — the two produce **identical** frames on the forward
strand, the reverse-complemented strand, and the amino-acid collapse.

The flow path reproduces the notebook's own `flow_summary.csv`
**numerically identically**, Sample 1 LSK = `6.78270` from `9.83 × 69.0`.

Worth asserting on new data:

```python
import trace as tr
from sanger_logo import sanger_weblogo
a = tr.base_proportions(ab1, start, end, rc)
b = sanger_weblogo(ab1, start, end, 'nc', rc, None, None)
assert a.round(10).equals(b.round(10)[['A', 'C', 'G', 'T']])
```

```{warning}
**No Sanger regression fixture ships with this skill, and that is not an
oversight.** The `.ab1` files in `hsc-ch-mutation/20260512-HSC-PE7-KP-PK-editing/`
are mouse KP/PK reads (Kras and p53, F and R). All 15 of the notebook's
invocation blocks point at human A549 TP53 reads under
`/Users/kexindong/Desktop/PhD/Experiments/05_TP53_comp_assay/...` that are **not
in the repo** — and its very first block names
`MEC_KP_PK_Genotyping_30-1331125458_ab1/2_L252P_FW_N-P2_FW.ab1`, a file that
does not exist in that folder (a stale path from the human dataset).

So **no known-good start/end exists for any file a user actually has.** The
published windows (118–138 for L252P, 138–152 for R248Q, 57–71 for R273C) are
coordinates in reads this skill cannot see. Derive your window with
`locate_window.py` and verify the WT residues; do not copy a number out of the
notebook. Recording one correct window for one KP/PK read would let a fixture
ship — until then there is nothing to regress against on the Sanger side.
```

## Scripts

| script | what it is |
|---|---|
| `sanger_logo.py` | `sanger_weblogo` **verbatim** from cell 2. The reference. Logos + proportion frames |
| `trace.py` | trace arithmetic re-expressed, **plus** the quality/indel guards that refuse to report |
| `call_edit.py` | the numeric readout: residue calls, control subtraction, genotype classes, plate batching. CLI |
| `locate_window.py` | derives start/end/frame/strand by aligning the read against an expected CDS. CLI |
| `flowjo.py` | FlowJo parser, gate chaining, count math, the three panels. CLI |
