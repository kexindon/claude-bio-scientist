---
name: library-design
description: >-
  Design prime-editing pegRNAs or base-editing sgRNAs from a variant table
  with PEGG — PAM sites, RTT/PBS lengths, ABE vs CBE, sensor cassettes, silent
  bystanders, synonymous-codon and safe-targeting controls, efficiency scores,
  and a per-variant dropout ledger. Use for "design guides for these
  mutations", "make a prime editing library", "ABE or CBE for this variant",
  "add a sensor", "add silent bystanders", "score/rank my pegRNAs". Ends at
  the design table; for an IDT Golden Gate order sheet use pegrna-oligo.
---

# library-design — variant table → pegRNA/sgRNA design table

## Where this skill stops

**This skill ends at the design table**: one row per candidate guide, with
sequences (`Protospacer`, `RTT`, `PBS`, `RTT_PBS`), geometry
(`Distance_to_nick`, `RHA_size`, `PAM_disrupted`), sensors, bystander columns,
and efficiency scores. It also emits the **pooled-synthesis** oligo column
(`pegRNA_oligo` / `gRNA_oligo`) — *one long sequence per guide*, for chip or
array synthesis of a whole pool.

| you want | skill |
|---|---|
| guides designed from a mutation list | **this skill** |
| one pooled-synthesis sequence per pegRNA | **this skill** |
| an **IDT order sheet** of individually ordered Golden Gate oligos — three annealed duplex pairs per pegRNA, `.dna` maps, cloning | **`pegrna-oligo`** |
| human variant → mouse orthologue modelling | `h2m-variant` |
| reading sensor amplicons out of screen FASTQs | `ngs-screen` |
| plot styling for any figure made here | `sci-plot` |

The discriminator against `pegrna-oligo` is the **synthesis route**, not the
file format. Pooled/array synthesis = one sequence per pegRNA, ordered as a
pool, which is what this skill emits. Golden Gate = protospacer + scaffold +
3'-extension oligos ordered individually and annealed, which is `pegrna-oligo`.
Hand it `*_designs.csv.gz`.

## Environment — check it first, every time

**Python 3.9 or 3.10, `scikit-learn==1.1.1`, `numpy<2`.** These are not soft
pins. pegg's scorers (`V3_model_full.pickle`, `V3_model_nopos.pickle`,
`RF_PE.pkl`) are pickles built against sklearn 1.1.1. On a newer sklearn they
raise an **unpickling error that does not explain itself** — nothing in the
traceback says "wrong sklearn" — or silently mis-score. That takes out
`ontarget_score`, and therefore `PEGG2_Score` and `RF_Score`, hours into a run.

```bash
conda env create -f envs/pegg_env.yml && conda activate pegg_env
python scripts/check_env.py            # exits non-zero on a wrong pin
python scripts/check_env.py --optiprime   # also require a working OptiPrime
```

Run `check_env.py` **before designing**, not after a failure. It also verifies
`pegg>=3.1.0` by probing `bystander.reverse_frame_anchor` for its
`genomic_span` argument — see the indel bug below.

OptiPrime scoring lives in a **separate Python 3.11 environment** (its jax
stack needs numpy>=2, so it cannot share `pegg_env`); pegg shells out to it via
`OPTIPRIME_SRC` and `OPTIPRIME_PYTHON`.

## Inputs

**A variant table in MAF column names** (pegg's `input_format='cBioPortal'` —
the only format that gives genomic coordinates, splice checks and a per-row
reading frame, so the only one to use for a real library):

`Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Variant_Type`
(`SNP`/`DNP`/`TNP`/`ONP`/`INS`/`DEL`/`INDEL`), `Reference_Allele`,
`Tumor_Seq_Allele2`.

`Chromosome` must be **int for 1–22 and the strings `'X'`/`'Y'`, no `chr`
prefix** — that is how `chrom_dict` is keyed. H2M output needs renaming
(`gene_name_m`→`Hugo_Symbol`, `start_m`→`Start_Position`, …);
`prep_input.py --from-h2m m` does it.

The other two formats (`'WT_ALT'` with a `WT`/`ALT` column pair,
`'PrimeDesign'` with bracket notation in `SEQ`) carry no coordinates: they need
`ORF_start` for bystanders, cannot be splice-checked, and **cannot be audited**
(see Gotchas).

**A reference genome and, for silent bystanders, a GENCODE `.db`** — both
reached through `$LAB_DATA_ROOT`, both the same build as the variant
coordinates (GRCh37↔v19, GRCh38↔v44, GRCm39↔vm33). Never hardcode a path; see
`docs/data-layout.md`.

## The procedure

### 1. Prepare and frame-annotate the input

```bash
python scripts/prep_input.py MUTATION_DESIGN.csv -o design_input.csv \
    --build GRCm39 --from-h2m m --tx-column tx_id_m
```

Order matters and `prep_input.prepare()` enforces it: frame annotation →
chromosome normalization → assert every chromosome is in `chrom_dict` →
`reset_index`. The assert is the point — without it a `'chr7'` becomes a
`KeyError` deep in the per-mutation design loop, minutes after the genome load
has already been paid for. `reset_index` matters because pegg's `mutation_idx`
*is* the input row index.

`--tx-column tx_id_m` takes the frame from the transcript H2M already chose, so
the design matches the orthologue that was actually modelled. **Always pass it
explicitly** — pegg's default is `tx_id_h` (human).

Read the printed `cds_annotation_report` before continuing. It catches a CDS
that is not a whole number of codons, overlapping blocks, and variants outside
their own CDS — each of which otherwise just yields zero bystanders with no
error to explain why.

### 2. Design

```bash
# prime editing, with silent bystanders and OptiPrime scoring
python scripts/design.py design_input.csv -o out/ --name ch_aml \
    --mode prime --build GRCm39 --tx-column tx_id_m \
    --silent-bystander --optiprime --optiprime-jobs 7

# base editing (MBESv2 profile: NG PAM, ABE+CBE, window 4-8)
python scripts/design.py design_input.csv -o out/ --name mbes \
    --mode base --build GRCm39
```

The profiles in `design.py` are the **validated** ones, lifted with their
justifications from the CH/AML and MBESv2 notebooks. For prime: `RTT_lengths =
range(10,31)`, `PBS_lengths = range(5,17)`, `min_RHA_size=8`,
`RE_sites=['CGTCTC']`, `sensor_length=55`, `seed=0`. Those PBS/RTT ranges match
OptiPrime's own design rule, and `min_RHA_size=8` is its homology-arm floor
(pegg's default of 1 allows arms far shorter than OptiPrime would score).

Writes, in this order — the order is load-bearing:

```
out/ch_aml_designs.csv.gz          every design, gzipped
out/ch_aml_oligos.csv              pegRNA_id + pooled-synthesis oligo
out/ch_aml_dropout_ledger.csv      one row per INPUT variant, naming what killed it
```

### 3. Audit the bystanders (whenever `--silent-bystander` was used)

```bash
python scripts/audit_bystanders.py out/ch_aml_designs.csv.gz design_input.csv \
    -o sanity_check/bystander_audit.csv --build GRCm39 --tx-column tx_id_m
```

Splices each RTT back onto the CDS, re-translates, and compares each bystander
design against the protein the **plain** designs for that same mutation agree
on. That baseline is the whole trick: "at most one amino acid changed vs
reference" only holds for missense, and frameshift/nonsense/in-frame-indel
variants legitimately rewrite many residues.

### 4. Matched synonymous controls (optional)

```bash
python scripts/silent_controls.py design_input.csv -o silent_input.csv \
    --build GRCm39 --seed 0
python scripts/design.py silent_input.csv -o out/ --name ch_aml_controls \
    --mode prime --build GRCm39 --tx-column tx_id_m \
    --silent-bystander --optiprime --optiprime-jobs 7
```

One synonymous-codon twin per variant, at the same codon. Design it with
**identical parameters**, bystanders included, so the controls span the same
efficiency range as the real designs. Variants with no comparable control
(indels, substitutions spanning >1 codon) are skipped and listed in
`*_skipped.csv`.

For arm-level controls instead, pegg has `library.safe_muts`,
`library.aavs1_muts`, `library.nontargeting_guides` and `library.library_maker`
— one-line calls, so use them directly rather than through a script here. Read
the Gotchas on all three first.

## Runtime

Design itself is **minutes**: genome load ~20 s and a few GB of RAM, then
roughly a second or two per variant. Azimuth is called once per unique
`Protospacer_30`, not per pegRNA, so `PEGG2_Score`/`RF_Score` are nearly free.

**OptiPrime is the only expensive step**: ~0.2155 s per pegRNA per worker, plus
~20 s of JAX JIT per worker, and it will not shard below 500 rows. The CH/AML
library (665k pegRNAs, bystanders on, `pegRNAs_per_mut='All'`,
`optiprime_prefilter=None`) took **~6.5 h on 7 workers**; the control library
another **~4.4 h**. Near-linear in cores, so set `--optiprime-jobs` to
cores − 1.

Memory and disk are the usual wall, not CPU: 665k rows × ~70 columns is ~1.45
GB in pandas, and the two libraries together are 2.9 GB of plain CSV. Always
write gzipped (~2–8% of raw) — and note the uncompressed row count is far past
Excel's 1,048,576-row limit, so these tables are for programmatic filtering.

Everything above ran on an 8-core macOS laptop. No cluster is needed; a cluster
helps only by giving OptiPrime more cores.

## Gotchas

**The PEGG 3.1 silent-bystander indel bug — regenerate affected libraries.**
Under pegg 3.0/3.0.1, `bystander.reverse_frame_anchor` had no `genomic_span`
argument. On an INS/DEL/INDEL whose RTT is read on the opposite strand, the
frame anchor was shifted by `(ref_len - alt_len) % 3`, and the shift is
*self-consistent*: the bystander options it produced passed pegg's own internal
translation check while being **non-synonymous in the real reading frame**. No
error was raised. If you have a library with `silent_bystander=True` on any
indel, generated before 3.1, **it must be regenerated** — and run
`audit_bystanders.py` on the result. `check_env.py` probes for this
specifically.

**`prime.other_filtration` never scans the protospacer for poly-T in prime
mode.** The source reads `if terminator_sequence in (pbs_rtt or proto)`
(prime.py:1456, again at 1555). Python evaluates `pbs_rtt or proto` to
`pbs_rtt`, because a non-empty string is truthy — so only the extension is ever
checked. This is a real bug in pegg: `contains_polyT_terminator == False` does
**not** mean the protospacer is terminator-free. `design.py:post_filter`
re-checks `Protospacer` itself. Base mode is unaffected (it has only a
protospacer).

**pegg flags restriction sites and terminators but does not drop those rows.**
`contains_RE_site` and `contains_polyT_terminator` are columns, not filters —
dropping them is the caller's job. `RE_sites` also accepts only A/T/C/G
(degenerate codes assert), and `other_filtration` adds each site's reverse
complement itself, so pass `['CGTCTC']`, not both strands.

**Checkpoint *after* `prime_oligo_generator()`, never before.** The oligo
generator is what adds `pegRNA_oligo`; a checkpoint written before it lacks
that column, and resuming from it dies with `KeyError: ['pegRNA_oligo'] not in
index`. `design.py` does oligo generation then the gzipped write with nothing
between them, because a `NameError` in a later notebook cell once destroyed a
completed 6.5 h scoring run.

**`has_silent_bystander` round-trips through CSV as the wrong dtype.** pegg's
not-designable branch pads every column with `None`, pushing the column to
`object`; `run()` casts it back but `pd.read_csv` does not. `~df['has_silent_
bystander']` then does *integer* negation instead of logical not, silently
selecting everything. Re-cast after loading any checkpoint —
`design.recast_bystander_flag()`.

**Bystander designs are chosen at random — always pass `seed`.** pegg prints a
note when `seed=None`, which is easy to miss in a long log. The scripts default
to `seed=0`.

**The bystander design space explodes with `bystander_window_nt`, not with
`max_bystander_muts`.** Measured: lifting `max_bystander_muts` from 2 to
unlimited adds ~4% (3,784→3,928 pegRNAs/variant), because within ±5 nt there
are only a couple of codons to work with. The window is what scales — per
variant: 5→3,928 / 10→19,440 / 15→87,984 / 30→1,334,448. At window=30 that is
~126 M pegRNAs and roughly **20 days** of OptiPrime scoring on 7 cores. Leave
it at 5 unless you have decided to pay for it.

**`optiprime_cutoff` drops whole mutations, not just weak designs.** Measured
on the CH/AML library: a 0.20 floor kept 23/40 mutations and lost Npm1,
Idh1/2, Tet2 and Trp53 outright; 0.15 kept 30/40; 0.10 kept 33/40. Those
variants are genuinely hard (Ppm1d's best design scores 0.006, Jak2 0.053).
Design with `optiprime_cutoff=None`, leave the score in the table, and
threshold afterwards with the coverage cost visible.

**`optiprime_prefilter` defaults to 4, which quietly reinstates `PEGG2_Score`
as the gatekeeper** — it OptiPrime-scores only the top `pegRNAs_per_mut × 4`
per mutation by PEGG2. Measured here, only 10 of 22 mutations kept an identical
top-50, and full scoring recovered a mutation the prefilter had lost. Set
`optiprime_prefilter=None` for anything that will be ordered.

**An unmatched `optiprime_group` is not an error.** `rx_model.load_rates`
leaves `group_factor=0` and logs only at INFO, so a fabricated `BOGUS_CELL`
returns a plausible 0.313. Cell type is how MMR status enters the model, so
that is a silent wrong answer. `design.py` validates against pegg's
`VALID_GROUPS` and hard-fails.

**Keep `optiprime_group='Liu_HeLa'` even when working in another cell type.**
The group dominates the score (0.174 HeLa vs 0.407 HEK293T on the same pegRNA).
HeLa is MMR-proficient and HEK293T partially MMR-deficient, so pegRNAs
optimized on HeLa parameters transfer to HEK293T but not the reverse — HeLa is
the conservative choice.

**A Jupyter kernel does not inherit the shell environment.** `OPTIPRIME_SRC` /
`OPTIPRIME_PYTHON` exported in `.zshrc` are invisible to a running kernel, and
`prime.run(optiprime=True)` then raises `OptiPrimeError` despite a working
install. Set them with `os.environ[...]` in the notebook, or pass
`optiprime_src=` / `optiprime_python=`.

**`add_cds_to_variants` defaults to `tx_column='tx_id_h'`** (bystander.py:1340)
— *human*. On a mouse-only table that is a `KeyError`; on an H2M table carrying
both `tx_id_h` and `tx_id_m` it silently annotates the **human** reading frame
and designs bystanders in it. Always pass `tx_column` explicitly.

**`GENE_CDS` is keyed by transcript id, not gene name**, when it came from
`add_cds_to_variants(tx_column=...)`. Grouping an audit by `Hugo_Symbol` looks
up a key that does not exist and skips every gene without saying so.

**`seq_start` exists only on the `'cBioPortal'` path** (set by
`df_formatter`). Anything that maps a design back to genomic coordinates — the
bystander audit included — works only for that input format. A library designed
from `'WT_ALT'` or `'PrimeDesign'` input cannot be audited;
`audit_bystanders.py` says so by name rather than guessing.

**`run_base` silently drops every non-SNP** with `auto_SNP_filter=True` (the
default), printing one line. Indels never reach base-editing design. Check the
input count against `mutation_idx.nunique()` — or read the dropout ledger,
which names this.

**A non-NGG PAM warns in prime mode and is silent in base mode.** prime.py:1831
prints `Warning: Efficiency predictions on non-NGG PAMs are un-tested`;
`base.py` has **no such warning at all** (only a dead comment at line 390). The
validated base profile here is `PAM='NG'` — the one configuration you will
actually run is the one with no warning. Treat `OnTarget_Azimuth_Score` and
`PEGG2_Score` as unvalidated on any non-NGG design, in either mode.

**`PBS_lengths` above 17 assert; `RTT_lengths` above 30 only warn** — and a
long RTT outruns the sensor: at `sensor_length=55`, `sensor_wt` comes back
`None` and `other_filtration` *raises* on it. Cap RTT at 30 and filter
`sensor_error` before anything that reads `sensor_wt`.

**`run()` scores and caps before the poly-T and sensor filters run**, so a
tight `pegRNAs_per_mut` can leave a variant with zero rows after filtering even
though `run()` returned 50. Check per-variant coverage after filtering — the
dropout ledger does.

**`pegRNAs_per_mut` does not bound the total when bystanders are on** unless
`cap_total_per_mut=True` (now the default). And do not cap by overall rank:
plain and bystander designs are not on one score scale (a bystander can disrupt
the PAM and shift `RTT_GC_content`, both positively weighted), so bystanders
sweep the top and a variant whose only viable designs are plain loses every
row. The budget is split between design types for exactly this reason.

**`library.safe_muts`, `library.aavs1_muts` and `library.library_maker` only
work on GRCh37 and GRCm38.** The safe regions and AAVS1 coordinates are
hardcoded per build. A library designed against GRCm39 needs a separate GRCm38
`chrom_dict` for its controls, or the controls regenerated in the target build.

**`library.nontargeting_guides` is human-only and silently clamps to 1000.**
There is no `organism` argument (library.py:403 reads
`non_targeting_human_prime.csv` unconditionally) and >1000 is clamped with a
print. The mouse source ships in the package as
`non_targeting_unformatted_mus.txt` but is **unprocessed and not wired up** —
every library this skill is built around is mouse, so a mouse non-targeting arm
has to be processed from that file by hand. Say so rather than shipping human
guides into a mouse pool.

**`prime.goldengate_oligos()` is a stub that returns `None`.** The Golden Gate
path is the `pegrna-oligo` skill, not pegg.

**`PEGG2_Score`'s weight dict in the source is dead code.** `peggscore2`
defines `factor_weights` with round numbers and immediately overwrites them
with the fitted vector `fw`. Do not read the dict as the model.

**`prime.genome_loader` returns a 2-tuple** `(chrom_dict,
filtered_descriptions)` and drops alternate, unplaced, patch and mitochondrial
records. `h2m.genome_loader` returns `(records, index_list)` — a different
shape, so the two are not interchangeable. Use
`prep_input.chrom_dict_from_h2m()` to convert.

## Verification

Before ordering anything:

1. **`check_env.py` exits 0.** Especially the `pegg>=3.1.0` probe.
2. **Read the dropout ledger.** Every intended variant should have `n_kept > 0`;
   anything with a `dropped_at` reason is absent from the library, not merely
   low-ranked.
3. **`audit_bystanders.py` reports zero non-silent designs**, and
   `plain_disagree` is 0 for every mutation — a non-zero `plain_disagree` means
   the splice is wrong for some geometry and the baseline itself is suspect.
4. **Silent controls**: `silent_controls.py` asserts every codon swap is
   synonymous and every `Reference_Allele` matches the genome. Both assertions
   must pass; a REF mismatch means pegg designed against the wrong sequence.
5. **Spot-check orientation** on any new data, the same check `pegrna-oligo`
   relies on: `protospacer.find(rc(PBS)) >= 0` and `protospacer.find(PBS) < 0`.

Figures made from these tables (score distributions, guides-per-mutation,
coverage-vs-cutoff curves) **defer to the `sci-plot` skill** for fonts, colors,
sizes and PDF export. Keep any plotting code here to the mechanics — what goes
on which axis — and let `sci-plot` own the style. A quick QC glance can be
emitted directly; anything headed for a paper hands off.

## Scripts

| script | what it is |
|---|---|
| `check_env.py` | pin and env-var guard; the `pegg>=3.1.0` indel probe. Run first. |
| `prep_input.py` | H2M→MAF rename, chromosome normalization, chrom_dict presence assert, CDS/frame annotation + report. Lifted from the CH/AML notebook cells 6, 11, 13. |
| `design.py` | CLI around `prime.run` / `base.run_base` with the validated profiles, the post-filters (including the protospacer poly-T re-check pegg misses), oligo-then-checkpoint order, and the dropout ledger. Profiles lifted from CH/AML cell 16 / cell 30 and MBESv2 cell 4. |
| `audit_bystanders.py` | per-transcript CDS cache + per-base RTT splice + re-translation, compared against the plain designs' consensus protein. Lifted verbatim from CH/AML cell 26. |
| `silent_controls.py` | one synonymous-codon twin per variant, with the synonymous and REF-matches-genome assertions. Lifted verbatim from CH/AML cell 29. |

Each script names its source notebook and cell at the top. pegg itself is
always called through — never vendored, never reimplemented — so a pegg upgrade
reaches this skill without edits.
