---
name: bio-scientist
description: >-
  Map of the lab's cancer-genomics and genome-editing toolkit — which tool
  answers which question, how they chain, and what data is on the machine.
  Use when the request spans several steps ("go from GENIE frequency to an
  IDT order", "plan the CH pilot"), when it is open-ended ("what can you do
  here?", "where do I start?", "what data do we have?"), when you are new to
  this setup, or when it is unclear which tool applies. Also use to find out
  why something is not possible. Not needed when the request already points
  at one step — ask for that directly and the right tool loads on its own.
---

# bio-scientist — the map

The lab's work runs as one loop: a clinically observed human variant becomes
a mouse model, becomes a guide library, becomes a screen, becomes a hit list
that informs the next design round. Each stage has its own skill. **You do not
have to name them** — asking for the stage loads it. This skill exists for
when the request spans stages, or when it is not yet clear which stage it is.

**If the request names one stage, stop reading and do that instead.** This
file is for routing and planning, not for doing the work.

## The loop

```
   cohort-analysis ──→ h2m-variant ──→ library-design ──→ pegrna-oligo
   (which variants      (can it be      (pegRNAs with      (IDT order,
    matter clinically)   a mouse?)       sensors)           cloning)
          ↑                                                      │
          │                                                      ↓
   screen-stats ←── ngs-screen ←───────────────────────  the bench
   (LFC, FDR,       (FASTQ → counts        genotyping    (transfect,
    hit calls)       → editing %)       (.ab1, FlowJo)    screen)

   sci-plot styles every figure that comes out of any of these.
```

## Which tool for which question

| You want to know | Skill | Starts from | Ends at |
| --- | --- | --- | --- |
| How often a mutation occurs in patients; what co-occurs; survival | `cohort-analysis` | Nothing (public cohorts on disk) | Frequency / co-mutation / Cox table |
| Whether a human variant can be made in mouse, and how | `h2m-variant` | Variant list or MAF | Mouse coordinates + H2M class 0–6 |
| pegRNAs or base-editing sgRNAs for those variants | `library-design` | MAF-column variant table | Design table with scores and sensors |
| Oligos to order and the cloning | `pegrna-oligo` | Design table | IDT sheet, duplex maps, product plasmids |
| Just the duplex maps / just the assembly | `anneal` / `gga` | An existing order sheet | `.dna` files |
| Screen sequencing → per-guide counts and editing rates | `ngs-screen` | Paired FASTQs on the cluster | Count matrix + editing table |
| Which guides are enriched or depleted, with FDR | `screen-stats` | Count matrix | LFC-FDR table, hit list |
| Did the edit work in these colonies | `genotyping` | `.ab1` traces, FlowJo `.xls` | Editing %, zygosity calls |
| Make the figure publication-ready | `sci-plot` | Any result table | Styled PDF |

Boundaries that are easy to get wrong:

- **`ngs-screen` vs `genotyping`** — pooled screen FASTQs vs arrayed
  single-locus Sanger. Both answer "did the edit work", on different data.
- **`ngs-screen` vs `screen-stats`** — the first ends at the count matrix,
  the second starts there.
- **`library-design` vs `pegrna-oligo`** — the first ends at the design
  table, the second turns it into something you can order.
- **anything vs `sci-plot`** — analysis skills emit tables and quick QC
  panels; `sci-plot` owns the styling of anything going into a paper.

## What is on this machine

Roughly 100 GB under `$LAB_DATA_ROOT` (`/Users/Shared/lab-data/Database` on
the lab server). Nothing needs downloading.

- **Cohorts** — AACR-GENIE v13.1–v17.0, MSK-IMPACT 2017, MSK-HEME 2022,
  COSMIC v99, ClinVar 2024-02-06, gnomAD v4.1.0, MMRRC catalog
- **Genomes** — human GRCh37 / GRCh38, mouse GRCm38 / GRCm39
- **Annotation** — GENCODE v19 (GRCh37), v44 (GRCh38), vM33 (GRCm39)
- **Effect prediction** — AlphaMissense, SIFT, dbNSFP v5.2a
- **Lab-generated** — H2M precomputed database, MGI homology tables

`AACR-GENIE is controlled access` under a data use agreement. It is on the
machine for approved work; check before results derived from it leave the lab.

Details and the exact tree: `docs/data-layout.md`.

## Planning a multi-step request

When the ask spans several stages, say the plan before starting — the steps
have very different costs and the user should get to choose where to stop.

1. **Name the stages** and which skill does each.
2. **Say what the expensive steps are.** A `library-design` run over
   thousands of variants is hours; `ngs-screen` is hours of SLURM time on the
   cluster, not a laptop job. A cohort query is seconds.
3. **Check the input exists** before the first long step. Most failures are a
   missing column or an unset `LAB_DATA_ROOT`, and they surface an hour in.
4. **Hand off a file between stages**, not a reconstructed table. Each skill
   writes a CSV the next one reads.

A worked example, the CH/AML pilot:

```
cohort-analysis   rank CH genes by frequency in MSK-HEME        seconds
h2m-variant       which are modelable in mouse, class 0-1       ~1 min + 55 s load
library-design    pegRNAs with sensors and silent bystanders    hours
pegrna-oligo      IDT order split across two backbones          seconds
```

## Environments

Several skills need `pegg`, pinned to Python 3.9/3.10 with
`scikit-learn==1.1.1` and `numpy<2` — its scoring models are pickles that
break on anything newer. A wrong environment gives an unpickling error that
does not explain itself.

```bash
conda env create -f ~/Documents/GitHub/claude-bio-scientist/envs/pegg_env.yml
conda activate pegg_env
```

`envs/sensor_env.yml` and `envs/crispresso_env.yml` cover the screen
pipeline, and are kept separate because CRISPResso2 conflicts with the rest.

## When the answer is "not with these tools"

Say so plainly rather than improvising something adjacent. These tools do not
cover: single-cell, bulk RNA-seq quantification, ChIP/ATAC peak calling,
structural variant calling from reads, protein structure, or anything outside
cancer genomics and genome editing. Some of that exists in the lab's older
repos as one-off notebooks — `h2m-public/methods/`,
`PhD-FSR-MH-Lab/14_histone_conservation/` — but not as skills.

If a skill nearly fits but the input is wrong, say which column or file is
missing. That is almost always the real problem.
