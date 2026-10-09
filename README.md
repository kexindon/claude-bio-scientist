# claude-bio-scientist

Claude Code skills for cancer-genomics and genome-editing research in the
Sánchez-Rivera lab — clinical cohort → mouse variant modeling → guide library
design → oligo ordering → screen analysis → figures.

Each skill wraps code that already works and is already validated in a
published analysis. The skills do not reimplement the science; they make it
runnable by someone who did not write it.

## Just ask — you never name the skill

Open Claude Code in any directory and type the question the way you would say
it to a labmate. The right skill loads on its own.

### Clinical and cohort data

```
How common is DNMT3A R882H in AML?
What fraction of AML patients carry FLT3-ITD?

Pull every B-ALL variant in GENIE v17 and rank them by patient count.
Which genes co-occur with TP53 in B-ALL, and which are mutually exclusive?
Are SRSF2 and SF3B1 mutually exclusive?

Fit a Cox model for DNMT3A + RAD21 double mutants vs either alone.
Plot Kaplan-Meier for the four-arm stratification and give me the HR.
Rank the CH genes by clinical frequency so I can pick pilot targets.
```

Those run against AACR-GENIE and MSK-HEME under `$LAB_DATA_ROOT`, or against
your own patient table for the survival half. **AACR-GENIE is controlled
access** — see [Reference data](#reference-data).

### From a clinical variant to a mouse model

```
Take the top 20 CH mutations by GENIE frequency and tell me which ones
can actually be modeled in mouse.

What's the mouse equivalent of TP53 R273H?
What do these human variants correspond to in mouse, and are they modelable?

Design prime editing pegRNAs for the class 0 and class 1 hits,
with sensors and silent bystanders.

Turn that design table into an IDT order for UPEC-tevo.
```

That chain is `cohort_analysis` → `h2m_variant` → `library_design` →
`pegrna_oligo`, and you never have to say any of those names.

### Screens and bench readouts

```
Process the sequencing from last week's screen — config is in the run folder.
Call hits for bone marrow and spleen, ABE only, T0 baseline.
Which guides are enriched in meninges? Give me LFC and FDR.
Why are all my p-values exactly 1e-4?

Did the edit work? Here are the .ab1 traces from the colonies.
What is the editing efficiency, and which clones are homozygous?
What percent of the population is LSK in this FlowJo export?
```

### Figures

Every skill hands plotting to `sci-plot`, so figures look the same wherever
they come from:

```
Make that a volcano for the paper, Nature single-column width.
Restyle this to match figure 3c.
```

## Install

```bash
git clone https://github.com/kexindon/claude-bio-scientist.git
cd claude-bio-scientist
./install.sh
```

`install.sh` symlinks every skill in `skills/` into `~/.claude/skills/`, so a
`git pull` updates them all and nothing is ever copied out of sync. Run it
again after pulling if new skills were added.

To remove: `./install.sh --uninstall`.

## Skills

| Skill | Does | Status |
| --- | --- | --- |
| `pegrna_oligo` | pegRNA design table → IDT order sheet, with BsmBI/poly-T QC and correct scales | ported |
| `anneal` | Oligo order sheet → annotated SnapGene `.dna` duplex maps | ported |
| `gga` | In-silico Golden Gate assembly → annotated product plasmid maps | ported |
| `sci-plot` | Publication figure style for matplotlib/seaborn | ported |
| `h2m_variant` | Human → mouse variant modeling via `bioh2m` | planned |
| `cohort_analysis` | GENIE / MSK / COSMIC frequency, co-mutation, Cox PH, Kaplan-Meier | planned |
| `library_design` | pegRNA and base-editing library design via `pegg` | planned |
| `ngs_screen` | Sensor screen FASTQ → counts → CRISPResso2 → aggregation | planned |
| `screen_stats` | LFC, bootstrap null, two-sided empirical FDR | planned |
| `genotyping` | Sanger `.ab1` edit calling, flow cytometry readouts | planned |

[`docs/pegrna-skills.md`](docs/pegrna-skills.md) is the detailed reference for
the four ported skills — their CLI flags, the QC they run, and the wet-lab
facts they encode (SnapGene sticky-end flags, PEGG's RTT/PBS orientation, the
three cloning cases).

Designs for the planned skills, with their critiques, are in
[`docs/design/skill-designs.json`](docs/design/skill-designs.json) — including
145 documented gotchas pulled out of the real source code.

## Two conventions that keep the skills from colliding

**One cohort skill, not two.** `cohort_analysis` covers both the public
database query (frequency, co-occurrence) and the patient-level outcome
analysis (Cox, Kaplan-Meier). Splitting them produced two skills that both
answered "co-mutation" questions, so Claude could not tell them apart.

**Figures go through `sci-plot`.** No skill defines its own figure style. A
skill that produces a plot says so in its body and defers to `sci-plot` for
fonts, colors, axes and export — so every figure in the lab looks the same
and the style lives in exactly one place.

## Reference data

No skill hardcodes a data path. Every skill reads reference data through
environment variables, so the same skill works on a laptop and on the lab
server, under any username:

```bash
export LAB_DATA_ROOT=/Users/Shared/lab-data/Database    # lab server
export LAB_DATA_ROOT=~/Documents/GitHub/Database        # personal machine
```

`LAB_DATA_ROOT` holds the reference genomes, GENCODE annotation databases and
public cohort data (GENIE, COSMIC, ClinVar, MSK-IMPACT). See
[`docs/data-layout.md`](docs/data-layout.md) for the expected tree and how to
provision it.

**AACR-GENIE is controlled-access.** It sits under `LAB_DATA_ROOT` but is
governed by a data use agreement. Confirm your access covers your use before
querying it, and never commit any of it to this or any repo.

## Environments

Several skills depend on `pegg`, which is pinned to Python 3.9/3.10 with
`scikit-learn==1.1.1` and `numpy<2` — its scoring models are pickles that
break on newer versions. Environment files live in [`envs/`](envs/).

## Repo conventions

- One skill per directory under `skills/`, each with a `SKILL.md`.
- Helper code goes in that skill's `scripts/`, lifted verbatim from the
  validated notebook rather than rewritten, so the two cannot drift.
- `SKILL.md` frontmatter needs a `description` rich in trigger phrases —
  that text is the only thing Claude sees when deciding whether to load it.
- Never commit reference data, FASTQs, or screen outputs. `.gitignore`
  blocks the common extensions; keep it that way.
