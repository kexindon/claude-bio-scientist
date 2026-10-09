# Adding a skill to this repo

## Layout

```
skills/<skill_name>/
├── SKILL.md              required — frontmatter + instructions
├── scripts/              helper code the skill runs
└── references/           optional longer reference the skill points Claude at
```

Directory name = skill name = the `name:` in frontmatter. Use
`lower_snake_case`.

## Frontmatter

```yaml
---
name: screen_stats
description: >
  Compute log-fold changes and empirical FDR for a pooled base- or
  prime-editing sensor screen — RPM normalization, replicate LFC against a
  T0 or barcoding baseline, a bootstrap null pooled from safe-targeting and
  non-targeting control guides, and two-sided empirical p-values combined
  across replicates. Use whenever asked to "call hits", "compute LFC",
  "run the FDR", "which guides are enriched", "make a volcano for the
  screen", or to analyze a count matrix from the sensor pipeline.
---
```

**The description is the whole trigger mechanism.** It is the only text
Claude sees when deciding whether to load the skill, so it must say both what
the skill does *and* the words a user would actually type. A description that
restates the title will not fire.

Write it in the third person. Pack in real trigger phrases — the ones from
your own lab Slack, not invented ones. If people ask in Chinese as well as
English, the English description still triggers, but naming the domain terms
both ways helps.

Check it against the neighbours: if two skills could both plausibly answer
"analyze my screen", say in each description what it is *not* for.

## Body

Write for a competent biologist who does not know this codebase. Procedural,
specific, and honest about what breaks.

Cover:

1. **When to use / when not to** — the boundary with neighbouring skills
2. **Inputs** — exact file formats and required column names
3. **The procedure** — numbered steps, with the real commands
4. **Gotchas** — version pins, silent failure modes, ordering requirements
5. **Verification** — how the user knows the output is right

Point at `scripts/` for anything longer than a few lines rather than
inlining it.

## Rules specific to this repo

**Never hardcode a data path.** Read `$LAB_DATA_ROOT`:

```python
import os
DATA = os.environ["LAB_DATA_ROOT"]
db = gffutils.FeatureDB(f"{DATA}/Genecode/release_m33_GRCm39/gencode_vm33_GRCm39.db")
```

Fail with a clear message if it is unset — not with a `KeyError`:

```python
DATA = os.environ.get("LAB_DATA_ROOT")
if not DATA:
    raise SystemExit("LAB_DATA_ROOT is not set — see docs/data-layout.md")
```

**Lift working code verbatim.** If a notebook already does this correctly and
has been used for a real order or a published figure, copy the functions into
`scripts/` unchanged rather than rewriting them. Note the source file at the
top so the two can be compared later. `pegrna_oligo` does this — its engine
came straight out of `pegrna_gga_oligo_order.ipynb`.

**State the environment.** If the skill needs `pegg`, say that it needs
Python 3.9/3.10 with `scikit-learn==1.1.1` and `numpy<2`, and point at
`envs/`. A skill that silently assumes an environment fails confusingly.

**No data in the repo.** No FASTQs, no genomes, no cohort files. Test
fixtures go in `tests/fixtures/` and stay small.

## Installing and testing

```bash
./install.sh          # symlink every skill into ~/.claude/skills/
```

Restart Claude Code. Then test the trigger by asking for the task in your own
words *without* naming the skill — if it does not load, the description is
the problem, not the body.
