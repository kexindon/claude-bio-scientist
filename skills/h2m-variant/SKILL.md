---
name: h2m_variant
description: >-
  Convert variants between human and mouse at codon resolution with the local
  bioh2m package, and triage each into H2M class 0-6 so you know which can
  actually become a mouse. Use for "model this human variant in mouse", "mouse
  equivalent of TP53 R273H", "mouse-ify these mutations", "can this be modeled
  in mice", "run H2M on this MAF", "which of these can I make a mouse for",
  "is this base-editable in mouse", "mouse ortholog and canonical transcript".
  Also does mouse-to-human and paralogs. Not for pegRNA design
  (library_design), cohort frequency (cohort_analysis), or
  genotyping primers (genotyping).
---

# h2m_variant — human variant in, mouse variant out, with a verdict

## The class system: read this first

Every modeled variant comes back with an integer `class`. **This is the output.**
The coordinates matter only once you know the class.

| class | statement | what it means for you |
|---|---|---|
| **0** | can be originally modeled | **Orderable.** The identical nucleotide change works in mouse. `*_m` and `*_m_ori` are the same variant. |
| **1** | can be alternatively modeled | **Orderable.** A *different* mouse codon gives the same amino acid change. Order `*_m`. **Not** `*_m_ori` — that one is the identical-nucleotide change and makes a different protein. |
| **2** | can be modeled, effect may not be consistent | **Judgement call.** Frameshift, splice-site and non-coding variants. H2M places a change at the orthologous position but does not vouch for the protein consequence, and no amino-acid-matched alternative exists. Read `HGVSp_m` yourself. |
| **3** | not originally modelable, no alternative found | **Drop it.** Model a different variant in the same gene. |
| **4** | mutated sequences are not identical | **Alignment problem, not your fault.** Human and mouse proteins do not align contiguously across this site, so the orthologous position is ambiguous. ~14% of real GENIE coding variants. |
| **5** | coordinate error, not in the query gene | **Your input is wrong.** Almost always GRCh37/GRCh38 confusion between your coordinates, `--ver`, and the GENCODE `.db`. |
| **6** | cannot be originally modeled | Only appears with `--no-alternative`, or under `--param BE` with the alternative search off. With defaults you get class 3 instead. |

**Buildable = class 0 or 1.** On a real AACR-GENIE v17 coding subset (TP53 / KIT /
DNMT3A, 49 variants, verified run) the split was **0:32, 1:7, 2:3, 4:7** — about
**80% orderable**. Budget for roughly that.

Two rules that are not obvious and will bite you:

- **Filter on `class`, never on `status`.** `status` is a module-level global
  that is not reset consistently: `module.py` assigns class 2 at three places
  (L861, L1689, L1786, L1907) and only the first sets `status=False`, so a
  class-2 variant usually arrives with `status == True`. `class` is reliable.
- **`match == False` has no class of its own.** If your reference allele
  disagrees with the genome, H2M sets `match=False`, *silently substitutes the
  genome's allele*, and still reports class 0. You get a clean-looking result
  for a variant you did not ask about. `h2m_triage.py` flags these and excludes
  them from the hand-off by default.

## When to use / when not to

Use for anything that is "the same variant, in the other species, at nucleotide
resolution" — or for the prior question of whether that is even possible.

Not for: pegRNA/sgRNA design and designability (`library_design`); cohort
frequency, recurrence, or the published **human2mouse.com** database
(`cohort_query` — "look this up in H2M" is ambiguous and may mean the website,
so ask); genotyping primers for the finished allele (`genotyping`); sequencing
readouts (`ngs_screen`, `screen_stats`).

This skill **ends by writing a CSV** and naming `library_design` as the next
step. "Model these in mouse and design the pegRNAs" is one task to the user but
two skills in two conda environments — say so rather than stopping silently.

## Environment

`bioh2m` is an installed pip package. **Do not vendor it; call it.**

```bash
pip install bioh2m          # 1.0.0 verified here; the lab cluster env pins python=3.10
export LAB_DATA_ROOT=~/Documents/GitHub/Database
python ~/.claude/skills/h2m_variant/scripts/h2m_refs.py      # verify the tree
```

This is **a separate environment from `pegg`** and cannot be merged with it:
pegg needs `scikit-learn==1.1.1` with `numpy<2` for its pickled scoring models,
bioh2m wants `numpy>=1.24.2`. The hand-off to `library_design` is therefore a
file on disk, not an import.

Reference data, all under `$LAB_DATA_ROOT` (verified present on disk):

```
RefGenome/ncbi-2023-09-12/GCF_000001405.25_GRCh37.p13_genomic.fna.gz
RefGenome/ncbi-2023-09-12/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
RefGenome/mouse-2023-09-13/GCF_000001635.27_GRCm39_genomic.fna.gz
Genecode/gencode_v19_GRCh37.db      <- note the spelling: Genecode
Genecode/gencode_v44_GRCh38.db
Genecode/gencode_vm33_GRCm39.db
```

`docs/data-layout.md` shows the `.db` files inside the per-release
subdirectories; on the machines checked they are at the top of `Genecode/`
instead. `h2m_refs.py` looks in both. The `.db` files are **not** bioh2m
downloads — someone built them with `gffutils.create_db` from a GENCODE
`.gtf.gz` (the `.gtf.gz` sources are in `Genecode/release_*/`, ~20 min each).

## The fast path — no genome, no network, ~2 seconds

Most "what's the mouse gene for X" questions need no reference data at all;
`h2m.query` and `h2m.get_tx_batch` read JSON packaged inside the wheel.
**Never route these through a genome load.**

```bash
python scripts/h2m_prepare_input.py --genes TP53,KIT,SMARCA4 -o /tmp/q
# -> gene_name_h, gene_name_m, homology_type, similarity, tx_id_h, tx_id_m
# TP53 -> Trp53 (one2one, 77.4%)   KIT -> Kit (one2one, 82.8%)
```

Add `--direction m2h` to go the other way (mouse gene name only — a mouse gene
*id* raises).

## The full procedure

### 1. Build the input table

```bash
python scripts/h2m_prepare_input.py VARIANTS.csv --source maf --ver 37 -o OUT
```

`--source` dispatches the bioh2m reader: `cbio` for cBioPortal / AACR-GENIE /
TCGA `data_mutations*.txt`, `gnomad` for a gnomAD browser CSV export (**add
`--gene TP53`** — a gnomAD export is a single-gene download with no gene
column), `clinvar` with `--clinvar-ids 375926,140953`, or `maf` for a plain CSV
you wrote yourself. It then chains `get_tx_batch(h)` → `query_batch` →
`get_tx_batch(m)` and writes `OUT_h2m_input.csv` plus `OUT_dropped.csv`.

**Read the dropped file.** Genes missing from the packaged dictionaries are
silently removed at each step, and the row count only goes down.

Required columns for the modeling step, exactly:
`start_h, end_h, type_h, ref_seq_h, alt_seq_h, tx_id_h, tx_id_m, index`.
Coordinates are 1-based, on the **positive strand of the chromosome** whatever
strand the transcript is on, with `start <= end`. `type_h` is one of
SNP/DNP/TNP/ONP/INS/DEL — not a MAF `Variant_Classification` like
`Missense_Mutation`. Empty alleles are `-`.

### 2. Model

```bash
python scripts/h2m_model.py OUT_h2m_input.csv -o RUN --ver 37
python scripts/h2m_model.py OUT_h2m_input.csv -o RUN --ver 37 --dry-run   # 1 s check
```

`--dry-run` validates columns and reference paths without loading anything —
use it first, so a typo costs a second instead of a minute. Flags:
`--direction {h2m,m2h,h2h,m2m}`, `--param BE`, `--coor aa`,
`--max-alternative N`, `--memory-size N`, `--no-alternative`, `--bind`,
`--show-sequence`.

Writes `RUN_result.csv` and `RUN_leftover.csv`. **Always open the leftover
file** — see the gotchas.

### 3. Triage and hand off

```bash
python scripts/h2m_triage.py RUN_result.csv -o FINAL
python scripts/h2m_triage.py RUN_result.csv --show 7     # every alternative for one variant
```

Prints the class table with counts, flags `match == False`, flags variants in
poorly conserved local context, and writes `FINAL_modeled_variants.csv` —
classes 0/1 only, `match == False` excluded, just the columns the design step
needs. Hand that path to `library_design`.

## Single variant, interactively

```python
import bioh2m as h2m, sys
sys.path.insert(0, "~/.claude/skills/h2m_variant/scripts")
import h2m_refs
rh, ih, rm, im, dbh, dbm = h2m_refs.load_refs(37)      # ~1 min, once per session

# TP53 R273H, GRCh37
r = h2m.model(rh, ih, rm, im, dbh, dbm,
              "ENST00000269305.4", "ENSMUST00000108658.10",
              7577120, 7577120, "C", "T", ty_h="SNP", ver=37)
r[0]["class"], r[0]["HGVSp_h"], r[0]["HGVSp_m"]        # (0, 'R273H', 'R270H')
```

Amino-acid input instead of coordinates — pass the residue number as both
`start` and `end`:

```python
h2m.model(..., 175, 175, "R", "H", coor="aa", ty_h="SNP", ver=37)   # -> R172H in mouse
```

The reverse question — "our line is Trp53 R172H, what human mutation is that?"
— is `m2h` and `coor='aa'` **together**, which is the common lab case:

```python
h2m.model(..., 172, 172, "R", "H", coor="aa", ty_h="SNP", ver=37, direction="m2h")
```

Paralogs within one species pass the same genome and db twice
(`direction='h2h'`, e.g. SMARCA2→SMARCA4); output columns are then suffixed
`_1`/`_2`, not `_h`/`_m`.

## Base editing

`--param BE` keeps only changes a base editor could install. The gate is
stricter than "it is a transition": `BE_check` needs equal lengths, *all*
reference bases identical, *all* alternate bases identical, and the pair in
`[(A,G),(T,C),(C,T),(G,A)]`. So `A>G`, `AA>GG` and `CCC>TTT` pass; `AC>GT` and
every indel do not. It also returns BE alternatives for variants that were
already class 0.

It says nothing about whether a PAM is in range — that is `library_design`.

## Figures

The per-variant alignment plot is `h2m.visualization(model_result, flank_size,
print_size)`, and the model call **must** have used `show_sequence=True` or it
raises `KeyError: 'seq_h'`. This is a fixed-layout plot from the package, not a
restyle target — it is **not** a `sci-plot` deliverable, so do not try to
restyle it.

Anything else you plot from these results — a class-distribution bar chart, a
modelability-by-gene panel — **is** a `sci-plot` job. Load that skill for fonts,
colors and sizes; keep your code to what goes on which axis.

## Gotchas

- **`h2m.model` returns a LIST, not a dict.** One entry per alternative, up to
  `max_alternative=5`. `pd.DataFrame(result)` gives you N rows for ONE variant,
  and `model_batch` flattens all of them, so counting result rows overcounts
  variants. A real run: 58 rows for 49 variants. Group on `index` first —
  `h2m_triage.py` does.

- **`model_batch` swallows every per-variant exception.** Both the alignment
  step and each `model` call sit inside a bare `try/except`; failures are
  appended to the second returned DataFrame and the function prints only
  "There were rows that could not be processed." No traceback, no reason. Read
  the leftover file and re-run those rows one at a time through `h2m.model` to
  see the real error.

- **There is a real bug behind some of those leftovers.** In `FindAlternative`
  (module.py L1637-1639) the class-3 branch passes `flank_size_left` /
  `flank_size_right`, which are only assigned inside the loop that just failed
  to find an alternative — so a genuine class 3 raises
  `UnboundLocalError: local variable 'flank_size_left' referenced before
  assignment` instead of being reported. Reproduced on GENIE: KIT
  chr17:55594199 G>T, GRCh37. **Workaround:** re-run the leftover rows with
  `--no-alternative`, which reports them as class 6. Do not "fix" this in the
  skill — it is upstream in bioh2m.

- **`flank_size_left` / `flank_size_right` are strings, not numbers**, and the
  unit changes with the variant: `'4aa'` for a coding variant, `'1nt'` for a
  splice-site one. Parse the digits off before comparing. They report identical
  flanking context shared by human and mouse, which is the field to cite when
  judging how far to trust a modeling — `KIT A533T` came back `1aa / 0aa`.

- **A GRCh37/38 mix-up gives you class 5 for every row, a minute in.** Nothing
  in bioh2m cross-checks your coordinates against `ver`, and the two assemblies
  differ by 67 kb (EGFR) to 1.7 Mb (PTEN). **gnomAD v4 exports are GRCh38;
  AACR-GENIE and the lab's ClinVar file are GRCh37.** `h2m_prepare_input.py`
  warns when every sampled coordinate is >50 kb from its own transcript.

- **`genome_loader` is the whole cost.** It parses the entire FASTA into
  Biopython records (~4 GB for human) and you need both species at once;
  measured ~55 s for the pair on a local SSD, so budget 8-16 GB of RAM. Load
  once per session via `h2m_refs.load_refs`, never inside a loop.

- **Chromosomes are handled positionally and X/Y differ by species** — human
  X=22, Y=23; mouse X=19, Y=20, indexing `records[index_list[chrom-1]]`.
  Passing the mouse genome where the human one belongs does not error; it
  returns sequence from the wrong chromosome. This is also why the FASTA must
  be the exact NCBI file (contigs are filtered by description keyword), not a
  differently sorted assembly.

- **One-to-many orthologs multiply your rows.** `query_batch` emits one row per
  human↔mouse pair and merges with `how='outer'`, so a gene with two mouse
  orthologs becomes two rows, each then up to `max_alternative` rows after
  modeling. The outer merge can also produce all-NaN variant rows;
  `h2m_prepare_input.py` drops those and says how many.

- **`index` must be a unique column, not the pandas index.** It is the join key
  for alternatives and for `--bind`; duplicates turn the bind merge into a
  cartesian product. `model_batch` raises "Error 2" if it is missing.

- **`query` raises instead of returning empty** for an unknown human gene, and
  in `m2h` mode it indexes `homo_dict_mouse` directly, so an unknown mouse gene
  is a bare `KeyError`. Both are caught in `h2m_prepare_input.lookup_genes`.

- **`get_tx_batch` is offline but frozen** at the installed wheel's version: it
  reads packaged `dict_tx_h.json` / `dict_tx_m.json`, picking index `[0]` for
  ver=37 and `[1]` for ver=38 — positionally, so a wrong `ver` silently returns
  a transcript from the other assembly. A gene renamed since the wheel was
  built lands in the failed frame. `get_tx_id` is the live Ensembl alternative;
  it needs internet, is one round-trip per gene, and returns
  `[chr, start, end, canonical_tx, all_tx, tx_dicts]` — **the canonical
  transcript is element 3**, and element 4 is the full list reordered with the
  canonical first.

- **A `.db` from a different GENCODE release fails as a bare `KeyError`** inside
  `GetTx` (`db_input[tx_id]`), because its transcript IDs and versions do not
  match what the wheel's dictionaries return. `h2m_refs.preflight()` catches
  this in a second: it takes a transcript from the packaged dict and asserts
  the `.db` contains it.

- **INS coordinates are MAF convention, not VCF.** `check_mut_info` shifts
  `start -= 1` / `end += 1` for an INS with a non-empty ref and `final_return`
  shifts back, so in and out you use flanking positions with `ref_seq='-'`.
  Never hand VCF-style indels to `model` directly — run them through
  `vcf_to_maf` or `clinvar_to_maf`, which do the base-trimming correctly.

- **`get_variant_type` mutates its input in place** and fills NaN alleles with
  `''`, so a zero-length ref reads as INS/DEL instead of erroring. Pass a
  `.copy()`. It also makes no frame check — In_Frame vs Frame_Shifting is
  decided later inside `model`.

- **`memory_protect=True` raises "Alarm 3"** when either protein exceeds
  `memory_size` (default 10000 aa). TTN, DMD and OBSCN hit this. The lab's
  production run uses `memory_size=2000` to fail fast rather than risk the
  node, and that is this skill's default too. Lowering the number is
  *stricter*; `memory_protect=False` is what removes the cap.

- **`model` is not thread-safe.** It carries its class decision between gate
  functions through module-level `global status, error, statement`. Scale out
  with the SLURM array pattern in `methods/run_h2m/task.sh` — split the input
  CSV, one process per shard — never with threads or a pool sharing the loaded
  genome.

- **Importing bioh2m calls `warnings.filterwarnings("ignore")` globally**, for
  the whole process. Any pandas or numpy warning you expected from unrelated
  analysis code in the same notebook is now silenced.

- **`methods/run_h2m/h2m_db.py` still says `import h2m`**, the pre-rename
  module name, which no longer resolves. The package is `bioh2m`, imported as
  `import bioh2m as h2m`. Lift the logic, not that line.

## Verifying a run

```bash
python scripts/h2m_triage.py RUN_result.csv -o FINAL
```

What to check, in order:

1. `RUN_leftover.csv` is empty. If not, re-run those rows singly, and with
   `--no-alternative` for the `flank_size_left` bug.
2. No `match == False` rows survive into the hand-off.
3. The class distribution is plausible — all class 5 means an assembly
   mismatch; all class 4 means the wrong transcript pair.
4. Spot-check one class 1 with `--show <index>` and confirm `HGVSp_m` equals
   the human amino acid change while `HGVSp_m_ori` does not. Verified example:
   TP53 R306Q → mouse `HGVSp_m = R303Q` via `c.907_908AG>CA`, while
   `HGVSp_m_ori = R303K`. Ordering `*_m_ori` there would build the wrong allele.

## Cluster scale

A few hundred variants across a handful of genes is seconds to a minute once
loaded — cost tracks the number of distinct `(tx_id_h, tx_id_m)` pairs, not the
variant count, because the protein alignment is computed once per pair. Thousands
of variants across thousands of genes is hours, and the lab runs it as a 16-way
SLURM array: `methods/run_h2m/task.sh` + `h2m_env.yml` (python=3.10) +
`h2m_config.txt`, splitting the input CSV and giving each shard its own process.
Those files hardcode `/net/bmc-lab2/data/lab/sanchezrivera` paths — read them
for the pattern, do not copy the paths.

## Citation

Dong, K., Gould, S.I., Li, M. et al. *Nat Biotechnol* (2025).
doi:10.1038/s41587-025-02925-0. Portal: human2mouse.com
