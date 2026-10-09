---
name: pegrna_oligo
description: Turn a PEGG/PRIDICT pegRNA design table into everything needed to order and clone — an IDT oligo order sheet (Golden Gate protospacer/scaffold/3'-extension oligos with BsmBI QC, dedup and correct IDT scales), the annealed oligo duplex .dna maps, and the in-silico assembly product plasmid .dna maps. Different groups can go into different backbones in one continuously numbered order (e.g. CH into UPEC-tevo, AML into UPEmS-tevo). Use whenever asked to design or order pegRNA oligos, build a pegRNA Golden Gate order, anneal pegRNA oligos, or simulate cloning pegRNAs into UPEmS/UPEC vectors.
---

# pegrna_oligo — designs → order sheet → annealed duplexes → plasmids

One command covers the whole path from a design table to the bench:

```bash
python ~/.claude/skills/pegrna_oligo/scripts/pegrna_oligo.py DESIGNS.csv -o OUTDIR \
    --file-name hsc_pilot_ordered \
    --group-col gene_group --order CH,AML \
    --group 'CH=UPEC-tevo:/path/UPEC_tevo.dna' \
    --group 'AML=UPEmS-tevo:/path/UPEmS_tevo.dna'
```

```
OUTDIR/
├── pegrna_orders/<name>_oligos.csv    IDT upload: Name,Sequence,Scale,Purification
├── pegrna_orders/<name>_report.csv    per design: Group, name_stem, spacer, QC
├── pegrna_orders/<name>_designs.csv   the designs actually ordered
├── oligos/<stem>_P.dna, _E.dna        annealed duplexes, real 4-nt sticky ends
└── plasmids/<stem>.dna                circular assembly products
```

| flag | meaning |
|---|---|
| `--group G=BACKBONE:/path/map.dna` | repeatable; routes each group to its own vector. A bare `BACKBONE:/path` applies to all designs |
| `--group-col` | column holding the group label (omit for a single group) |
| `--order CH,AML` | index block order, so CH numbers come before AML |
| `--name` | oligo name prefix; `""` (default) → `index_gene_mut_P_T` |
| `--gene-col` / `--mut-col` | columns used to label oligos (`None` to omit) |
| `--epeg oligo\|none` | required for plain UPEmS/UPEC; leave unset for `-tevo` |
| `--plate` | add Well Position as the first column |
| `--no-anneal` / `--no-plasmids` | stop after the order sheet |

## Oligo naming

`NAME_index_gene_mut_CODE_TAG`, e.g. `7_Trp53_R245W_E_T`. `P`=protospacer,
`E`=3' extension, `S`=scaffold; `T`/`B` = top/bottom. Mutation labels are
sanitized for IDT (`p.R451*` → `R451X`, `D594_F595insNS` → `D594F595insNS`) and
trimmed to IDT's 30-character limit, mutation first then gene — the index alone
keeps names unique, so trimming is safe.

## The three cloning cases

| backbone | `--epeg` | up OH | down OH | 3' ext top | 3' ext bottom |
|---|---|---|---|---|---|
| `UPEmS-tevo`, `UPEC-tevo` | *(unset)* | CACC | CGCG | `gtgc`+ext | `cgcg`+rc(ext) |
| `UPEmS`, `UPEC` | `oligo` | CACC | TTTT | `gtgc`+ext+tevo | `aaaa`+rc(tevo)+rc(ext) |
| `UPEmS`, `UPEC` | `none` | CACC | TTTT | `gtgc`+ext | `aaaa`+rc(ext) |

Impossible backbone/epeg combinations raise rather than producing dead oligos.
`CGCG` is the first 4 nt of tevopreQ1, which is why the `-tevo` backbones join
there.

## What it checks

- **BsmBI scan** of protospacer, extension and every junction — a `CGTCTC`/
  `GAGACG` inside the RTT/PBS silently destroys the assembly, so it raises.
- **Overhang compatibility** at every junction, against the vector's real map.
- **Dedup** — designs sharing a protospacer share `P_T`/`P_B`, ordered once; the
  report's `reused` column maps the rest back.
- **IDT scale by length** (25nm ≤60 nt, 100nm ≤90, 250nm ≤100, else Ultramer).
- **In-silico assembly** — the product must contain
  `CACC spacer GTTTTAGAG scaffold GTGC RTT PBS [tevopreQ1]` with zero surviving
  BsmBI sites, scanned across the circular join.

## PEGG orientation (verified, do not "fix")

PEGG's `RTT` and `PBS` are **already in the right orientation to concatenate**:

- `RTT_PBS == RTT + PBS` — 5'-RTT-PBS-3', the natural pegRNA 3' extension order.
- `RTT` appears **as-is** in `sensor_alt` (not reverse-complemented).
- `rc(PBS)` — not `PBS` — is found inside the protospacer, ending 3 nt from its
  3' end (the PE2 nick). PBS anneals to the nicked DNA, so this is correct.

So `E_T = "gtgc" + RTT_PBS` with **no** reverse-complementing. Sanity check on
new data: `protospacer.find(rc(PBS)) >= 0` and `protospacer.find(PBS) < 0`.

## Colors in the .dna output

spacer cyan `#00ccff`, RTT amber `#ffcc00`, PBS green `#99cc00`, **all overhangs
and junctions red** `#ff0000`, scaffold grey, tevopreQ1 mauve, cassette bisque.

Two conventions that matter:

- **`GTTTTAGAG` is sgRNA scaffold nt 1–9, not a linker.** `GTTTT` is
  double-stranded, `AGAG` is the sticky end. The scaffold feature spans all 9 nt
  and the red overhang is drawn over its last 4; they overlap deliberately.
- **Duplex files must be double-stranded with real sticky ends.** Flags byte
  `0x02` (bit 1 = double-stranded; `0x00` opens single-stranded), and block 8
  `UpstreamStickiness`/`DownstreamStickiness` = 4. Without block 8 both default
  to 0 = **blunt**, and the overhangs are only red paint. Plasmids are `0x03`
  (circular + double-stranded) with stickiness 0.
  Reference: `~/Desktop/HSC x Order/Assembly/pegRNA_scaffold.dna`.

## Verifying

`scripts/pegoligo.py` is lifted verbatim from the validated notebook
`pegrna_gga_oligo_order.ipynb`, so the two cannot drift. Output of this skill was
checked **byte-identical** against that notebook's order sheet plus the `anneal`
and `gga` skills' .dna files.

Worth asserting after a run:

```python
seq, feats = read_dna(f'{stem}_P.dna')
assert read_stickiness(f'{stem}_P.dna') == (4, 4)      # not (0, 0)
assert seq[:len(T)].upper() == T.upper()               # top strand == ordered _T
assert rc(seq)[:len(B)].upper() == B.upper()           # bottom strand == ordered _B
```

## Golden Gate bench sheet (`--gga-plan`)

Adds `gga_plan/gga_plan.tsv` — the cloning planning sheet, laid out like the
lab's "02 Golden Gate Assembly (pegRNAs)" workbook:

```bash
python scripts/pegrna_oligo.py DESIGNS.csv -o OUTDIR ... \
    --gga-plan --start-tube 79 --plan-title "HSC Project Pilot pegRNAs"
```

or standalone from an existing order:

```bash
python scripts/gga_plan.py REPORT.csv -o OUTDIR --start-tube 79 \
    --title "HSC Pilot" --vector "UPEmS-tevo (uL)"
```

It writes the annealing master mix and per-tube layout (one `P` and one `E` row
per pegRNA, in racks of 8), then the Golden Gate master mix and layout. **Open
`gga_plan.tsv`, select all, copy, and paste into a new tab** — tab-separated text
lands in the right cells. `gga_plan_anneal.csv` / `gga_plan_gga.csv` are the same
two tables as flat CSV for scripting.

`--start-tube N` is the inventory index of this order's **first** oligo; the rest
are numbered from there in order-sheet row order (`P_T, P_B, E_T, E_B` per
design), so tube *k*'s four oligos are `N+4k .. N+4k+3`. Omit it and the sheet
labels wells by full oligo name instead.

```{note}
When designs go into more than one backbone, the Golden Gate master mix is split
into one column per vector (`UPEC-tevo (n=10)`, `UPEmS-tevo (n=6)`) — they are
physically separate mixes, so a single pooled column would be wrong at the bench.
```

The Drive connector can read a sheet and create files but cannot write cells or
add a tab, so this skill cannot push into an existing Google Sheet. Paste the TSV,
or use a Sheets-API script / Apps Script if you want it fully hands-off.

## Related skills

`anneal` (duplex maps only) and `gga` (plasmids only) do the second and third
steps standalone, for when an order sheet already exists. This skill runs all
three and shares their writer.

## Wet-lab hand-off

1. Resuspend to 100 µM. 2. Anneal P and E pairs with T4 DNA Ligase (no PNK);
scaffold stock is annealed with T4 PNK. 37 °C 30 min → 95 °C 5 min → ramp to
25 °C at 5 °C/min. 3. Dilute 1:100. 4. Golden Gate: 1.5 µL vector (75 ng/µL) +
1 µL each duplex + 2 µL 10× T4 buffer + 1 µL NEB Golden Gate Mix (BsmBI) +
12.5 µL H₂O; `(42 °C 1 min → 16 °C 1 min) × 60 → 60 °C 5 min`.

```{note}
The `CGCG` junction of the `-tevo` backbones is palindromic, so the 3' extension
duplex can ligate to itself. High background of small colonies on a `-tevo`
cloning is likely that; the `TTTT`/`AAAA` junction of plain UPEmS/UPEC does not
have this issue.
```
