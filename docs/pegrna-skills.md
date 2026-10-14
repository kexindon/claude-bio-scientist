# The pegRNA and figure skills

Reference for the four skills that were written before this repo existed, and
for the wet-lab facts they encode. The [README](../README.md) lists every
skill; this file is the detail for these four.

Skills load automatically when a request matches — you do not have to name
them, though `/skill-name` works.

| skill | does | reach for it when |
|---|---|---|
| [`pegrna-oligo`](#pegrna-oligo) | pegRNA designs → IDT order sheet + annealed duplex maps + product plasmids + Golden Gate bench sheet | ordering or cloning pegRNAs |
| [`anneal`](#anneal) | oligo order sheet → annealed duplex `.dna` | you already have an order sheet |
| [`gga`](#gga) | oligo order sheet → assembly product plasmid `.dna` | simulating/verifying the cloning |
| [`sci-plot`](#sci-plot) | publication figure style for matplotlib/seaborn | any figure panel for a paper |

`pegrna-oligo` runs the whole pegRNA path; `anneal` and `gga` are its second and
third steps standalone. All three share one SnapGene writer.

The script paths below use `~/.claude/skills/<name>/`, which is where
`install.sh` links them. They resolve to `skills/<name>/` in this repo.

---

## pegrna-oligo

PEGG/PRIDICT design table → everything needed to order and clone.

```bash
python ~/.claude/skills/pegrna-oligo/scripts/pegrna_oligo.py DESIGNS.csv -o OUTDIR \
    --file-name my_order \
    --group-col gene_group --order CH,AML \
    --group 'CH=UPEC-tevo:/path/UPEC_tevo.dna' \
    --group 'AML=UPEmS-tevo:/path/UPEmS_tevo.dna'
```

```
OUTDIR/pegrna_orders/*_oligos.csv    IDT upload (Name,Sequence,Scale,Purification)
                    /*_report.csv    per design: Group, name_stem, spacer, QC
                    /*_designs.csv   designs actually ordered
       oligos/<stem>_P.dna, _E.dna   annealed duplexes, real 4-nt sticky ends
       plasmids/<stem>.dna           circular assembly products
       gga_plan/gga_plan.tsv         bench sheet (with --gga-plan)
```

Groups can go into different backbones in **one continuously numbered order** —
CH into UPEC-tevo and AML into UPEmS-tevo, checked separately, exported as a
single sheet. `--no-anneal` / `--no-plasmids` stop early.

**Checks before anything is written:** BsmBI scan of protospacer, extension and
every junction; overhang compatibility against the vector's real `.dna` map;
dedup of shared protospacers; IDT scale by length; and a full in-silico assembly
whose product must contain `CACC spacer GTTTTAGAG scaffold GTGC RTT PBS
[tevopreQ1]` with zero surviving BsmBI sites.

### Naming

Oligo names are `index[arch]_gene_mut_CODE_TAG` — `7_Trp53_R245W_E_T`,
`17t_Tet2_S1412*_E_B`.

**The index identifies the mutation, not the row.** `--index-by mutation` (the
default) keeps one number per gene+mut, so a pegRNA and its epegRNA are both
design 17. `--index-seed EARLIER_REPORT.csv` carries the numbering across orders,
so a mutation ordered in September keeps its index in October. `--index-by row`
restores the old per-row counter.

Because the index is shared, the **architecture** marker is what distinguishes
the two arms: `17_…` plain, `17t_…` with the tevopreQ1 motif. Without it the
plain and tevo `E_B` oligos — which differ only in their downstream overhang
(`aaaa` vs `cgcg`) — would collide on one name, and two different molecules must
never share a name on an order sheet. `--vector-code` additionally inserts the
backbone letter (`17Ct`, C=UPEC / S=UPEmS); the September order uses it because
its tubes are already labelled that way, new orders do not.

**`*` is kept, not rewritten to `X`.** `Tet2_L35*` is the mutation's actual
HGVS name. No IDT documentation restricts the name charset, and IDT's own
bulk-upload parser takes the Name field raw — no validation, no truncation. Set
`STOP_AS_X = True` (or `stop_as_x=True`) only for a vendor form that rejects it.

Note the asterisk *is* reserved in IDT's **sequence** field, where it denotes a
phosphorothioate bond — so an `*` in a sequence silently orders a different,
costlier molecule. That hazard is sequence-side only and is most likely where
the "IDT forbids `*`" folklore came from.

`IDT_NAME_LIMIT = 30` is a **tube-label budget, not an IDT rule** (the web form
allows 80). It is kept so old and new sheets agree; raise it per call via
`limit=` if a label needs the room.

**Plasmid maps get the full vector name, oligos stay compact.** Filed on the
Desktop as `2_UPEC_Dnmt3a_R878H_tevo.dna`, because a plasmid name is read in
SnapGene where there is room for it; the same name on an oligo
(`12_UPEmS_Flt3_D594F595insNS_tevo_P_T` = 36 chars) blows the label budget, so
oligos keep `12t_Flt3_D594F595insNS_P_T`. Both carry the same index, so they
cross-reference.

The engine `scripts/pegoligo.py` is lifted **verbatim** from the validated
notebook `pegrna_gga_oligo_order.ipynb`, so the two cannot drift.

### Golden Gate bench sheet

`--gga-plan --start-tube 79 --plan-title "HSC Pilot"` additionally writes
`gga_plan/gga_plan.tsv`, laid out like the lab's *02 Golden Gate Assembly
(pegRNAs)* workbook: annealing master mix + per-tube layout (one `P` and one `E`
row per pegRNA, in racks of 8), then the Golden Gate master mix + layout. **Open
the TSV, select all, copy, paste into a new tab** — it lands in the right cells.

`--start-tube N` is the inventory index of this order's **first** oligo; the rest
follow in order-sheet row order (`P_T, P_B, E_T, E_B`), so design *k* uses
`N+4k .. N+4k+3`. Omit it to label by oligo name instead.

Also available standalone from an existing order:
`python ~/.claude/skills/pegrna-oligo/scripts/gga_plan.py REPORT.csv -o OUTDIR --start-tube 79`

When designs span two backbones the GGA master mix is split into one column per
vector — they are physically separate mixes, so one pooled column would be wrong
at the bench.

The Drive connector can read a sheet and create files but **cannot write cells or
add a tab**, so no skill can push into an existing Google Sheet. Paste the TSV, or
use a Sheets-API / Apps Script if you want it hands-off.

## anneal

Oligo order sheet → one `.dna` per annealed duplex, each with its own sticky
ends. Not the ligated product.

```bash
python ~/.claude/skills/anneal/scripts/anneal.py ORDER.csv -o OUTDIR \
    --report REPORT.csv --designs DESIGNS.csv
```

Pass both optional flags — they matter for correctness, not convenience:
`--report` restores oligos that dedup removed (otherwise those designs are
skipped entirely), and `--designs` splits RTT from PBS (the order sheet alone
cannot say where the boundary is).

## gga

Oligo order sheet → circular product plasmids. Digests the vector with BsmBI,
ligates by overhang, verifies, and carries the vector's own features (Cre,
mScarlet, AmpR, ori…) onto the product's coordinates.

```bash
python ~/.claude/skills/gga/scripts/gga.py ORDER.csv -o OUTDIR \
    --report REPORT.csv --designs DESIGNS.csv \
    --vector CH=/path/UPEC_tevo.dna --vector AML=/path/UPEmS_tevo.dna
```

A design that fails verification is reported and skipped, never written.

## sci-plot

Publication figure style for matplotlib/seaborn: Nature-style panels at final
print size (mm widths, 7 pt text), open axes, black-edged marks, fixed color
meanings (tissue/sample type, ABE vs CBE, COSMIC roles, enriched vs depleted),
PDF export with editable text. Tested recipes for volcano, waterfall, bar +
replicate dots, box + exact-P brackets, correlation scatter, z-score heatmap,
PCA and more, in `scripts/kdplot.py`.

Loads for any figure work — "plot the LFCs", "make figure 3c", "clean this up
for the paper" — even when style isn't mentioned.

**Every other skill in this repo defers to `sci-plot` for figures** rather than
defining its own style, so the lab's figures stay consistent and the style
lives in one place. A skill that produces a plot says so and hands off.

---

## Conventions shared by the pegRNA skills

**Colors in `.dna` output** — spacer cyan `#00ccff`, RTT amber `#ffcc00`, PBS
green `#99cc00`, **all overhangs and junctions red** `#ff0000`, scaffold grey,
tevopreQ1 mauve, cassette bisque.

**`GTTTTAGAG` is sgRNA scaffold nt 1–9, not a linker.** `GTTTT` is
double-stranded and `AGAG` is the sticky end; both are scaffold. The scaffold
feature spans all 9 nt with the red overhang drawn over its last 4. The overlap
is deliberate.

**SnapGene strandedness and sticky ends are two separate fields**, and both have
to be right or the file looks plausible and is wrong:

| | DNA flags byte | block 8 stickiness |
|---|---|---|
| annealed duplex | `0x02` linear + double-stranded | `(4, 4)` |
| product plasmid | `0x03` circular + double-stranded | `(0, 0)` |

Bit 0 alone is "circular", so `0x00` for "linear" silently drops the second
strand. Omitting block 8 defaults stickiness to 0, i.e. **blunt** — the
overhangs then exist only as red annotations and SnapGene will not ligate the
molecule. Reference file:
`~/Desktop/HSC x Order/Assembly/pegRNA_scaffold.dna` (`flags=2`, 71 nt,
`up=4 down=4`).

**PEGG's `RTT`/`PBS` orientation is already correct — do not reverse-complement
it.** Verified across every design in the CH/AML pilot:

- `RTT_PBS == RTT + PBS` (5'-RTT-PBS-3', the natural 3' extension order)
- `RTT` appears **as-is** in `sensor_alt`
- `rc(PBS)` — not `PBS` — is found inside the protospacer, ending 3 nt from its
  3' end (the PE2 nick), because PBS anneals to the nicked DNA strand

So `E_T = "gtgc" + RTT_PBS` with no rc. Sanity check on new data:
`protospacer.find(rc(PBS)) >= 0` and `protospacer.find(PBS) < 0`. If `PBS`
itself is found, the orientation is reversed.

**The three cloning cases** (backbone + epeg decide the 3' extension oligos):

| backbone | epeg | up OH | down OH | 3' ext top | 3' ext bottom |
|---|---|---|---|---|---|
| `UPEmS-tevo`, `UPEC-tevo` | *(unset)* | CACC | CGCG | `gtgc`+ext | `cgcg`+rc(ext) |
| `UPEmS`, `UPEC` | `oligo` | CACC | TTTT | `gtgc`+ext+tevo | `aaaa`+rc(tevo)+rc(ext) |
| `UPEmS`, `UPEC` | `none` | CACC | TTTT | `gtgc`+ext | `aaaa`+rc(ext) |

Impossible combinations raise rather than producing dead oligos.

---

## Where the work using these lives

- `~/Documents/GitHub/early-leukemia-evolution/20260930_ORDERED/` — the finalized
  16-mutation CH/AML pilot order (10 CH → UPEC-tevo, 6 AML → UPEmS-tevo), with
  `gga_plan/` for the 20261005 cloning and `viz_pilot/` figures.
- Cloning planning sheet: *02 Golden Gate Assembly (pegRNAs)* on Google Sheets
  (tabs by date).
- `~/Documents/GitHub/early-leukemia-evolution/archive/20260929_order_oligos/` —
  the earlier 73-design pilot, kept for reference.
- `~/Documents/GitHub/prime-editing-visualization/` — Flask pegRNA design
  visualizer, `python run.py` → <http://127.0.0.1:5050> (needs the `pegg_env`
  conda environment; see [`envs/pegg_env.yml`](../envs/pegg_env.yml)).

## Adding or changing a skill

See [`writing-a-skill.md`](writing-a-skill.md) for the full convention.

These skills are **no longer machine-local**. They live in this repo and
`install.sh` symlinks them into `~/.claude/skills/`, so editing either path
edits the same files and a `git pull` updates what Claude loads. Moving to a
new machine means cloning this repo and running `./install.sh` — nothing to
copy by hand.
