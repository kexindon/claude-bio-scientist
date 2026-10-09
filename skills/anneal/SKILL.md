---
name: anneal
description: Build annotated SnapGene .dna files from a pegRNA oligo order sheet — one file per annealed oligo duplex (_P protospacer, _E extension, _S scaffold), each showing its own 4-nt sticky ends, with spacer, PBS, RTT and the overhangs annotated in color (overhangs red). Use whenever asked to "anneal" oligos, turn an oligo order csv into .dna or SnapGene files, visualize or annotate ordered pegRNA oligos, or make maps of the annealed duplexes before Golden Gate. Works on the order sheets written by the pegRNA GGA notebook (pegrna_gga_oligo_order.ipynb, columns Name/Sequence/Scale/Purification with names like 1_Asxl1_G643WfsX14_P_T).
---

# anneal — oligo order sheet → annotated SnapGene maps

Turns an IDT pegRNA oligo order into one `.dna` file **per annealed oligo
duplex**, each carrying its own sticky ends, with each functional part colored.
These are the separate annealed products that go into the Golden Gate reaction —
**not** the ligated insert.

A design yields two files (three when scaffold oligos are in the sheet):

```
1_Asxl1_G643WfsX14_P.dna    CACC | spacer | GTTTTAGAG        33 bp
1_Asxl1_G643WfsX14_E.dna    GTGC | RTT | PBS | CGCG          37 bp
1_Asxl1_G643WfsX14_S.dna    AGAG | scaffold | GTGC           (if ordered)
```

## Usage

```bash
python ~/.claude/skills/anneal/scripts/anneal.py ORDER.csv -o OUTDIR \
    --report REPORT.csv --designs DESIGNS.csv
```

| flag | why it matters |
|---|---|
| `ORDER.csv` | the IDT sheet (`*_oligos.csv`): `Name,Sequence,Scale,Purification` |
| `-o OUTDIR` | created if absent; one `<stem>.dna` per design |
| `--report` | **pass this** — resolves deduplicated oligos (see below) |
| `--designs` | **pass this** — splits the 3' extension into true RTT and PBS |
| `--only REGEX` | build a subset, e.g. `--only 'Trp53'` |

Both optional flags matter for correctness, not convenience:

- **`--report`**: designs sharing a protospacer have their `P_T`/`P_B` ordered
  only once, so those oligos are *absent* from the order sheet for the second
  design. The report's `reused` column (`29_..._P_T->28_..._P_T`) maps them back.
  Without it, every deduplicated design is skipped.
- **`--designs`**: the order sheet concatenates RTT+PBS into one extension
  oligo and cannot say where the boundary is. With the designs csv the two are
  annotated separately; without it they collapse into a single
  `3' extension (RTT+PBS)` block. The script matches on RTT+PBS sequence
  content, not on names, and only splits when `RTT + PBS` reproduces the
  extension exactly — so a mismatched csv degrades to the single block rather
  than drawing a wrong boundary.

## Colors

| part | color | |
|---|---|---|
| spacer | `#00ccff` | cyan |
| RTT | `#ffcc00` | amber |
| PBS | `#99cc00` | green |
| **all overhangs / junctions** | `#ff0000` | **red** |
| sgRNA scaffold | `#a6acb3` | grey (context) |
| tevopreQ1 | `#cc99b2` | mauve (context) |

Change them in `COLORS` at the top of `scripts/anneal.py`.

## What each file contains

One annealed duplex, written as `top strand + the 4 nt only the bottom strand
supplies`. Both oligos carry a 4-nt 5' overhang, so every file is
`sticky end | double-stranded core | sticky end`:

```
5'-NNNN[------- core -------]      -3'   top oligo
3'      [------- core -------]NNNN-5'   bottom oligo
   ^^^^                       ^^^^
   red                        red
```

Files are **linear and double-stranded**, so SnapGene shows both strands and the
bottom strand reads back as the ordered `_B` oligo. `_P` is 33 bp with 4
features; `_E` runs 37–95 bp depending on RTT/PBS length and whether tevopreQ1
is on the oligo.

```{warning}
Two separate things must both be right, or the file looks plausible and is wrong:

1. **The DNA block's flags byte needs bit 1 set** (`0x02`) or SnapGene opens the
   file as **single-stranded**. Bit 0 alone is "circular", so a linear molecule is
   `0x02`, not `0x00`. A real plasmid reads `31` (`0b11111`: circular +
   double-stranded + three methylation bits).
2. **Block 8 must declare the sticky ends.** `UpstreamStickiness` /
   `DownstreamStickiness` hold the number of unpaired 5' nt at each end — 4 for
   these duplexes. With block 8 absent both default to `0`, i.e. a **blunt**
   fragment: the overhangs then exist only as red annotations, the strands are
   drawn flush, and SnapGene will not ligate the molecule.

`~/Desktop/HSC x Order/Assembly/pegRNA_scaffold.dna` is the reference: `flags=2`,
71 nt, `up=4 down=4`, with red `Overhang` features at 1-4 and 68-71.
```

### Three geometry details the annotation gets right

- **The right sticky end lives only on the bottom oligo.** The top strand stops
  at the end of the core, so those 4 nt must be appended to the sequence. A
  naive `end-3..end` span over the top strand alone would recolor the last 4 nt
  of the PBS red instead.
- **`GTTTTAGAG` in the `_P` duplex is scaffold, all 9 nt.** It is sgRNA scaffold
  nt 1–9, not a linker: `GTTTT` is double-stranded (on both oligos) and `AGAG`
  is the sticky end. The scaffold feature therefore spans the **whole**
  `GTTTTAGAG` and the red `AGAG overhang` is drawn **over its last 4 nt**. A
  scaffold feature that stopped at `GTTTT` would wrongly imply `AGAG` is
  something other than scaffold. The two features overlap deliberately;
  SnapGene and Biopython both render that fine.
- **The `_E` right overhang is read from the bottom oligo**, so it self-labels
  `CGCG` for `-tevo` backbones and `AAAA` for plain UPEmS/UPEC rather than being
  hard-coded.

## Verifying output

`snapgene_writer.read_dna(path)` returns `(sequence, [(name, start, end, color)])`
for round-trip checks. Biopython also reads the result
(`SeqIO.read(path, "snapgene")`), which is a good independent check that the
binary is valid. Worth asserting after a run:

- the DNA block's flags byte is `2` (linear + double-stranded)
- `read_stickiness(path)` returns `(4, 4)` — **not** `(0, 0)`, which would mean
  blunt ends dressed up as red annotations
- the single-stranded regions equal the real overhangs: `seq[:4] == _T[:4]` and
  `seq[-4:] == rc(_B[:4])`
- each file's sequence equals `top + rc(bottom)[-4:]` for its oligo pair, and
  `rc(sequence)` equals `bottom + rc(top)[-4:]` — i.e. **both** ordered oligos
  are recoverable from the duplex, which is the real test that the file
  represents the annealed molecule
- the `spacer` feature equals the report's `spacer` column
- `RTT` and `PBS` equal the designs csv columns
- exactly two red 4-nt overhang features per file, one at each end
- no overhang overlaps `spacer`, `RTT` or `PBS` — but the `AGAG overhang` **is**
  expected to sit inside the `sgRNA scaffold 5' end` feature
- every feature ends within the sequence

When filtering features by name, match `startswith("sgRNA scaffold")` rather
than `"scaffold" in name`: the AGAG overhang's label also contains the word.

## Why a custom writer

Biopython **reads** SnapGene but cannot write it — `snapgene` is not in
`SeqIO._FormatToWriter`. `scripts/snapgene_writer.py` emits the binary directly:
a block stream of `[1-byte type][4-byte big-endian length][body]`, of which only
three blocks are needed — `9` (cookie), `0` (DNA, 1 flags byte + lowercase
sequence), `10` (features XML). Per-segment `color="#rrggbb"` is what drives the
display color. Writing GenBank instead would lose the colors, which are the
point here.

## Scope

Built for the `P`/`E` (+ optional `S`) oligo naming of
`pegrna_gga_oligo_order.ipynb`. A design needs at least `P_T` and `E_T` to be
built; anything unbuildable is reported at the end with a reason rather than
failing the whole run.
