#!/usr/bin/env python3
"""Build annotated SnapGene .dna files from a pegRNA oligo order sheet.

Each *annealed oligo duplex* becomes its own .dna file, with its 4-nt sticky
ends shown — not the ligated product. A design yields

    1_Asxl1_G643WfsX14_P.dna    protospacer duplex   CACC | spacer | AGAG
    1_Asxl1_G643WfsX14_E.dna    extension duplex     GTGC | RTT | PBS | CGCG

plus `_S.dna` when scaffold oligos are in the sheet. spacer / PBS / RTT and the
overhangs are annotated in color.

    python anneal.py ORDER.csv -o OUTDIR [--designs DESIGNS.csv] [--report REPORT.csv]

ORDER.csv is the IDT sheet written by the pegRNA GGA notebook (Name, Sequence,
Scale, Purification). Names carry the oligo role, e.g.

    1_Asxl1_G643WfsX14_P_T   protospacer top        (cacc + spacer + gtttt)
    1_Asxl1_G643WfsX14_P_B   protospacer bottom
    1_Asxl1_G643WfsX14_E_T   3' extension top       (gtgc + RTT+PBS [+ tevopreQ1])
    1_Asxl1_G643WfsX14_E_B   3' extension bottom

Deduplicated oligos are shared between designs; pass --report so the reused
names can be resolved back to the tube that was actually ordered. Pass
--designs to split the 3' extension into its true RTT and PBS parts (the order
sheet alone cannot tell where one ends and the other begins).
"""

from __future__ import annotations

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from snapgene_writer import Feature, read_dna, write_dna  # noqa: E402

# --- fixed cassette parts (must match the notebook) --------------------------

SCAFFOLD = "CTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG"
TEVOPREQ1 = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"
_COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# --- colors ------------------------------------------------------------------
# Distinct hues for the three functional parts; overhangs red, as asked.
COLORS = {
    "spacer": "#00ccff",    # cyan
    "PBS": "#99cc00",       # green
    "RTT": "#ffcc00",       # amber
    "overhang": "#ff0000",  # red
    "scaffold": "#a6acb3",  # grey — context, not one of the three
    "tevopreQ1": "#cc99b2", # mauve — context
}


def parse_order(path: str) -> dict[str, dict[str, str]]:
    """Group the order sheet by design stem -> {'P_T': seq, 'P_B': ..., ...}."""
    import csv

    designs: dict[str, dict[str, str]] = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            name, seq = row["Name"].strip(), row["Sequence"].strip()
            m = re.match(r"^(.*)_([PSE])_([TB])$", name)
            if not m:
                continue
            stem, code, tag = m.groups()
            designs.setdefault(stem, {})[f"{code}_{tag}"] = seq
    return designs


def resolve_reused(designs: dict[str, dict[str, str]], report_path: str) -> None:
    """Fill in oligos that were deduped away, using the report's `reused` column.

    Without this, a design that shares a protospacer with an earlier one is
    missing its P_T/P_B in the order sheet and cannot be assembled.
    """
    import csv

    with open(report_path, newline="") as handle:
        for row in csv.DictReader(handle):
            reused = (row.get("reused") or "").strip()
            if not reused:
                continue
            for mapping in reused.split(","):
                if "->" not in mapping:
                    continue
                missing, ordered = (s.strip() for s in mapping.split("->"))
                m_miss = re.match(r"^(.*)_([PSE])_([TB])$", missing)
                m_ord = re.match(r"^(.*)_([PSE])_([TB])$", ordered)
                if not (m_miss and m_ord):
                    continue
                stem, code, tag = m_miss.groups()
                src_stem, src_code, src_tag = m_ord.groups()
                src = designs.get(src_stem, {}).get(f"{src_code}_{src_tag}")
                if src:
                    designs.setdefault(stem, {}).setdefault(f"{code}_{tag}", src)


def load_design_parts(path: str) -> dict[str, dict[str, str]]:
    """Map design_id -> {'RTT':..., 'PBS':...} from the designs csv."""
    import csv

    parts: dict[str, dict[str, str]] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"RTT", "PBS"} <= set(reader.fieldnames or []):
            return parts
        for row in reader:
            key = (row.get("pegRNA_id") or "").strip()
            if key:
                parts[key] = {"RTT": row["RTT"].strip().upper(),
                              "PBS": row["PBS"].strip().upper()}
    return parts


def duplex_sequence(top: str, bot: str) -> tuple[str, int, int]:
    """Lay out one annealed oligo pair as a single span.

    Both oligos carry a 4-nt 5' overhang, so the duplex is

        5'-NNNN[----- core -----]      -3'   top
        3'      [----- core -----]NNNN-5'   bottom

    The returned sequence is left overhang + core + right overhang, i.e. the top
    strand plus the 4 nt that only the bottom strand supplies. The core is the
    double-stranded part; the flanks are the sticky ends.
    """
    top, bot = top.strip(), bot.strip()
    if len(top) != len(bot):
        raise ValueError(f"top ({len(top)} nt) and bottom ({len(bot)} nt) differ")
    core = top[4:]
    if rc(bot)[:-4] != core:
        raise ValueError("oligos do not anneal — core is not complementary")
    right = rc(bot)[-4:]
    return (top + right).upper(), 4, len(top)  # seq, core_start-1, core_end


def annotate_protospacer(oligos: dict[str, str]) -> tuple[str, list[Feature]]:
    """The P duplex: CACC overhang | spacer | GTTTTAGAG (scaffold 5' end).

    The GTTTTAGAG tail is scaffold nt 1-9, not a linker: GTTTT is the
    double-stranded part and AGAG is the 4-nt sticky end. Both are scaffold, so
    the whole motif is annotated as such and the overhang is drawn over its last
    4 nt — a scaffold feature that stopped at GTTTT would imply AGAG is
    something else.
    """
    top, bot = oligos["P_T"], oligos["P_B"]
    seq, _, core_end = duplex_sequence(top, bot)
    spacer_start, spacer_end = 5, len(top) - 5
    tail = seq[spacer_end:]  # GTTTTAGAG
    feats = [
        Feature("CACC overhang", 1, 4, COLORS["overhang"]),
        Feature("spacer", spacer_start, spacer_end, COLORS["spacer"], directionality=1),
        Feature(f"sgRNA scaffold 5' end ({tail})", spacer_end + 1, len(seq),
                COLORS["scaffold"]),
        Feature("AGAG overhang (scaffold nt 6-9)", core_end + 1, len(seq),
                COLORS["overhang"]),
    ]
    return seq, feats


def annotate_extension(oligos: dict[str, str], parts: dict[str, str] | None) -> tuple[str, list[Feature]]:
    """The E duplex: GTGC | RTT | PBS | [tevopreQ1] | CGCG-or-AAAA."""
    top, bot = oligos["E_T"], oligos["E_B"]
    seq, _, core_end = duplex_sequence(top, bot)
    body = top[4:].upper()

    tevo_len = 0
    if body.endswith(TEVOPREQ1):
        tevo_len = len(TEVOPREQ1)
        body = body[:-tevo_len]

    feats = [Feature("GTGC overhang", 1, 4, COLORS["overhang"])]
    at = 5  # first base of the duplex core
    placed = False
    if parts:
        rtt, pbs = parts.get("RTT", ""), parts.get("PBS", "")
        if rtt and pbs and body == rtt + pbs:
            feats.append(Feature("RTT", at, at + len(rtt) - 1, COLORS["RTT"], directionality=1))
            feats.append(Feature("PBS", at + len(rtt), at + len(rtt) + len(pbs) - 1,
                                 COLORS["PBS"], directionality=1))
            placed = True
    if not placed:
        # Without the designs csv the RTT/PBS boundary is unknowable; one block
        # beats an invented split.
        feats.append(Feature("3' extension (RTT+PBS)", at, at + len(body) - 1,
                             COLORS["RTT"], directionality=1))
    if tevo_len:
        feats.append(Feature("tevopreQ1", at + len(body), at + len(body) + tevo_len - 1,
                             COLORS["tevopreQ1"]))

    right = seq[core_end:]
    feats.append(Feature(f"{right} overhang", core_end + 1, len(seq), COLORS["overhang"]))
    return seq, feats


def annotate_scaffold(oligos: dict[str, str]) -> tuple[str, list[Feature]]:
    """The optional S duplex: AGAG | scaffold | GTGC overhang."""
    top, bot = oligos["S_T"], oligos["S_B"]
    seq, _, core_end = duplex_sequence(top, bot)
    feats = [
        Feature("AGAG overhang", 1, 4, COLORS["overhang"]),
        Feature("sgRNA scaffold", 5, core_end, COLORS["scaffold"]),
        Feature(f"{seq[core_end:]} overhang", core_end + 1, len(seq), COLORS["overhang"]),
    ]
    return seq, feats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("order", help="IDT order csv (Name,Sequence,...)")
    ap.add_argument("-o", "--outdir", required=True, help="output directory for .dna files")
    ap.add_argument("--report", help="report csv, to resolve deduplicated oligos")
    ap.add_argument("--designs", help="designs csv, to split RTT from PBS")
    ap.add_argument("--only", help="regex; build only designs whose name matches")
    args = ap.parse_args()

    designs = parse_order(args.order)
    if args.report and os.path.exists(args.report):
        resolve_reused(designs, args.report)

    parts_by_id = load_design_parts(args.designs) if args.designs else {}
    # The order sheet's stem is `index_gene_mut`; the designs csv is keyed by
    # pegRNA_id, so match on RTT+PBS content instead of trusting either name.
    parts_by_seq = {v["RTT"] + v["PBS"]: v for v in parts_by_id.values()}

    os.makedirs(args.outdir, exist_ok=True)
    written, skipped = [], []
    for stem in sorted(designs, key=lambda s: (int(m.group(1)) if (m := re.match(r"^(\d+)", s)) else 0, s)):
        if args.only and not re.search(args.only, stem):
            continue
        oligos = designs[stem]

        # The RTT/PBS split is looked up by extension content, not by name.
        body = oligos.get("E_T", "")[4:].upper()
        if body.endswith(TEVOPREQ1):
            body = body[: -len(TEVOPREQ1)]
        parts = parts_by_seq.get(body)

        # One file per annealed duplex, each with its own sticky ends.
        for code, builder in (("P", lambda o: annotate_protospacer(o)),
                              ("S", lambda o: annotate_scaffold(o)),
                              ("E", lambda o: annotate_extension(o, parts))):
            if f"{code}_T" not in oligos or f"{code}_B" not in oligos:
                continue  # scaffold is usually not ordered; P/B may be deduped away
            name = f"{stem}_{code}"
            try:
                seq, feats = builder(oligos)
                # Both ends carry a real 4-nt 5' overhang; block 8 is what makes
                # them single-stranded rather than just red annotations.
                write_dna(os.path.join(args.outdir, f"{name}.dna"), seq, feats,
                          up_sticky=4, down_sticky=4)
                written.append((name, len(seq), len(feats)))
            except Exception as exc:  # noqa: BLE001 - report and keep going
                skipped.append((name, str(exc)))

    for name, length, nfeat in written:
        print(f"  {name}.dna  {length} bp, {nfeat} features")
    print(f"\nwrote {len(written)} .dna files to {args.outdir}")
    if skipped:
        print(f"skipped {len(skipped)}:")
        for name, why in skipped:
            print(f"  {name}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
