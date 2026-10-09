#!/usr/bin/env python3
"""Build the Golden Gate cloning bench sheet from an oligo order.

Reproduces the layout of "02 Golden Gate Assembly (pegRNAs)": the annealing
master mix, the per-tube annealing layout (one P and one E row per pegRNA,
grouped into racks), the Golden Gate master mix, and the GGA layout.

    python gga_plan.py REPORT.csv -o OUTDIR --start-tube 79 --title "HSC Pilot"

`--start-tube` is the index of this order's FIRST oligo in your IDT inventory;
oligos are numbered from there in order-sheet row order (P_T, P_B, E_T, E_B per
design). Omit it to label wells by oligo name instead of inventory number.

Output: gga_plan.tsv (paste into a new sheet tab), plus gga_plan.csv and the
same tables as *_anneal.csv / *_gga.csv for scripting.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import sys

# Volumes, from the lab's sheet. Per-tube µL.
ANNEAL = [("H2O (uL)", 6.5), ("T4 Buffer (uL)", 1), ("Top Oligo (uL)", 1),
          ("Bottom Oligo (uL)", 1), ("T4 Ligase", 0.5)]
ANNEAL_TOTAL = 10
GGA = [("H2O (uL)", 12.5), ("T4 Buffer (uL)", 2), ("Vector (uL)", 1.5),
       ("Scaffold (uL)", 1), ("GGA Enzyme Mix", 1), ("Annealed P (uL)", 1),
       ("Annealed E (uL)", 1)]
GGA_TOTAL = 20

# The master mix is H2O + buffer + ligase; top and bottom oligos go in per tube.
PREMIX_PARTS = {"H2O (uL)", "T4 Buffer (uL)"}
PER_RACK = 8          # annealing rows per rack block, as in the lab's sheet


def read_report(path: str) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def oligo_numbers(n_designs: int, start: int | None) -> dict[tuple[int, str, str], object]:
    """Map (design_index, 'P'|'E', 'T'|'B') -> inventory number, or None.

    Order-sheet row order is P_T, P_B, E_T, E_B per design, so the k-th design's
    four oligos are start+4k .. start+4k+3.
    """
    if start is None:
        return {}
    out = {}
    for k in range(n_designs):
        base = start + 4 * k
        out[(k, "P", "T")] = base
        out[(k, "P", "B")] = base + 1
        out[(k, "E", "T")] = base + 2
        out[(k, "E", "B")] = base + 3
    return out


def build(rows: list[dict], start_tube: int | None, title: str,
          date: str, vector_label: str) -> tuple[list[list], dict]:
    n = len(rows)
    n_anneal = 2 * n            # one P duplex and one E duplex per design
    n_oligo = 4 * n
    nums = oligo_numbers(n, start_tube)

    def cell(k: int, code: str, tag: str) -> str:
        """Inventory number when known, else the oligo's full name."""
        if nums:
            return str(nums[(k, code, tag)])
        return f"{rows[k]['name_stem']}_{code}_{tag}"

    sheet: list[list] = []
    add = sheet.append

    add([f"{date} cloning: {title}"])
    add([f"{n} GGAs, {n_anneal} annealings, {n_oligo} oligos."])
    add([])

    # --- step 1: annealing -------------------------------------------------
    # Master mix is made for the annealing tubes only; ligase is added per tube
    # on ice, so it is excluded from the premix volume.
    premix_per_tube = sum(v for k, v in ANNEAL if k in PREMIX_PARTS)
    steps = [
        "1. Anneal oligos.",
        "1.1. Make master mix of H2O and T4 buffer",
        f"1.2. Label tubes 1P, 1E, ..., {n}P, {n}E",
        f"1.3. Split into {n_anneal} tubes, {premix_per_tube:g}uL each",
        "1.4. Add top and bottom oligos to each tube",
        f"1.5. Add {dict(ANNEAL)['T4 Ligase']:g}uL Ligase to each tube",
        "1.6. Annealing Program",
        "1.7. Dilute 2:200uL and stock",
    ]
    mix_rows = [["Anneal Tube", "N = ? >", n_anneal]]
    for label, per in ANNEAL:
        mix_rows.append([label, per, round(per * n_anneal, 2)])
    mix_rows.append(["Total (uL)", ANNEAL_TOTAL, round(ANNEAL_TOTAL * n_anneal, 2)])

    for i in range(max(len(steps), len(mix_rows))):
        left = [steps[i]] if i < len(steps) else [""]
        right = mix_rows[i] if i < len(mix_rows) else ["", "", ""]
        add([""] + left + ["", "", "", "", ""] + list(right))
    add([])

    # --- annealing layout, in racks of PER_RACK rows -----------------------
    anneal_tbl: list[list] = []
    rack, row_in_rack = 1, 0
    add(["", f"TUBE {rack}", "", "Mix", "Top 1uL", "Bottom 1uL"])
    anneal_tbl.append(["rack", "anneal_tube", "part", "mix_uL", "top", "bottom", "design"])
    for k in range(n):
        for code in ("P", "E"):
            if row_in_rack == PER_RACK:
                add([])
                rack += 1
                row_in_rack = 0
                add(["", f"TUBE {rack}", "", "Mix", "Top 1uL", "Bottom 1uL"])
            add(["", k + 1, code, premix_per_tube,
                 cell(k, code, "T"), cell(k, code, "B")])
            anneal_tbl.append([rack, k + 1, code, premix_per_tube,
                               cell(k, code, "T"), cell(k, code, "B"),
                               rows[k]["name_stem"]])
            row_in_rack += 1
    add([])
    add([])

    # --- step 2: Golden Gate ----------------------------------------------
    gga_steps = ["2. Golden Gate.",
                 "2.1. Make master mix of H2O, T4 buffer, vector, scaffold, enzyme",
                 f"2.2. Split into {n} tubes",
                 "2.3. Add 1uL annealed P and 1uL annealed E to each",
                 "2.4. (42C 1min -> 16C 1min) x60 -> 60C 5min",
                 "2.5. Transform the whole 20uL"]
    multimix = sum(v for label, v in GGA
                   if label not in {"Annealed P (uL)", "Annealed E (uL)"})
    # Designs going into different vectors need SEPARATE master mixes, so the
    # mix table is per vector group rather than one column for everything.
    by_vec: dict[str, int] = {}
    for r in rows:
        by_vec[r.get("vector") or vector_label] = by_vec.get(r.get("vector") or vector_label, 0) + 1
    gmix = [["GG Tube", "N = ? >"] + [f"{v} (n={c})" for v, c in by_vec.items()]]
    for label, per in GGA:
        name = vector_label if label == "Vector (uL)" else label
        gmix.append([name, per] + [round(per * c, 2) for c in by_vec.values()])
    gmix.append(["Total (uL)", GGA_TOTAL] + [round(GGA_TOTAL * c, 2) for c in by_vec.values()])

    for i in range(max(len(gga_steps), len(gmix))):
        left = [gga_steps[i]] if i < len(gga_steps) else [""]
        right = gmix[i] if i < len(gmix) else [""] * (2 + len(by_vec))
        add([""] + left + ["", "", "", "", ""] + list(right))
    add([])

    gga_tbl = [["gga_tube", "group", "vector", "multimix_uL",
                "annealed_P", "annealed_E", "design"]]
    add(["", "GG Tube", "Group", "Multimix", "Annealed P 1uL", "Annealed E 1uL"])
    for k, r in enumerate(rows):
        grp = r.get("Group", "")
        add(["", k + 1, grp, multimix, f"{k + 1}P", f"{k + 1}E", r["name_stem"]])
        gga_tbl.append([k + 1, grp, r.get("vector", vector_label), multimix,
                        f"{k + 1}P", f"{k + 1}E", r["name_stem"]])

    return sheet, {"anneal": anneal_tbl, "gga": gga_tbl,
                   "n_designs": n, "n_anneal": n_anneal, "n_oligo": n_oligo}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", help="*_report.csv from the order (Group, name_stem, ...)")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--start-tube", type=int, default=None,
                    help="inventory index of this order's FIRST oligo (e.g. 79)")
    ap.add_argument("--title", default="pegRNA cloning")
    ap.add_argument("--date", default=_dt.date.today().strftime("%Y%m%d"))
    ap.add_argument("--vector", default="Vector (uL)",
                    help="label for the vector row, e.g. 'UPEmS-tevo (uL)'")
    ap.add_argument("--name", default="gga_plan")
    args = ap.parse_args()

    rows = read_report(args.report)
    if not rows:
        print("report is empty", file=sys.stderr)
        return 1
    if "name_stem" not in rows[0]:
        print("report has no name_stem column", file=sys.stderr)
        return 1

    sheet, tbl = build(rows, args.start_tube, args.title, args.date, args.vector)
    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir, args.name)

    with open(f"{base}.tsv", "w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(sheet)
    with open(f"{base}.csv", "w", newline="") as handle:
        csv.writer(handle).writerows(sheet)
    for key in ("anneal", "gga"):
        with open(f"{base}_{key}.csv", "w", newline="") as handle:
            csv.writer(handle).writerows(tbl[key])

    groups = {}
    for r in rows:
        groups[r.get("Group", "")] = groups.get(r.get("Group", ""), 0) + 1
    print(f"{tbl['n_designs']} GGAs, {tbl['n_anneal']} annealings, "
          f"{tbl['n_oligo']} oligos")
    for g, c in groups.items():
        print(f"  {g or 'all'}: {c}")
    if args.start_tube:
        print(f"inventory tubes {args.start_tube}-{args.start_tube + tbl['n_oligo'] - 1}")
    print(f"\nwrote {base}.tsv   <- select all, copy, paste into a new sheet tab")
    print(f"      {base}.csv, {base}_anneal.csv, {base}_gga.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
