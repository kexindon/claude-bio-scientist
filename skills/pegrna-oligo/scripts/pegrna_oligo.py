#!/usr/bin/env python3
"""pegRNA designs -> IDT order sheet + annealed duplex maps + product plasmids.

One command covers the whole path from a PEGG/PRIDICT design table to everything
you need at the bench:

    python pegrna_oligo.py DESIGNS.csv -o OUTDIR \
        --group-col gene_group \
        --group CH=UPEC-tevo:/path/UPEC_tevo.dna \
        --group AML=UPEmS-tevo:/path/UPEmS_tevo.dna

produces, under OUTDIR:

    pegrna_orders/<name>_oligos.csv     IDT upload (Name,Sequence,Scale,Purification)
    pegrna_orders/<name>_report.csv     one row per design + Group + QC
    pegrna_orders/<name>_designs.csv    the designs actually ordered
    oligos/<stem>_P.dna, _E.dna         annealed duplexes, real 4-nt sticky ends
    plasmids/<stem>.dna                 circular assembly products

Each group is checked against its own backbone and vector map but numbered and
exported as ONE order, so CH and AML share a continuous index and a single sheet.

Use --no-anneal / --no-plasmids to stop after the order sheet.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pandas as pd  # noqa: E402

import pegoligo as P  # noqa: E402
from snapgene_writer import Feature, write_dna  # noqa: E402

COLORS = {
    "spacer": "#00ccff",
    "PBS": "#99cc00",
    "RTT": "#ffcc00",
    "overhang": "#ff0000",
    "scaffold": "#a6acb3",
    "tevopreQ1": "#cc99b2",
    "cassette": "#ffe4c4",
}


# --------------------------------------------------------------------------
# annealed duplexes
# --------------------------------------------------------------------------

def duplex(top: str, bot: str) -> tuple[str, int]:
    """Annealed pair as `left overhang + core + right overhang`.

    Both oligos carry a 4-nt 5' overhang, so the stored sequence is the top
    strand plus the 4 nt only the bottom strand supplies.
    """
    if len(top) != len(bot):
        raise ValueError(f"top {len(top)} nt vs bottom {len(bot)} nt")
    core = top[4:]
    if P.rc(bot)[:-4] != core:
        raise ValueError("oligos do not anneal — core is not complementary")
    return (top + P.rc(bot)[-4:]).upper(), len(top)


def write_duplexes(outdir: str, stem: str, oligos: dict[str, str],
                   rtt: str, pbs: str) -> int:
    """Write <stem>_P.dna and <stem>_E.dna. Returns how many were written."""
    os.makedirs(outdir, exist_ok=True)
    written = 0

    if "P_T" in oligos and "P_B" in oligos:
        seq, core_end = duplex(oligos["P_T"], oligos["P_B"])
        spacer_end = len(oligos["P_T"]) - 5
        # GTTTTAGAG is sgRNA scaffold nt 1-9, not a linker: GTTTT is
        # double-stranded, AGAG is the sticky end. Annotate the whole motif as
        # scaffold and draw the overhang over its last 4 nt.
        feats = [
            Feature("CACC overhang", 1, 4, COLORS["overhang"]),
            Feature("spacer", 5, spacer_end, COLORS["spacer"], directionality=1),
            Feature(f"sgRNA scaffold 5' end ({seq[spacer_end:]})",
                    spacer_end + 1, len(seq), COLORS["scaffold"]),
            Feature("AGAG overhang (scaffold nt 6-9)", core_end + 1, len(seq),
                    COLORS["overhang"]),
        ]
        write_dna(os.path.join(outdir, f"{stem}_P.dna"), seq, feats,
                  up_sticky=4, down_sticky=4)
        written += 1

    if "E_T" in oligos and "E_B" in oligos:
        seq, core_end = duplex(oligos["E_T"], oligos["E_B"])
        body = oligos["E_T"][4:].upper()
        tevo = len(P.TEVOPREQ1) if body.endswith(P.TEVOPREQ1) else 0
        if tevo:
            body = body[:-tevo]
        feats = [Feature("GTGC overhang", 1, 4, COLORS["overhang"])]
        at = 5
        if rtt and pbs and body == (rtt + pbs).upper():
            feats.append(Feature("RTT", at, at + len(rtt) - 1, COLORS["RTT"],
                                 directionality=1))
            feats.append(Feature("PBS", at + len(rtt),
                                 at + len(rtt) + len(pbs) - 1, COLORS["PBS"],
                                 directionality=1))
        else:
            # Without RTT/PBS columns the boundary is unknowable; one block
            # beats an invented split.
            feats.append(Feature("3' extension (RTT+PBS)", at,
                                 at + len(body) - 1, COLORS["RTT"], directionality=1))
        if tevo:
            feats.append(Feature("tevopreQ1", at + len(body),
                                 at + len(body) + tevo - 1, COLORS["tevopreQ1"]))
        feats.append(Feature(f"{seq[core_end:]} overhang", core_end + 1, len(seq),
                             COLORS["overhang"]))
        write_dna(os.path.join(outdir, f"{stem}_E.dna"), seq, feats,
                  up_sticky=4, down_sticky=4)
        written += 1
    return written


# --------------------------------------------------------------------------
# product plasmids
# --------------------------------------------------------------------------

def shift_features(record, cut_offset: int, backbone_len: int) -> list[Feature]:
    """Re-express the vector's own features on the product's coordinates."""
    orig_len = len(record.seq)
    out: list[Feature] = []
    for feat in record.features:
        label = (feat.qualifiers.get("label") or feat.qualifiers.get("note")
                 or ["feature"])[0]
        color = (feat.qualifiers.get("ApEinfo_fwdcolor")
                 or feat.qualifiers.get("ApEinfo_revcolor") or ["#a6acb3"])[0]
        if not (isinstance(color, str) and color.startswith("#") and len(color) == 7):
            color = "#a6acb3"
        segments, ok = [], True
        for part in feat.location.parts:
            start, end = int(part.start) - cut_offset, int(part.end) - cut_offset
            if start < 0:
                start += orig_len
                end += orig_len
            if end > backbone_len:          # lies in the excised stuffer
                ok = False
                break
            segments.append((start + 1, end))
        if not ok or not segments:
            continue
        # A reverse-strand joined feature lists its parts in descending order,
        # so sort before taking the outer bounds or start > end.
        segments.sort()
        strand = feat.location.strand
        out.append(Feature(label, segments[0][0], segments[-1][1], color,
                           type=str(feat.type),
                           directionality=2 if strand == -1 else (1 if strand == 1 else 0),
                           segments=segments if len(segments) > 1 else None))
    return out


def write_plasmid(outdir: str, stem: str, vector_path: str,
                  spacer: str, ext: str, rtt: str, pbs: str) -> int:
    from Bio import SeqIO

    os.makedirs(outdir, exist_ok=True)
    record = SeqIO.read(vector_path, "snapgene")
    vec = str(record.seq).upper()
    backbone_top, up, down = P.digest_bsmbi(vec)
    cut_offset = vec.find("CGTCTC") + 7

    insert = ("CACC" + spacer + "GTTTT" + "AGAG" + P.SCAFFOLD + "GTGC" + ext).upper()
    product = backbone_top + insert

    cassette = ("CACC" + spacer + "GTTTTAGAG" + P.SCAFFOLD + "GTGC" + ext).upper()
    if cassette not in product:
        raise ValueError("assembled cassette does not match the design")
    if P.bsmbi_sites(product + product[:5]):
        raise ValueError("BsmBI site survives in the product")

    feats = shift_features(record, cut_offset, len(backbone_top))
    at = len(backbone_top)
    spacer_start = at + 5
    spacer_end = spacer_start + len(spacer) - 1
    scaf_end = at + len("CACC" + spacer + "GTTTT") + len("AGAG" + P.SCAFFOLD)
    ext_start = scaf_end + 5
    body = ext.upper()
    tevo = len(P.TEVOPREQ1) if body.endswith(P.TEVOPREQ1) else 0
    if tevo:
        body = body[:-tevo]

    feats += [
        Feature("pegRNA cassette", at + 1, len(product), COLORS["cassette"]),
        Feature("CACC junction", at + 1, at + 4, COLORS["overhang"]),
        Feature("spacer", spacer_start, spacer_end, COLORS["spacer"], directionality=1),
        Feature("sgRNA scaffold", spacer_end + 1, scaf_end, COLORS["scaffold"]),
        Feature("AGAG junction", spacer_end + 6, spacer_end + 9, COLORS["overhang"]),
        Feature("GTGC junction", scaf_end + 1, scaf_end + 4, COLORS["overhang"]),
    ]
    if rtt and pbs and body == (rtt + pbs).upper():
        feats.append(Feature("RTT", ext_start, ext_start + len(rtt) - 1,
                             COLORS["RTT"], directionality=1))
        feats.append(Feature("PBS", ext_start + len(rtt),
                             ext_start + len(rtt) + len(pbs) - 1, COLORS["PBS"],
                             directionality=1))
    else:
        feats.append(Feature("3' extension (RTT+PBS)", ext_start,
                             ext_start + len(body) - 1, COLORS["RTT"], directionality=1))
    if tevo:
        feats.append(Feature("tevopreQ1", ext_start + len(body),
                             ext_start + len(body) + tevo - 1, COLORS["tevopreQ1"]))

    write_dna(os.path.join(outdir, f"{stem}.dna"), product, feats, circular=True)
    return len(product)


# --------------------------------------------------------------------------

def parse_group(spec: str) -> tuple[str, str, str]:
    """`CH=UPEC-tevo:/path/map.dna`, or `UPEC-tevo:/path/map.dna` for one group."""
    group, _, rest = spec.partition("=")
    if not rest:
        group, rest = "", spec
    backbone, _, vector = rest.partition(":")
    return group, backbone, vector


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", help="PEGG/PRIDICT design csv")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--name", default="", help="oligo name prefix; '' -> index_gene_mut_P_T")
    ap.add_argument("--file-name", default="pegrna_order", help="stem for the csv files")
    ap.add_argument("--group", action="append", required=True,
                    help="GROUP=BACKBONE:/path/map.dna (repeatable)")
    ap.add_argument("--group-col", default=None, help="column holding the group label")
    ap.add_argument("--order", default=None, help="comma-separated group order for numbering")
    ap.add_argument("--columns", default="pegg", choices=["pegg", "pridict"])
    ap.add_argument("--gene-col", default="Hugo_Symbol")
    ap.add_argument("--mut-col", default="HGVSp_m")
    ap.add_argument("--epeg", default=None, help="'oligo'|'none' for plain backbones")
    ap.add_argument("--plate", action="store_true", help="add Well Position")
    ap.add_argument("--no-anneal", action="store_true")
    ap.add_argument("--no-plasmids", action="store_true")
    ap.add_argument("--gga-plan", action="store_true",
                    help="also write the Golden Gate bench sheet (gga_plan.tsv)")
    ap.add_argument("--start-tube", type=int, default=None,
                    help="with --gga-plan: inventory index of this order's FIRST oligo")
    ap.add_argument("--plan-title", default=None,
                    help="with --gga-plan: title line for the sheet")
    args = ap.parse_args()

    groups = {}
    for spec in args.group:
        g, backbone, vector = parse_group(spec)
        if vector and not os.path.exists(vector):
            print(f"vector map not found: {vector}", file=sys.stderr)
            return 1
        groups[g] = {"backbone": backbone, "vector": vector}

    designs = pd.read_csv(args.designs)
    gcol = args.group_col
    if gcol and gcol not in designs.columns:
        print(f"--group-col {gcol!r} is not a column of the input", file=sys.stderr)
        return 1
    if gcol:
        labels = set(designs[gcol].dropna().unique())
        missing = labels - set(groups)
        if missing:
            print(f"no --group given for {sorted(missing)}", file=sys.stderr)
            return 1
        order = ([g for g in (args.order.split(",") if args.order else []) if g in labels]
                 or sorted(labels))
        order += sorted(labels - set(order))
        designs = (designs.assign(_g=pd.Categorical(designs[gcol], order, ordered=True))
                          .sort_values("_g", kind="stable").drop(columns="_g")
                          .reset_index(drop=True))
    else:
        order = [""]

    # --- order sheet: per group, one continuous numbering, shared dedup ------
    order_parts, report_parts = [], []
    seen: dict[str, str] = {}
    for grp in order:
        cfg = groups[grp] if gcol else groups[list(groups)[0]]
        subset = designs[designs[gcol] == grp] if gcol else designs
        o, r = P.build_oligos(
            subset, name=args.name, backbone=cfg["backbone"], epeg=args.epeg,
            columns=args.columns, gene_col=args.gene_col, mut_col=args.mut_col,
            include_scaffold=False, g_start="prepend", dedup=True, plate=False,
            strict=True, start_index=sum(len(p) for p in report_parts) + 1, seen=seen,
        )
        r.insert(0, "Group", grp)
        order_parts.append(o)
        report_parts.append(r)
        print(f"  {grp or 'all':>4}: {len(subset):>3} pegRNAs -> {len(o):>3} oligos "
              f"into {cfg['backbone']}")

    orders = pd.concat(order_parts, ignore_index=True)
    report = pd.concat(report_parts, ignore_index=True)
    if args.plate:
        orders.insert(0, "Well Position", P.well_positions(len(orders)))

    orders_dir = os.path.join(args.outdir, "pegrna_orders")
    os.makedirs(orders_dir, exist_ok=True)
    base = os.path.join(orders_dir, args.file_name)
    idt_cols = [c for c in ["Well Position", "Name", "Sequence", "Scale", "Purification"]
                if c in orders.columns]
    orders[idt_cols].to_csv(f"{base}_oligos.csv", index=False)
    report.to_csv(f"{base}_report.csv", index=False)
    designs.to_csv(f"{base}_designs.csv", index=False)
    print(f"\n{len(orders)} oligos -> {base}_oligos.csv  ({', '.join(idt_cols)})")

    # --- annealed duplexes and product plasmids -----------------------------
    by_name = dict(zip(orders["Name"], orders["Sequence"]))
    for row in report.itertuples():                    # restore deduped oligos
        for mapping in (getattr(row, "reused", "") or "").split(","):
            if "->" in mapping:
                miss, src = (s.strip() for s in mapping.split("->"))
                by_name.setdefault(miss, by_name[src])

    parts = {}
    if args.columns == "pegg" and {"RTT", "PBS"} <= set(designs.columns):
        id_col = "pegRNA_id" if "pegRNA_id" in designs.columns else None
        if id_col:
            parts = {r[id_col]: (str(r["RTT"]).upper(), str(r["PBS"]).upper())
                     for _, r in designs.iterrows()}

    n_dup = n_plas = 0
    for row in report.itertuples():
        stem = row.name_stem
        rtt, pbs = parts.get(row.design_id, ("", ""))
        oligos = {k: by_name[f"{stem}_{k}"] for k in ("P_T", "P_B", "E_T", "E_B")
                  if f"{stem}_{k}" in by_name}
        if not args.no_anneal:
            n_dup += write_duplexes(os.path.join(args.outdir, "oligos"), stem,
                                    oligos, rtt, pbs)
        if not args.no_plasmids:
            cfg = groups[row.Group] if gcol else groups[list(groups)[0]]
            if cfg["vector"]:
                write_plasmid(os.path.join(args.outdir, "plasmids"), stem,
                              cfg["vector"], row.spacer.upper(),
                              row.extension.upper(), rtt, pbs)
                n_plas += 1
    if n_dup:
        print(f"{n_dup} annealed duplex .dna -> {os.path.join(args.outdir, 'oligos')}")
    if n_plas:
        print(f"{n_plas} product plasmid .dna -> {os.path.join(args.outdir, 'plasmids')}")

    if args.gga_plan:
        import gga_plan as GP
        rows = report.to_dict("records")
        for r in rows:                       # let the plan name each tube's vector
            cfg = groups[r.get("Group", "")] if gcol else groups[list(groups)[0]]
            r["vector"] = cfg["backbone"]
        vecs = {g["backbone"] for g in groups.values()}
        vlabel = (f"{vecs.pop()} (uL)" if len(vecs) == 1 else "Vector (uL)")
        sheet, tbl = GP.build(rows, args.start_tube,
                              args.plan_title or args.file_name,
                              __import__("datetime").date.today().strftime("%Y%m%d"),
                              vlabel)
        plan_dir = os.path.join(args.outdir, "gga_plan")
        os.makedirs(plan_dir, exist_ok=True)
        import csv as _csv
        base_p = os.path.join(plan_dir, "gga_plan")
        with open(f"{base_p}.tsv", "w", newline="") as h:
            _csv.writer(h, delimiter="\t").writerows(sheet)
        with open(f"{base_p}.csv", "w", newline="") as h:
            _csv.writer(h).writerows(sheet)
        for key in ("anneal", "gga"):
            with open(f"{base_p}_{key}.csv", "w", newline="") as h:
                _csv.writer(h).writerows(tbl[key])
        print(f"GGA plan: {tbl['n_designs']} GGAs, {tbl['n_anneal']} annealings "
              f"-> {base_p}.tsv (paste into a new sheet tab)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
