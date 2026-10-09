"""Triage an H2M result table by class, and write the hand-off CSV.

The class 0-6 vocabulary and its statement strings are lifted verbatim from
  h2m-public/h2m_package_tutorial.ipynb cell 140 (list_of_class / list_of_staterment)
and cross-checked against where each one is actually assigned in
  h2m-public/h2m-package/src/bioh2m/module.py
    class 5  L459-461      class 4  L472-474, L478-480 (GateFlankH)
    class 6  L510-526, L915, L1944
    class 0  L514-534, L1196, L1939
    class 2  L861, L1689, L1786, L1907
    class 1  L1581
    class 3  L1637
Do not reword these without re-reading those lines.

Nothing here recomputes any H2M decision. It groups, counts and filters the
table bioh2m produced.
"""

import argparse
import os

import pandas as pd


# Verbatim from tutorial cell 140, plus what each one means for the bench.
CLASS_INFO = {
    0: ("This mutation can be originally modeled.",
        "ORDERABLE -- the identical nucleotide change works in mouse. "
        "*_m and *_m_ori are the same variant."),
    1: ("This mutation can be alternatively modeled.",
        "ORDERABLE -- a different mouse codon gives the same amino acid change. "
        "Order *_m (NOT *_m_ori, which is the identical-nucleotide version "
        "and makes the wrong protein)."),
    2: ("This mutation can be modeled, but the effect may not be consistent.",
        "JUDGEMENT CALL -- frameshift, splice-site and non-coding variants. "
        "H2M can place a change at the orthologous position but does not vouch "
        "for the protein consequence, and there is no amino-acid-matched "
        "alternative. Read HGVSp_m yourself before ordering."),
    3: ("This mutation cannot be originally modeled and no alternative is found.",
        "NOT MODELABLE -- drop it, or model a different variant in the same gene."),
    4: ("Mutated sequences are not identical.",
        "ALIGNMENT PROBLEM -- the human and mouse proteins do not align "
        "contiguously across this variant, so the orthologous position is "
        "ambiguous. This is about homology, NOT about your input."),
    5: ("Coordinate error. This mutation is not in the query gene.",
        "YOUR INPUT IS WRONG -- the coordinate is not inside the transcript. "
        "Almost always a GRCh37/GRCh38 mismatch between your coordinates, "
        "--ver, and the GENCODE .db."),
    6: ("This mutation cannot be originally modeled.",
        "NOT MODELABLE AS GIVEN -- only appears when the alternative search is "
        "switched off (--no-alternative) or under --param BE with "
        "search_alternative=False. With defaults you get class 3 instead."),
}

ORDERABLE = (0, 1)

# What library_design / pegg needs, and nothing else.
HANDOFF_COLUMNS = [
    "index", "class", "match",
    "gene_name_h", "tx_id_h", "HGVSp_h", "classification_h",
    "gene_name_m", "tx_id_m", "chr_m", "strand_m",
    "start_m", "end_m", "ref_seq_m", "alt_seq_m", "type_m",
    "HGVSc_m", "HGVSp_m", "classification_m", "exon_m",
    "flank_size_left", "flank_size_right",
]


def _flank_aa(value):
    """flank_size_left/right come back as strings like '4aa', not integers.

    Returns an int, or None if it cannot be parsed.
    """
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return None
    digits = "".join(ch for ch in s if ch.isdigit())
    return int(digits) if digits else None


def per_variant(df):
    """One row per input variant, choosing the best alternative available.

    h2m.model returns a LIST of dicts -- one per alternative -- and model_batch
    flattens all of them into the result frame keyed on 'index'. Counting rows
    therefore overcounts variants. Alternatives for one variant can also be
    literal duplicates of each other.

    Picks the lowest class per index, because 0 beats 1 beats 2 beats 3.
    """
    if "index" not in df.columns:
        raise SystemExit("result table has no 'index' column -- is this an H2M result CSV?")
    d = df.copy()
    d["_cls"] = pd.to_numeric(d["class"], errors="coerce")
    d = d.sort_values(["index", "_cls"], kind="mergesort")
    out = d.drop_duplicates(subset="index", keep="first").drop(columns=["_cls"])
    return out.reset_index(drop=True)


def alternatives(df, idx):
    """Every distinct alternative for one variant, deduplicated."""
    sub = df[df["index"] == idx]
    cols = [c for c in ("class", "HGVSp_m_ori", "HGVSp_m", "HGVSc_m",
                        "start_m", "ref_seq_m", "alt_seq_m") if c in sub.columns]
    return sub[cols].drop_duplicates().reset_index(drop=True)


def report(df, min_flank=2):
    """Print the triage. Returns the per-variant frame."""
    pv = per_variant(df)
    n = len(pv)
    print("=" * 72)
    print(str(len(df)) + " result rows  ->  " + str(n) + " distinct variants")
    print("(rows exceed variants because each alternative modeling is its own row)")
    print("=" * 72)
    print("")

    counts = pv["class"].value_counts().sort_index()
    for cls in sorted(CLASS_INFO):
        k = int(counts.get(cls, 0))
        if k == 0:
            continue
        statement, meaning = CLASS_INFO[cls]
        print("class " + str(cls) + "  n=" + str(k) + "   " + statement)
        for line in meaning.split(" -- "):
            print("          " + line)
        print("")

    n_ord = int(pv["class"].isin(ORDERABLE).sum())
    print("-" * 72)
    print("MODELABLE (class 0 or 1): " + str(n_ord) + "/" + str(n)
          + "  (" + str(round(100.0 * n_ord / n, 1)) + "%)" if n else "")
    print("-" * 72)

    # match == False is NOT a class. Verified in module.py: match_h is set at
    # L1233/1235/1251 and never touches `error`, so a variant whose reference
    # allele disagrees with the genome still comes back class 0 with
    # status True. H2M quietly modeled the genome's allele, not yours.
    if "match" in pv.columns:
        bad = pv[pv["match"] == False]  # noqa: E712
        if len(bad):
            print("")
            print("!! " + str(len(bad)) + " variant(s) have match == False.")
            print("   Your reference allele disagrees with the genome at that position.")
            print("   H2M does NOT give these a class of their own -- they are reported")
            print("   as class 0/1 having modeled the genome's allele instead of yours.")
            print("   Usually a strand mistake (coordinates and alleles must both be on")
            print("   the chromosome's POSITIVE strand) or a wrong assembly.")
            for _, r in bad.head(10).iterrows():
                print("     index " + str(r["index"]) + "  " + str(r.get("gene_name_h"))
                      + " " + str(r.get("HGVSp_h")))

    # Low local homology: flank_size_* are strings like '4aa'.
    if "flank_size_left" in pv.columns:
        lo = []
        for _, r in pv[pv["class"].isin(ORDERABLE)].iterrows():
            left = _flank_aa(r.get("flank_size_left"))
            right = _flank_aa(r.get("flank_size_right"))
            if left is not None and right is not None and min(left, right) < min_flank:
                lo.append((r["index"], r.get("gene_name_h"), r.get("HGVSp_h"), left, right))
        if lo:
            print("")
            print(str(len(lo)) + " modelable variant(s) sit in poorly conserved local context")
            print("(fewer than " + str(min_flank) + " identical flanking amino acids on one side):")
            for idx, gene, hgvsp, left, right in lo[:10]:
                print("     index " + str(idx) + "  " + str(gene) + " " + str(hgvsp)
                      + "   flank " + str(left) + "aa / " + str(right) + "aa")
            print("Trust these less; the orthologous position is real but the")
            print("surrounding sequence is not, which matters for guide design.")

    return pv


def main():
    ap = argparse.ArgumentParser(
        description="Triage an H2M result CSV by class 0-6 and write the hand-off table.",
    )
    ap.add_argument("result_csv")
    ap.add_argument("-o", "--out", help="write modeled_variants.csv next to this prefix")
    ap.add_argument("--min-flank", type=int, default=2,
                    help="flag class 0/1 variants with fewer identical flanking aa than this")
    ap.add_argument("--keep-class", default="0,1",
                    help="classes to put in the hand-off CSV (default 0,1)")
    ap.add_argument("--keep-mismatched", action="store_true",
                    help="keep rows with match == False (default: excluded)")
    ap.add_argument("--show", type=str,
                    help="print every alternative for this 'index' value and stop")
    args = ap.parse_args()

    df = pd.read_csv(args.result_csv)

    if args.show is not None:
        # 'index' may be int64 (renumbered) or str (e.g. a ClinVar ID), so try
        # the declared dtype before falling back to the literal string.
        idx = args.show
        if df["index"].dtype.kind in "iu":
            try:
                idx = int(idx)
            except ValueError:
                raise SystemExit("--show expects an integer for this table's 'index'")
        sub = alternatives(df, idx)
        if sub.empty:
            raise SystemExit("no rows with index == " + repr(idx))
        print(sub.to_string(index=False))
        return

    pv = report(df, min_flank=args.min_flank)

    if args.out:
        keep = set(int(x) for x in args.keep_class.split(","))
        hand = pv[pv["class"].isin(keep)]
        if not args.keep_mismatched and "match" in hand.columns:
            n0 = len(hand)
            hand = hand[hand["match"] != False]  # noqa: E712
            if len(hand) < n0:
                print("")
                print("excluded " + str(n0 - len(hand)) + " match == False variant(s) from the "
                      "hand-off (--keep-mismatched to override)")
        cols = [c for c in HANDOFF_COLUMNS if c in hand.columns]
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        path = args.out + "_modeled_variants.csv"
        hand[cols].to_csv(path, index=False)
        print("")
        print("wrote " + path + "  (" + str(len(hand)) + " variants, classes "
              + args.keep_class + ")")
        print("Hand this path to the library_design skill. It runs in the pegg "
              "environment, so the hand-off is this file on disk -- not an import.")


if __name__ == "__main__":
    main()
