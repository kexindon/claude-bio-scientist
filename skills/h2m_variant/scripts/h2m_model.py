"""Batch H2M modeling: a variant table in, a result table + a leftover table out.

The modeling call and the two-file output convention are lifted verbatim from
  h2m-public/methods/run_h2m/h2m_db.py  (the lab's production cluster driver)
Specifically its line 25 and the result[0]/result[1] write-out, including
`memory_size=2000` -- the deliberate choice to fail fast on very long proteins
rather than risk the node. Do not edit here without editing there; these two
must not drift.

Everything else in this file is a wrapper: a $LAB_DATA_ROOT resolver (h2m_refs),
an argparse CLI, and column validation that runs BEFORE the genomes load so a
missing column costs one second instead of one minute.

bioh2m does the science. This file does not reimplement any of it.
"""

import argparse
import os
import sys

import pandas as pd

import bioh2m as h2m

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h2m_refs  # noqa: E402


# model_batch's own requirement, copied from module.py L2083-L2092 so we can
# check it up front instead of letting it raise "Error 2" after the load.
_REQUIRED = {
    "h": ["start_h", "end_h", "type_h", "ref_seq_h", "alt_seq_h", "tx_id_h", "tx_id_m", "index"],
    "m": ["start_m", "end_m", "type_m", "ref_seq_m", "alt_seq_m", "tx_id_h", "tx_id_m", "index"],
    "": ["start", "end", "type", "ref_seq", "alt_seq", "tx_id_h", "tx_id_m", "index"],
}

_VALID_TYPES = ("SNP", "DNP", "TNP", "ONP", "INS", "DEL")


def validate(df, coor="nc"):
    """Fail loudly on an input table model_batch would reject or mis-read.

    model_batch picks its expected column names by sniffing for 'ref_seq_h'
    then 'ref_seq_m' then the unsuffixed form, so we do the same.
    """
    cols = list(df.columns)
    if "ref_seq_h" in cols:
        flavor = "h"
    elif "ref_seq_m" in cols:
        flavor = "m"
    else:
        flavor = ""
    need = _REQUIRED[flavor]

    missing = [c for c in need if c not in cols]
    if missing:
        raise SystemExit(
            "Input is missing required column(s): " + ", ".join(missing) + ".\n"
            "model_batch needs exactly " + ", ".join(need) + ".\n"
            "If your identifier is in the pandas index rather than a column, "
            "run df.reset_index() first -- 'index' must be a real column."
        )

    if df["index"].duplicated().any():
        dups = df["index"][df["index"].duplicated()].unique()[:5]
        raise SystemExit(
            "Column 'index' has duplicate values (" + ", ".join(str(x) for x in dups) + " ...). "
            "It is the join key for alternatives and for --bind; duplicates turn the "
            "bind merge into a cartesian product. Make it unique."
        )

    suffix = "_" + flavor if flavor else ""
    tcol = "type" + suffix
    bad = sorted(set(str(x) for x in df[tcol]) - set(_VALID_TYPES))
    if bad:
        raise SystemExit(
            "Column " + tcol + " has value(s) H2M does not accept: " + ", ".join(bad) + ".\n"
            "Allowed: " + "/".join(_VALID_TYPES) + ". Use h2m.get_variant_type() to infer "
            "them from allele lengths (h2m_prepare_input.py does this for you)."
        )

    scol, ecol = "start" + suffix, "end" + suffix
    if coor == "nc":
        flipped = df[df[scol] > df[ecol]]
        if len(flipped):
            raise SystemExit(
                str(len(flipped)) + " row(s) have " + scol + " > " + ecol + ". Coordinates are "
                "1-based positions on the POSITIVE strand of the chromosome with start <= end, "
                "regardless of which strand the transcript is on."
            )
    return flavor


def run(df, ver=37, param="default", direction="h2m", coor="nc",
        max_alternative=5, memory_size=2000, search_alternative=True,
        bind=False, show_sequence=False):
    """Load the references once, then call model_batch. Returns (result, failed)."""
    validate(df, coor=coor)
    records_h, index_list_h, records_m, index_list_m, db_h, db_m = h2m_refs.load_refs(ver)

    # ---- lifted from methods/run_h2m/h2m_db.py line 25 --------------------
    # result = h2m.model_batch(df_run, records_h, index_list_h, records_m,
    #                          index_list_m, db_h, db_m, 37, memory_size=2000)
    result = h2m.model_batch(
        df, records_h, index_list_h, records_m, index_list_m, db_h, db_m, ver,
        param=param, direction=direction, coor=coor,
        search_alternative=search_alternative, max_alternative=max_alternative,
        memory_size=memory_size, show_sequence=show_sequence, bind=bind,
    )
    # -----------------------------------------------------------------------
    return result[0], result[1]


def main():
    ap = argparse.ArgumentParser(
        description="Model a table of variants across species with bioh2m.",
        epilog="Writes <out>_result.csv and <out>_leftover.csv. Always read the "
               "leftover file: model_batch swallows per-variant exceptions silently.",
    )
    ap.add_argument("input_csv")
    ap.add_argument("-o", "--out", required=True, help="output prefix (path, no extension)")
    ap.add_argument("--ver", type=int, default=37, choices=(37, 38),
                    help="human assembly; must match your coordinates AND the .db")
    ap.add_argument("--direction", default="h2m", choices=("h2m", "m2h", "h2h", "m2m"))
    ap.add_argument("--param", default="default", choices=("default", "BE"),
                    help="BE restricts output to base-editable changes")
    ap.add_argument("--coor", default="nc", choices=("nc", "aa"),
                    help="aa: start/end are residue numbers, ref/alt are 1-letter codes")
    ap.add_argument("--max-alternative", type=int, default=5)
    ap.add_argument("--memory-size", type=int, default=2000,
                    help="max aligned protein length; 2000 is the lab's production value")
    ap.add_argument("--no-alternative", action="store_true",
                    help="original modeling only (no codon-level alternatives)")
    ap.add_argument("--bind", action="store_true",
                    help="merge the input's non-coordinate columns back on 'index'")
    ap.add_argument("--show-sequence", action="store_true",
                    help="add the 14 sequence/index columns visualization() needs")
    ap.add_argument("--dry-run", action="store_true",
                    help="validate columns and reference paths, then stop")
    args = ap.parse_args()

    df = pd.read_csv(args.input_csv)
    flavor = validate(df, coor=args.coor)
    print("input OK: " + str(len(df)) + " rows, "
          + str(df["index"].nunique()) + " unique index values, "
          + (flavor + "-suffixed" if flavor else "unsuffixed") + " columns")

    if args.direction in ("h2h", "m2m"):
        print("note: paralog mode -- output columns are suffixed _1/_2, not _h/_m")

    if args.dry_run:
        for label, path in (
            ("human genome", h2m_refs.genome_path("h", args.ver)),
            ("mouse genome", h2m_refs.genome_path("m")),
            ("human anno", h2m_refs.anno_path("h", args.ver)),
            ("mouse anno", h2m_refs.anno_path("m")),
        ):
            print(label + ": " + path)
        h2m_refs.preflight(args.ver)
        print("dry run OK")
        return

    result, failed = run(
        df, ver=args.ver, param=args.param, direction=args.direction,
        coor=args.coor, max_alternative=args.max_alternative,
        memory_size=args.memory_size,
        search_alternative=not args.no_alternative,
        bind=args.bind, show_sequence=args.show_sequence,
    )

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    res_path = args.out + "_result.csv"
    left_path = args.out + "_leftover.csv"
    result.to_csv(res_path, index=False)
    failed.to_csv(left_path, index=False)

    n_in = df["index"].nunique()
    n_out = result["index"].nunique() if len(result) else 0
    print("")
    print("wrote " + res_path + "  (" + str(len(result)) + " rows for "
          + str(n_out) + " variants -- rows > variants because each "
          "alternative is its own row)")
    print("wrote " + left_path + "  (" + str(len(failed)) + " rows)")
    if n_out < n_in:
        print("WARNING: " + str(n_in - n_out) + " of " + str(n_in) + " input variants produced "
              "no result. Re-run those rows one at a time through h2m.model() to see "
              "the real exception -- model_batch catches it and tells you nothing.")
    if "class" in result.columns and len(result):
        counts = result.drop_duplicates(subset="index")["class"].value_counts().sort_index()
        print("class distribution (one row per variant): "
              + ", ".join(str(k) + ":" + str(v) for k, v in counts.items()))
    if "match" in result.columns and len(result):
        n_bad = int((result["match"] == False).sum())  # noqa: E712
        if n_bad:
            print("WARNING: " + str(n_bad) + " row(s) have match == False -- the reference allele "
                  "you supplied disagrees with the genome. H2M modeled the GENOME's allele "
                  "instead and still reported class 0/1. Fix the input; do not order these.")


if __name__ == "__main__":
    main()
