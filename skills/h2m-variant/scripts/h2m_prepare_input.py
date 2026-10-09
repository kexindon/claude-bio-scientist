"""Turn a public-database export or a gene list into a model-ready H2M table.

No genome is loaded anywhere in this file. Every function here reads either a
small text file or the JSON dictionaries packaged inside the installed bioh2m
wheel, so the whole script runs in about a second.

The four-step chain is lifted from h2m-public/h2m_package_tutorial.ipynb:
  cell 21/23  h2m.cbio_reader(path)                    cBioPortal / GENIE / TCGA
  cell 25/26  h2m.vcf_reader(path) + h2m.vcf_to_maf()  gnomAD browser CSV
  cell 28     h2m.clinvar_reader(path, ids) + clinvar_to_maf(), ID -> 'index'
  cell 32     h2m.get_tx_batch(df, 'h', ver)           human canonical tx
  cell 35     h2m.query_batch(df, 'h2m')               mouse ortholog
  cell 38     h2m.get_tx_batch(df, 'm')                mouse canonical tx
Do not edit the call sequence here without editing the tutorial -- they must
not drift. The wrapper adds: an argparse CLI, a unique 'index' column, and a
dropped-rows file at every step so no variant vanishes silently.
"""

import argparse
import os
import sys

import pandas as pd

import bioh2m as h2m

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# Fast path: ortholog + canonical transcript, no genome, no network
# ---------------------------------------------------------------------------

def lookup_genes(genes, ver=37, direction="h2m"):
    """Ortholog + canonical transcript for a handful of gene names.

    Sub-second. This is the right answer to "what's the mouse gene for KIT"
    and it must not route through load_refs().

    h2m.query raises rather than returning empty for an unknown gene, and in
    m2h mode it indexes homo_dict_mouse directly so an unknown mouse gene is a
    bare KeyError. Both are caught here.
    """
    rows = []
    for g in genes:
        rec = {"query": g, "direction": direction}
        try:
            hits = h2m.query(g, direction=direction, show=False)
        except (ValueError, KeyError) as exc:
            rec["error"] = str(exc).split(".")[0]
            rows.append(rec)
            continue
        if not hits or hits[0] is None:
            rec["error"] = "no ortholog in the packaged homology dictionary"
            rows.append(rec)
            continue
        for hit in hits:
            out = dict(rec)
            out.update(hit)
            rows.append(out)
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # attach canonical transcripts for whatever resolved
    ok = df[df.get("error").isna()] if "error" in df.columns else df
    if len(ok):
        if direction == "h2m":
            tmp = pd.DataFrame({"gene_name_h": ok["query"].values,
                                "gene_name_m": ok["gene_name_m"].values})
            tmp["index"] = range(len(tmp))
            tmp, _ = h2m.get_tx_batch(tmp, "h", ver=ver)
            tmp, _ = h2m.get_tx_batch(tmp, "m")
        else:
            tmp = pd.DataFrame({"gene_name_m": ok["query"].values,
                                "gene_name_h": ok["gene_name_h"].values})
            tmp["index"] = range(len(tmp))
            tmp, _ = h2m.get_tx_batch(tmp, "h", ver=ver)
            tmp, _ = h2m.get_tx_batch(tmp, "m")
        keep = [c for c in ("gene_name_h", "gene_name_m", "tx_id_h", "tx_id_m",
                            "ref_genome_h") if c in tmp.columns]
        df = df.merge(tmp[keep].drop_duplicates(),
                      on=[c for c in ("gene_name_h", "gene_name_m") if c in df.columns
                          and c in keep],
                      how="left")
    return df


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_source(source, path, clinvar_ids=None, keep=False):
    """Dispatch to the right bioh2m reader and return a MAF-convention frame."""
    if source == "cbio":
        # cbio_reader sets tx_id_h from the HGVSc field, so these rows already
        # carry a transcript and must NOT be overwritten by get_tx_batch.
        df = h2m.cbio_reader(path, keep=keep)
    elif source == "gnomad":
        df = h2m.vcf_reader(path, keep=keep)
        df = h2m.vcf_to_maf(df)
    elif source == "clinvar":
        if not clinvar_ids:
            raise SystemExit("--source clinvar needs --clinvar-ids (comma-separated "
                             "ClinVar variation IDs); reading the whole VCF is a bad idea.")
        df = h2m.clinvar_reader(path, clinvar_ids, keep=keep)
        df = h2m.clinvar_to_maf(df)
        if "ID" in df.columns:
            df = df.rename(columns={"ID": "index"})
    elif source == "maf":
        # A plain CSV the user wrote by hand or exported from somewhere else.
        df = pd.read_csv(path)
    else:
        raise SystemExit("unknown --source " + repr(source))
    return df


def prepare(df, ver=37, has_tx=False, gene=None):
    """Add type_h, a unique index, and both canonical transcripts.

    Returns (model_ready, dropped) where `dropped` records every row lost and
    why, so a shrinking row count is never a mystery.
    """
    dropped = []

    # get_variant_type mutates its argument in place, so copy first.
    df = df.copy()

    # A gnomAD browser export is a single-gene download and carries no gene
    # column at all; the tutorial (cell 25) hardcodes df['gene_name_h'] = 'TP53'.
    if gene:
        df["gene_name_h"] = gene
    if "type_h" not in df.columns:
        df = h2m.get_variant_type(df, "ref_seq_h", "alt_seq_h", col_name="type_h")

    # Empty alleles are '' or '-' in MAF convention. get_variant_type fills NaN
    # with '' which then reads as INS/DEL, so normalize explicitly.
    for col in ("ref_seq_h", "alt_seq_h"):
        if col in df.columns:
            df[col] = df[col].fillna("-").replace("", "-")

    if "index" not in df.columns:
        df = df.reset_index(drop=True)
        df["index"] = range(len(df))
    if df["index"].duplicated().any():
        df = df.reset_index(drop=True)
        df["index"] = range(len(df))
        print("note: input 'index' had duplicates; renumbered 0..n-1")

    if "gene_name_h" not in df.columns:
        raise SystemExit(
            "Input needs a 'gene_name_h' column (HUGO symbol) to find the mouse "
            "ortholog. cbio_reader supplies it from Hugo_Symbol; a gnomAD browser "
            "export does NOT -- it is a single-gene download with no gene column, "
            "so pass --gene TP53."
        )

    # human canonical transcript, unless the reader already gave us one
    if not has_tx or "tx_id_h" not in df.columns:
        n0 = len(df)
        df, bad = h2m.get_tx_batch(df, "h", ver=ver)
        if len(bad):
            bad = bad.copy()
            bad["dropped_at"] = "get_tx_batch(h): gene not in the packaged transcript dict"
            dropped.append(bad)
        print("human transcript: " + str(len(df)) + "/" + str(n0) + " rows kept")

    # mouse ortholog -- note this FANS OUT: one row per human/mouse pair
    n0 = len(df)
    df, bad = h2m.query_batch(df, direction="h2m")
    if len(bad):
        bad = bad.copy()
        bad["dropped_at"] = "query_batch: no mouse ortholog"
        dropped.append(bad)
    if len(df) > n0:
        print("mouse ortholog: " + str(n0) + " rows became " + str(len(df))
              + " -- some genes have more than one mouse ortholog. Each original "
                "variant is still identified by its 'index'.")
    else:
        print("mouse ortholog: " + str(len(df)) + "/" + str(n0) + " rows kept")

    # mouse canonical transcript
    n0 = len(df)
    df, bad = h2m.get_tx_batch(df, "m")
    if len(bad):
        bad = bad.copy()
        bad["dropped_at"] = "get_tx_batch(m): mouse gene not in the packaged transcript dict"
        dropped.append(bad)
    print("mouse transcript: " + str(len(df)) + "/" + str(n0) + " rows kept")

    # query_batch merges with how='outer', so a gene present in the homology
    # dict but absent from the variant table arrives as an all-NaN variant row.
    before = len(df)
    df = df.dropna(subset=["start_h", "end_h", "tx_id_h", "tx_id_m"]).reset_index(drop=True)
    if len(df) < before:
        print("dropped " + str(before - len(df)) + " row(s) with no coordinates "
              "(outer-merge artifacts from query_batch)")

    for col in ("start_h", "end_h"):
        df[col] = df[col].astype(int)

    _warn_assembly(df, ver)

    df_dropped = pd.concat(dropped, ignore_index=True) if dropped else pd.DataFrame()
    return df, df_dropped


def _warn_assembly(df, ver):
    """Cheap assembly sanity check against the GENCODE .db -- no genome needed.

    Nothing in bioh2m cross-checks your coordinates against your `ver`. A
    GRCh38 coordinate run with ver=37 comes back as class 5 for every row,
    three minutes after you started. gnomAD v4 is GRCh38-only, so this is a
    real and common mistake: a TP53 variant is at chr17:7,577,xxx on GRCh37
    and chr17:7,661,xxx on GRCh38.

    Compares each sampled variant to its transcript's span and complains only
    when EVERY sampled row is more than 50 kb outside it. The threshold is
    measured, not guessed: across 15 genes checked on these two .db files the
    GRCh37->GRCh38 shift puts a coordinate 67 kb (EGFR) to 1.7 Mb (PTEN)
    outside the other assembly's span, while a variant that merely belongs to
    a non-canonical isoform of the same locus sits within a few kb of it. A
    gnomAD gene export covers every transcript at the locus, so those
    legitimately-outside rows are routine and must not trigger a warning.
    """
    try:
        import h2m_refs as _refs
    except ImportError:
        return
    try:
        db = h2m.anno_loader(_refs.anno_path("h", ver))
    except SystemExit:
        return

    sample = df.drop_duplicates(subset="tx_id_h").head(25)
    far = 0
    checked = 0
    for _, row in sample.iterrows():
        try:
            feat = db[row["tx_id_h"]]
        except Exception:
            continue
        checked += 1
        pos = int(row["start_h"])
        distance = max(feat.start - pos, pos - feat.end, 0)
        if distance > 50000:
            far += 1
    if checked and far == checked:
        print("")
        print("!! every sampled coordinate is more than 50 kb from its own transcript "
              "in GRCh" + str(ver) + ".")
        print("   These are almost certainly coordinates from the other assembly; "
              "every row would come back class 5.")
        print("   Re-run with --ver " + ("38" if ver == 37 else "37") + ". "
              "gnomAD v4 exports are GRCh38; AACR-GENIE and the lab's ClinVar "
              "file are GRCh37.")


def main():
    ap = argparse.ArgumentParser(
        description="Build a model-ready H2M input table, or look up orthologs fast.",
    )
    ap.add_argument("path", nargs="?", help="input file (omit when using --genes)")
    ap.add_argument("-o", "--out", required=True, help="output prefix (path, no extension)")
    ap.add_argument("--source", default="maf",
                    choices=("cbio", "gnomad", "clinvar", "maf"))
    ap.add_argument("--genes", help="comma-separated gene names: ortholog + canonical "
                                    "transcript only, no genome load, sub-second")
    ap.add_argument("--direction", default="h2m", choices=("h2m", "m2h"),
                    help="only affects --genes")
    ap.add_argument("--clinvar-ids", help="comma-separated ClinVar variation IDs")
    ap.add_argument("--gene", help="HUGO symbol to stamp on every row; required for a "
                                   "gnomAD export, which has no gene column")
    ap.add_argument("--ver", type=int, default=37, choices=(37, 38))
    ap.add_argument("--keep-columns", action="store_true",
                    help="keep every original column from the source file")
    args = ap.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        df = lookup_genes(genes, ver=args.ver, direction=args.direction)
        path = args.out + "_orthologs.csv"
        df.to_csv(path, index=False)
        print(df.to_string(index=False))
        print("wrote " + path)
        return

    if not args.path:
        raise SystemExit("give an input file, or use --genes for a quick ortholog lookup")

    ids = None
    if args.clinvar_ids:
        ids = [int(x) for x in args.clinvar_ids.split(",")]

    df = read_source(args.source, args.path, clinvar_ids=ids, keep=args.keep_columns)
    print("read " + str(len(df)) + " rows from " + args.path)

    ready, dropped = prepare(df, ver=args.ver, has_tx=(args.source == "cbio"),
                             gene=args.gene)

    ready_path = args.out + "_h2m_input.csv"
    ready.to_csv(ready_path, index=False)
    print("")
    print("wrote " + ready_path + "  (" + str(len(ready)) + " rows, "
          + str(ready["index"].nunique()) + " variants)")
    if len(dropped):
        drop_path = args.out + "_dropped.csv"
        dropped.to_csv(drop_path, index=False)
        print("wrote " + drop_path + "  (" + str(len(dropped)) + " rows -- read this)")


if __name__ == "__main__":
    main()
