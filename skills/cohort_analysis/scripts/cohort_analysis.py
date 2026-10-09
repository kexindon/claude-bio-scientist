#!/usr/bin/env python3
"""cohort_analysis CLI -- one command per question.

Three subcommands, split by DATA SOURCE:

  freq     public cohort (GENIE / MSK-HEME / MSK-IMPACT) -> per-variant or
           per-gene frequency table, each number carrying its denominator
  comut    patient x gene binary matrix -> Fisher co-occurrence / mutual
           exclusivity with BH FDR reported both within-focal and globally
  surv     LOCAL patient-level cohort table -> Cox PH + log-rank + KM frames,
           under a named contrast

Examples:

    python cohort_analysis.py freq --cohort heme \
        --disease-detailed 'Acute Myeloid Leukemia' \
        --genes DNMT3A,TET2,ASXL1 --classifications LOF -o OUT

    python cohort_analysis.py freq --cohort genie --disease 'Leukemia' \
        --genes DNMT3A --protein-change R882H -o OUT

    python cohort_analysis.py comut --cohort genie \
        --disease 'B-Lymphoblastic Leukemia/Lymphoma' \
        --focal TP53 --vs top50 -o OUT

    python cohort_analysis.py surv --table /path/outside/repo/final_data_matrix.RData \
        --subset de_novo --focal-group ch --vs top20 \
        --contrast four_arm --min-per-arm 10 -o OUT

The cohort file for `surv` is never copied into the repo and its path is never
hardcoded; pass it each time. Public cohorts are reached through
$LAB_DATA_ROOT -- see docs/data-layout.md.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pandas as pd  # noqa: E402

import cohorts as C  # noqa: E402
import comut as M  # noqa: E402


def _split(value):
    if not value:
        return []
    return [v.strip() for v in value.split(",") if v.strip()]


def _write(df, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    df.to_csv(path, index=False)
    print("wrote " + path + "  (%d rows)" % len(df))
    return path


def _provenance(outdir, payload):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "provenance.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    print("wrote " + path)


# --------------------------------------------------------------------------
# freq
# --------------------------------------------------------------------------

def cmd_freq(args):
    genes = _split(args.genes)
    classifications = None
    if args.classifications:
        classifications = (C.LOF_CLASSES if args.classifications.upper() == "LOF"
                           else _split(args.classifications))

    mut, samp, paths = C.load_cohort(args.cohort, args.version, genes=genes or None)
    mut, prov = C.restrict_to_disease(mut, samp, cancer_type=args.disease,
                                      cancer_type_detailed=args.disease_detailed)
    if genes:
        mut = mut[mut["Hugo_Symbol"].isin(genes)].reset_index(drop=True)

    cosmic = None
    if args.cosmic_role:
        cosmic = C.load_census(args.census)

    n_cohort = prov["n_samples"]
    if args.per == "patient" and prov["column"] is not None:
        ids = samp[samp[prov["column"]].isin(prov["value"])]["SAMPLE_ID"]
        n_cohort = int(ids.map(C.patient_id).nunique())
    elif args.per == "patient":
        n_cohort = int(samp["SAMPLE_ID"].map(C.patient_id).nunique())

    pc = C.normalize_protein_change(_split(args.protein_change)) \
        if args.protein_change else None

    variants = C.variant_frequency(mut, n_cohort, per=args.per,
                                   classifications=classifications,
                                   protein_change=pc)
    _write(variants, args.out, "variant_frequency.csv")

    by_gene = C.gene_frequency(mut, n_cohort, per=args.per, cosmic=cosmic)
    _write(by_gene, args.out, "gene_frequency.csv")

    if args.fallback:
        fb_mut, fb_samp, fb_paths = C.load_cohort(args.fallback, None,
                                                  genes=genes or None)
        hits = []
        for gene in genes:
            wanted = pc or sorted(set(
                mut.loc[mut["Hugo_Symbol"] == gene, "HGVSp_Short"].dropna()))
            hit = C.select_variants(gene, wanted, mut, fb_mut,
                                    classifications=classifications,
                                    primary_name=args.cohort,
                                    fallback_name=args.fallback)
            if len(hit):
                hits.append(hit)
        if hits:
            _write(pd.concat(hits, ignore_index=True), args.out,
                   "variants_with_fallback.csv")

    _provenance(args.out, {
        "command": "freq", "cohort": paths["name"], "version": paths["version"],
        "genome_build": paths["build"], "disease_column": prov["column"],
        "disease_value": prov["value"], "n_samples_in_cohort": prov["n_samples"],
        "denominator_used": n_cohort, "counted_per": args.per,
        "classifications": classifications, "genes": genes,
        "scope": "point mutations and indels only -- fusions live in "
                 "data_sv.txt and copy number in the CNA files, neither covered",
    })


# --------------------------------------------------------------------------
# comut
# --------------------------------------------------------------------------

def _binary_from_maf(mut, per):
    work = mut
    col = "Tumor_Sample_Barcode"
    if per == "patient":
        work = C.add_patient_col(mut)
        col = "Patient_ID"
    return M.build_binary_matrix(work, patient_col=col, gene_col="Hugo_Symbol")


def cmd_comut(args):
    focal = _split(args.focal)
    if not focal:
        raise SystemExit("--focal is required (one or more genes, comma separated)")

    if args.matrix:
        mat = pd.read_csv(args.matrix, index_col=0)
        prov = {"column": None, "value": None, "n_samples": int(len(mat))}
        paths = {"name": "user matrix", "version": None, "build": None}
    else:
        mut, samp, paths = C.load_cohort(args.cohort, args.version)
        mut, prov = C.restrict_to_disease(mut, samp, cancer_type=args.disease,
                                          cancer_type_detailed=args.disease_detailed)
        if args.snp_only:
            mut = mut[mut["Variant_Type"] == "SNP"].reset_index(drop=True)
            print("WARNING: --snp-only drops every indel, including the "
                  "clinically important ones. co-mutation-leukemia.ipynb did "
                  "this before testing; it is not the default here.")
        mat = _binary_from_maf(mut, args.per)

    if args.vs and args.vs.startswith("top"):
        n = int(args.vs[3:])
        ranked = mat.sum(axis=0).sort_values(ascending=False)
        comparators = [g for g in ranked.index.tolist() if g not in focal][:n]
        print("comparators: top %d genes by mutated %ss" % (n, args.per))
    else:
        comparators = _split(args.vs)
    if not comparators:
        raise SystemExit("--vs is required: a gene list, or topN")

    res = M.pairwise_grid(mat, focal, comparators, q_thresh=args.q_thresh)
    if len(res) == 0:
        raise SystemExit("no pairs could be tested -- none of the genes are in "
                         "the matrix")
    _write(res, args.out, "cooccurrence.csv")
    _write(M.comut_volcano_frame(res, correction=args.correction), args.out,
           "cooccurrence_volcano_frame.csv")
    _write(mat.reset_index(), args.out, "patient_gene_binary.csv")

    n_floored = int(res["or_floored"].sum())
    if n_floored:
        print("%d of %d pairs have a floored odds ratio (0 double mutants): "
              "log_or is a floor, not an effect size -- the or_floored column "
              "marks them" % (n_floored, len(res)))

    _provenance(args.out, {
        "command": "comut", "cohort": paths["name"], "version": paths["version"],
        "disease_column": prov["column"], "disease_value": prov["value"],
        "n_patients_in_matrix": int(len(mat)), "counted_per": args.per,
        "q_thresh": args.q_thresh, "correction_for_figures": args.correction,
        "fdr_families": "q_within_focal = BH within each focal gene "
                        "(CH.ipynb); q_global = BH across the whole grid "
                        "(Figure_2.R line 451). They differ.",
        "scope": "point mutations and indels only; a negative co-occurrence "
                 "result for a CNA-driven gene (CDKN2A, IKZF1) is "
                 "uninterpretable from this table alone",
    })


# --------------------------------------------------------------------------
# surv
# --------------------------------------------------------------------------

def cmd_surv(args):
    import survival as S  # imported here so freq/comut need no lifelines

    df = S.load_cohort_table(args.table, object_name=args.rdata_object,
                             subset=args.subset)
    rep = S.cohort_qc(df)
    S.print_qc(rep)

    sample_gene = S.unique_sample_gene(df)
    base_surv = S.base_survival(df)

    if args.focal_group == "ch":
        focal = [g for g in S.CH_UNIQUE if g in set(sample_gene["Gene"])]
        absent = [g for g in S.CH_UNIQUE if g not in set(sample_gene["Gene"])]
        if absent:
            print("CH genes absent from this cohort: " + ", ".join(absent))
    else:
        focal = _split(args.focal)
    if not focal:
        raise SystemExit("--focal is required (or --focal-group ch)")

    freq = (sample_gene.groupby("Gene")["Sample"].nunique()
            .sort_values(ascending=False).rename("n_patients").reset_index())
    _write(freq, args.out, "gene_frequency.csv")

    if args.vs and args.vs.startswith("top"):
        n = int(args.vs[3:])
        comparators = freq.loc[~freq["Gene"].isin(focal)].head(n)["Gene"].tolist()
    else:
        comparators = _split(args.vs)
    if not comparators:
        raise SystemExit("--vs is required: a gene list, or topN")

    covariates = _split(args.covariates) or None
    if covariates:
        missing = [c for c in covariates if c not in df.columns]
        if missing:
            raise SystemExit("covariate column(s) not in the table: " +
                             ", ".join(missing))
        extra = df[["Sample"] + covariates].drop_duplicates(subset=["Sample"])
        base_surv = base_surv.merge(extra, on="Sample", how="left")

    res, skipped, plot_data = S.pair_survival_grid(
        base_surv, sample_gene, focal, comparators,
        contrast=args.contrast, min_per_arm=args.min_per_arm,
        covariates=covariates)

    _write(res, args.out, "pair_survival_summary.csv")
    if len(skipped):
        _write(skipped, args.out, "pairs_skipped.csv")

    if args.km_pair:
        if args.km_pair not in plot_data:
            raise SystemExit("pair " + repr(args.km_pair) + " was not tested; "
                             "see pairs_skipped.csv")
        pdf = plot_data[args.km_pair]
        _write(S.logrank_matrix(pdf), args.out, "logrank_%s.csv"
               % args.km_pair.replace(" + ", "_").replace(" ", ""))
        arms = []
        for label, times, events in S.km_frame(pdf):
            arms.append(pd.DataFrame({"arm": label, "time_years": times,
                                      "event": events}))
        _write(pd.concat(arms, ignore_index=True), args.out,
               "km_frame_%s.csv" % args.km_pair.replace(" + ", "_").replace(" ", ""))

    _provenance(args.out, {
        "command": "surv", "table": os.path.basename(args.table),
        "note": "patient-level cohort data stays outside the repo; only this "
                "basename is recorded, never the contents or the full path",
        "subset": args.subset, "contrast": args.contrast,
        "contrast_meaning": S._CONTRAST_DOC[args.contrast],
        "min_per_arm": args.min_per_arm,
        "min_per_arm_note": "CH.ipynb used 5 per arm; Figure_2.R used "
                            "n_pts_both >= 10, so the two disagree on which "
                            "pairs exist at all",
        "covariates": covariates or "none -- UNIVARIATE, unadjusted",
        "qc": rep,
    })


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("freq", help="public cohort variant/gene frequency")
    p.add_argument("--cohort", default="genie", help="genie | heme | impact")
    p.add_argument("--version", help="GENIE version, default v17.0")
    p.add_argument("--disease", help="CANCER_TYPE value")
    p.add_argument("--disease-detailed", help="CANCER_TYPE_DETAILED value")
    p.add_argument("--genes", help="comma-separated Hugo_Symbol list")
    p.add_argument("--protein-change", help="comma-separated HGVSp_Short")
    p.add_argument("--classifications", help="'LOF' or a comma-separated list")
    p.add_argument("--per", default="patient", choices=["patient", "sample"])
    p.add_argument("--fallback", help="cohort to fall back to when a named "
                                      "variant is absent, e.g. genie")
    p.add_argument("--cosmic-role", action="store_true",
                   help="annotate genes with COSMIC Role in Cancer")
    p.add_argument("--census", help="explicit path to Census_all.csv")
    p.add_argument("-o", "--out", required=True)
    p.set_defaults(func=cmd_freq)

    p = sub.add_parser("comut", help="Fisher co-occurrence / mutual exclusivity")
    p.add_argument("--cohort", default="genie")
    p.add_argument("--version")
    p.add_argument("--disease")
    p.add_argument("--disease-detailed")
    p.add_argument("--matrix", help="skip the cohort; a patient x gene 0/1 csv")
    p.add_argument("--focal", help="comma-separated focal genes")
    p.add_argument("--vs", help="comma-separated comparators, or topN")
    p.add_argument("--per", default="patient", choices=["patient", "sample"])
    p.add_argument("--q-thresh", type=float, default=0.05)
    p.add_argument("--correction", default="q_within_focal",
                   choices=["q_within_focal", "q_global", "p"],
                   help="which column drives the volcano frame's y and category")
    p.add_argument("--snp-only", action="store_true",
                   help="drop every indel before testing (not recommended)")
    p.add_argument("-o", "--out", required=True)
    p.set_defaults(func=cmd_comut)

    p = sub.add_parser("surv", help="Cox PH / log-rank / KM on a local cohort")
    p.add_argument("--table", required=True,
                   help="path to the cohort table (.RData/.rds/.csv/.tsv), "
                        "OUTSIDE the repo")
    p.add_argument("--rdata-object", help="object name inside a multi-object .RData")
    p.add_argument("--subset", help="value of the Subset column to keep, e.g. de_novo")
    p.add_argument("--focal", help="comma-separated focal genes")
    p.add_argument("--focal-group", choices=["ch"],
                   help="use the 10-gene CH list from CH.ipynb cell 1")
    p.add_argument("--vs", help="comma-separated comparators, or topN")
    # kept literal rather than imported from survival.CONTRASTS so that
    # --help works without lifelines installed
    p.add_argument("--contrast", default="double_vs_driver_only",
                   choices=["double_vs_driver_only", "double_vs_rest",
                            "four_arm"])
    p.add_argument("--min-per-arm", type=int, default=5,
                   help="5 reproduces CH.ipynb; Figure_2.R used 10")
    p.add_argument("--covariates", help="comma-separated columns to adjust for, "
                                        "e.g. Age,Sex,Cytogenetics")
    p.add_argument("--km-pair", help="'DNMT3A + RAD21' -- emit its KM frame "
                                     "and pairwise log-rank table")
    p.add_argument("-o", "--out", required=True)
    p.set_defaults(func=cmd_surv)

    args = ap.parse_args(argv)
    if not getattr(args, "func", None):
        ap.print_help()
        return 1
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
