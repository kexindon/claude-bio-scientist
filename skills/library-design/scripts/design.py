"""Design pegRNAs or base-editing sgRNAs through pegg, with the validated
parameter sets, the standard post-filters, the checkpoint in the right place,
and a per-variant dropout ledger.

pegg is the ENGINE and is called through, never vendored: prime.run(),
base.run_base(), prime.prime_oligo_generator(), base.base_oligo_generator().

Parameter profiles and the filter/checkpoint ORDER are lifted from:
    early-leukemia-evolution/ch_aml_library_design.ipynb cell 16 (prime.run
        call and its justifications), cell 18 (oligo-then-checkpoint), cell 30
        (the matched silent-control run)
    B-ALL-in-vivo-base-editing/library-generation/MBESv2_library_generation.ipynb
        cell 4 (the validated run_base call)
Do not edit the profiles here without editing there -- these must not drift.

Added as wrappers, not edits: argparse, LAB_DATA_ROOT resolution (via
prep_input), the dropout ledger, and the OptiPrime group validation.
"""

import argparse
import os

import pandas as pd

import prep_input

# ---------------------------------------------------------------------------
# Validated parameter profiles.
#
# PRIME -- cell 16. PBS/RTT matched to OptiPrime's own design rule
# (DESIGN_PE.py): PBS 5-16, and RTT enumerated absolutely to 30 with
# min_RHA_size=8 doing the homology-arm work (OptiPrime's hom_len lower bound;
# pegg's old min_RHA_size=1 allowed arms far shorter than OptiPrime would
# score). RTT is capped at 30 because OptiPrime's hom_len stops at 29, pegg
# warns above 30, and a 55-nt sensor cannot hold a longer RTT -- sensor_wt comes
# back None and other_filtration() raises on it.
PRIME_PROFILE = dict(
    PAM="NGG",
    RTT_lengths=list(range(10, 31)),
    PBS_lengths=list(range(5, 17)),
    min_RHA_size=8,
    RE_sites=["CGTCTC"],          # Esp3I / BsmBI
    sensor=True,
    sensor_length=55,
    before_proto_context=5,
    pegRNAs_per_mut="All",
    seed=0,                        # bystanders are chosen at RANDOM
)

# BASE -- MBESv2 cell 4. PAM='NG' deliberately (MBESv2 is an NG-PAM editor);
# note prime warns on a non-NGG PAM and base does NOT warn at all, so nothing
# tells you the efficiency predictions are unvalidated there.
BASE_PROFILE = dict(
    PAM="NG",
    filtration="ABE+CBE",
    ideal_edit_window=[4, 8],
    proto_size=19,
    before_proto_context=10,
    sensor_length=42,
    RE_sites=["CGTCTC", "GAATTC"],   # Esp3I, EcoRI
)

SENSOR_TOO_SMALL = ("sensor too small; increase sensor size in parameters; "
                    "or decrease before_proto_context")


def validate_optiprime_group(group):
    """An unmatched optiprime_group is NOT an error inside OptiPrime.

    rx_model.load_rates() leaves group_factor at 0 and logs only at INFO, so a
    typo returns a plausible-looking score from an undefined MMR context. Cell
    type is how MMR status enters the model, so that is a silent correctness
    failure. Validate against pegg's own list and hard-fail.
    """
    from pegg import optiprime as _op
    if group not in _op.VALID_GROUPS:
        raise SystemExit(
            "optiprime_group %r is not one of OptiPrime's groups. A name that "
            "is not in this list is NOT rejected by OptiPrime -- it silently "
            "scores with group_factor=0. Valid: %s"
            % (group, ", ".join(_op.VALID_GROUPS)))
    return group


def post_filter(df, mode):
    """The standard library filters, in the notebook's order.

    pegg sets contains_polyT_terminator and contains_RE_site but does NOT drop
    those rows -- filtering is the caller's job. Both filters run AFTER run()
    has already scored and capped, so a tight pegRNAs_per_mut can leave a
    variant with zero rows here even though run() returned 50.

    The protospacer poly-T re-check is NOT in the notebook: it is here because
    prime.other_filtration never scans the protospacer in prime mode. The
    source reads `if terminator_sequence in (pbs_rtt or proto)` (prime.py:1456
    and again at 1555) -- Python evaluates `pbs_rtt or proto` to pbs_rtt, since
    a non-empty string is truthy, so only the extension is ever scanned. That
    is a real bug in pegg, not a convention. Base mode is unaffected (it only
    has a protospacer).
    """
    n0 = len(df)
    stages = []

    if "contains_polyT_terminator" in df.columns:
        df = df[df["contains_polyT_terminator"] == False]  # noqa: E712
        stages.append(("contains_polyT_terminator", n0 - len(df)))

    if mode == "prime":
        n = len(df)
        polyt = "T" * 4
        df = df[~df["Protospacer"].astype(str).str.upper().str.contains(polyt)]
        stages.append(("polyT in Protospacer (pegg misses this)", n - len(df)))

    if "sensor_error" in df.columns:
        n = len(df)
        df = df[df["sensor_error"] != SENSOR_TOO_SMALL]
        stages.append(("sensor too small", n - len(df)))

    if "contains_RE_site" in df.columns:
        n = len(df)
        df = df[df["contains_RE_site"] == False]  # noqa: E712
        stages.append(("contains_RE_site", n - len(df)))

    for name, dropped in stages:
        if dropped:
            print("  filtered %7d rows: %s" % (dropped, name))
    print("  %d -> %d rows after filtering" % (n0, len(df)))
    return df


def recast_bystander_flag(df):
    """has_silent_bystander can come back as object dtype.

    The not-designable branch pads every column with None, so the column goes
    object; run() casts it back, but a round-trip through CSV does not.
    pd.read_csv gives strings or bools depending on content, and
    ~df['has_silent_bystander'] then does INTEGER negation instead of logical
    not -- quietly inverting nothing and selecting everything. Re-cast after
    loading any checkpoint.
    """
    if "has_silent_bystander" not in df.columns:
        return df
    col = df["has_silent_bystander"]
    if col.dtype == bool:
        return df
    df = df.copy()
    df["has_silent_bystander"] = col.map(
        {True: True, False: False, "True": True, "False": False,
         "true": True, "false": False, 1: True, 0: False}).astype(bool)
    return df


def dropout_ledger(design_input, designed, filtered, mode, gene_column="Hugo_Symbol"):
    """One row per INPUT variant, naming the stage that killed it.

    "Why did my variant drop out of the library" has at least eight independent
    silent-drop stages, so the answer has to be a table. mutation_idx is the
    input row index, which is why prepare() reset_index()es.
    """
    rows = []
    designed_idx = (set(designed["mutation_idx"]) if "mutation_idx" in designed
                    else set())
    kept_idx = (set(filtered["mutation_idx"]) if "mutation_idx" in filtered
                else set())
    n_by_idx = filtered.groupby("mutation_idx").size() if len(filtered) else None

    for i, row in design_input.iterrows():
        vt = row.get("Variant_Type")
        if mode == "base" and vt != "SNP":
            reason = "not a SNP (run_base drops every non-SNP)"
        elif i not in designed_idx:
            reason = ("no pegRNA designed: no eligible PAM within reach, or "
                      "min_RHA_size, or the variant is outside its own CDS")
        elif i not in kept_idx:
            reason = "all designs removed by the post-filters"
        else:
            reason = ""
        rows.append(dict(
            mutation_idx=i,
            gene=row.get(gene_column),
            Variant_Type=vt,
            Chromosome=row.get("Chromosome"),
            Start_Position=row.get("Start_Position"),
            cds_valid=row.get("cds_valid"),
            n_designed=int((designed["mutation_idx"] == i).sum())
            if "mutation_idx" in designed else 0,
            n_kept=int(n_by_idx.get(i, 0)) if n_by_idx is not None else 0,
            dropped_at=reason,
        ))
    led = pd.DataFrame(rows)
    lost = led[led["dropped_at"] != ""]
    if len(lost):
        print("\n%d of %d input variants have NO guide in the final table:"
              % (len(lost), len(led)))
        print(lost[["gene", "Variant_Type", "dropped_at"]].to_string(index=False))
    else:
        print("\nevery one of the %d input variants kept at least one guide"
              % len(led))
    return led


def run_prime(design_input, chrom_dict, args):
    from pegg import prime

    kwargs = dict(PRIME_PROFILE)
    kwargs["seed"] = args.seed
    if args.pegrnas_per_mut != "All":
        kwargs["pegRNAs_per_mut"] = int(args.pegrnas_per_mut)
    if args.pam:
        kwargs["PAM"] = args.pam

    if args.silent_bystander:
        # No cap on either: MEASURED on this library, lifting
        # max_bystander_muts 2 -> unlimited adds only ~4% (3784 -> 3928
        # pegRNAs/variant). The binding constraint is bystander_window_nt.
        kwargs.update(silent_bystander=True,
                      silent_per_mut=10 ** 6,
                      max_bystander_muts=10 ** 6,
                      bystander_window_nt=args.bystander_window_nt)

    if args.optiprime:
        kwargs.update(
            optiprime=True,
            optiprime_group=validate_optiprime_group(args.optiprime_group),
            optiprime_jobs=args.optiprime_jobs,
            # No floor and no prefilter, on purpose. A cutoff drops whole
            # MUTATIONS (0.20 lost 17 of 40 here, Npm1/Idh1/Idh2/Tet2/Trp53
            # among them); the default prefilter=4 reinstates PEGG2_Score as a
            # gatekeeper and only 10 of 22 mutations kept an identical top-50.
            # Keep the score in the table and threshold afterwards, with the
            # coverage cost visible.
            optiprime_cutoff=None,
            optiprime_prefilter=None,
            rankby="OptiPrime_Score",
        )

    print("prime.run(%s)" % ", ".join(
        "%s=%r" % (k, v) for k, v in sorted(kwargs.items())))
    designed = prime.run(design_input, "cBioPortal", chrom_dict=chrom_dict,
                         **kwargs)
    df = post_filter(designed, "prime")
    # prime_oligo_generator BEFORE the checkpoint: it adds pegRNA_oligo, and a
    # checkpoint taken before it lacks that column, so resuming from it dies
    # with KeyError: ['pegRNA_oligo'] not in index.
    df = prime.prime_oligo_generator(df, epeg=args.epeg,
                                     epeg_motif=args.epeg_motif)
    return designed, df


def run_base(design_input, chrom_dict, args):
    from pegg import base

    kwargs = dict(BASE_PROFILE)
    if args.pam:
        kwargs["PAM"] = args.pam
    if args.filtration:
        kwargs["filtration"] = (args.filtration.split(",")
                                if "," in args.filtration else args.filtration)

    print("base.run_base(%s)" % ", ".join(
        "%s=%r" % (k, v) for k, v in sorted(kwargs.items())))
    designed = base.run_base(design_input, "cBioPortal", chrom_dict=chrom_dict,
                             **kwargs)
    df = post_filter(designed, "base")
    df = base.base_oligo_generator(df)
    return designed, df


def main():
    ap = argparse.ArgumentParser(
        description="Design pegRNAs or base-editing sgRNAs with pegg. Ends at "
                    "the design table plus the pooled-synthesis oligo column; "
                    "for an IDT sheet of Golden Gate oligos use pegrna_oligo.")
    ap.add_argument("variants", help="prepared CSV from prep_input.py")
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("--name", default="library", help="output file stem")
    ap.add_argument("--mode", choices=["prime", "base"], required=True)
    ap.add_argument("--build", required=True,
                    choices=sorted(prep_input.GENOMES),
                    help="genome build; must match the variant coordinates")
    ap.add_argument("--pam", default=None, help="override the profile's PAM")
    ap.add_argument("--seed", type=int, default=0,
                    help="bystander designs are chosen at random; pegg prints "
                         "a note when seed is None, easy to miss in a long log")
    ap.add_argument("--pegrnas-per-mut", default="All")
    ap.add_argument("--silent-bystander", action="store_true")
    ap.add_argument("--bystander-window-nt", type=int, default=5,
                    help="LEAVE AT 5 unless you have decided to pay for it: "
                         "measured pegRNAs per variant are 5->3,928 "
                         "10->19,440 15->87,984 30->1,334,448")
    ap.add_argument("--tx-column", default=None,
                    help="transcript column for the reading frame (e.g. "
                         "tx_id_m); required with --silent-bystander unless "
                         "the prepared table already carries start_end_cds")
    ap.add_argument("--gene-column", default="Hugo_Symbol")
    ap.add_argument("--filtration", default=None,
                    help="base mode: 'ABE+CBE', or a comma list for another "
                         "editor, e.g. 'C>G,T>A'")
    ap.add_argument("--epeg", action="store_true", default=True)
    ap.add_argument("--no-epeg", dest="epeg", action="store_false")
    ap.add_argument("--epeg-motif", default="tevopreQ1")
    ap.add_argument("--optiprime", action="store_true",
                    help="add OptiPrime_Score. SLOW: ~0.2155 s per pegRNA per "
                         "worker, which was 6.5 h for a 665k-row library")
    ap.add_argument("--optiprime-group", default="Liu_HeLa",
                    help="keep Liu_HeLa even in another cell type: it is "
                         "MMR-proficient, so designs transfer to HEK293T but "
                         "not the reverse")
    ap.add_argument("--optiprime-jobs", type=int, default=None,
                    help="cores - 1; each worker pays ~20 s of JAX JIT and "
                         "score() will not shard below 500 rows")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    design_input = pd.read_csv(args.variants)

    chrom_dict, _ = prep_input.load_genome(args.build)
    db = None
    needs_frame = (args.mode == "prime" and args.silent_bystander
                   and "start_end_cds" not in design_input.columns)
    if needs_frame:
        db = prep_input.load_annotation(args.build)
    design_input, _gene_cds, _key = prep_input.prepare(
        design_input, chrom_dict, db=db, tx_column=args.tx_column,
        gene_column=args.gene_column)

    if args.mode == "prime":
        designed, df = run_prime(design_input, chrom_dict, args)
    else:
        designed, df = run_base(design_input, chrom_dict, args)

    # Checkpoint the MOMENT the result exists. Scoring is the only expensive
    # step and until this line its result lives in memory alone; a NameError in
    # a later cell once destroyed a completed 6.5 h run. Gzipped, because this
    # data compresses to ~2-8% of raw.
    stem = os.path.join(args.out_dir, args.name)
    ckpt = stem + "_designs.csv.gz"
    df.to_csv(ckpt, index=False, compression="gzip")
    print("checkpoint: %s (%d rows, %d cols, %.1f MB)"
          % (os.path.basename(ckpt), df.shape[0], df.shape[1],
             os.path.getsize(ckpt) / 1e6))

    oligo_col = "pegRNA_oligo" if args.mode == "prime" else "gRNA_oligo"
    id_col = "pegRNA_id" if args.mode == "prime" else "gRNA_id"
    if oligo_col in df.columns:
        opath = stem + "_oligos.csv"
        df[[id_col, oligo_col]].to_csv(opath, index=False)
        print("wrote %s (%d pooled-synthesis oligos)"
              % (os.path.basename(opath), len(df)))

    led = dropout_ledger(design_input, designed, df, args.mode,
                         gene_column=args.gene_column)
    lpath = stem + "_dropout_ledger.csv"
    led.to_csv(lpath, index=False)
    print("wrote %s" % os.path.basename(lpath))

    if "has_silent_bystander" in df.columns:
        df = recast_bystander_flag(df)
        print("plain %d, bystander %d"
              % (int((~df["has_silent_bystander"]).sum()),
                 int(df["has_silent_bystander"].sum())))


if __name__ == "__main__":
    main()
