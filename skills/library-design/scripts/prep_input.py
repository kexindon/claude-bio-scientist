"""Turn a variant table into something prime.run()/run_base() will accept.

Lifted from early-leukemia-evolution/ch_aml_library_design.ipynb:
    cell 6  -> load_genome / load_annotation (paths, via LAB_DATA_ROOT)
    cell 11 -> chrom_dict_from_h2m   (the re.search('chromosome (\\w+)') walk)
    cell 13 -> maf_from_h2m, normalize_chromosome, annotate_frame,
               assert_chroms_present
Do not edit here without editing there -- these two must not drift.

What is NOT lifted (added as wrappers): the LAB_DATA_ROOT resolver, argparse,
and the named errors. The notebook hardcodes
/Users/kexindong/Documents/GitHub/Database as a fallback with a "remember to
replace the paths with yours" comment; that fallback is deliberately gone here
so the same code runs for anyone.
"""

import argparse
import os
import re

import pandas as pd

# Builds this skill knows how to resolve under $LAB_DATA_ROOT. The genome and
# the annotation must be the SAME build as the variant coordinates: GRCh37 goes
# with GENCODE v19, GRCh38 with v44, GRCm39 with vm33.
GENOMES = {
    "GRCh37": "RefGenome/ncbi-2023-09-12/GCF_000001405.25_GRCh37.p13_genomic.fna.gz",
    "GRCm39": "RefGenome/mouse-2023-09-13/GCF_000001635.27_GRCm39_genomic.fna.gz",
}
ANNOTATIONS = {
    "GRCh37": "Genecode/gencode_v19_GRCh37.db",
    "GRCh38": "Genecode/gencode_v44_GRCh38.db",
    "GRCm39": "Genecode/gencode_vm33_GRCm39.db",
}

# H2M column names -> the MAF names pegg's 'cBioPortal' reader wants.
# '_m' is the mouse orthologue H2M modelled onto; '_h' is the human original.
H2M_RENAME = {
    "gene_name": "Hugo_Symbol",
    "chr": "Chromosome",
    "start": "Start_Position",
    "end": "End_Position",
    "ref_seq": "Reference_Allele",
    "alt_seq": "Tumor_Seq_Allele2",
    "type": "Variant_Type",
}

MAF_COLUMNS = ("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
               "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2")


def data_root():
    """$LAB_DATA_ROOT, or a named error. Never a hardcoded fallback."""
    root = os.environ.get("LAB_DATA_ROOT")
    if not root:
        raise SystemExit("LAB_DATA_ROOT is not set -- see docs/data-layout.md")
    return root


def _resolve(table, build, what):
    if build not in table:
        raise SystemExit("no %s known for build %r; have: %s"
                         % (what, build, ", ".join(sorted(table))))
    path = os.path.join(data_root(), table[build])
    if not os.path.exists(path):
        raise SystemExit("%s for %s not found at %s -- see docs/data-layout.md"
                         % (what, build, path))
    return path


def genome_path(build):
    return _resolve(GENOMES, build, "reference genome")


def annotation_path(build):
    return _resolve(ANNOTATIONS, build, "GENCODE annotation .db")


def load_genome(build):
    """chrom_dict via pegg's own loader.

    prime.genome_loader returns a 2-TUPLE (chrom_dict, filtered_descriptions)
    and drops alternate, unplaced, unlocalized, patch and mitochondrial
    records. h2m.genome_loader returns (records, index_list) instead -- a
    different shape, so the two are not interchangeable; use
    chrom_dict_from_h2m() for the H2M form.
    """
    from pegg import prime
    chrom_dict, filtered = prime.genome_loader(genome_path(build))
    print("chromosomes loaded: %d (%s)" % (len(chrom_dict), build))
    return chrom_dict, filtered


def load_annotation(build):
    """gffutils.FeatureDB for the GENCODE .db of this build."""
    import gffutils
    return gffutils.FeatureDB(annotation_path(build))


def chrom_dict_from_h2m(records, index_list):
    """chrom_dict from records H2M already loaded, so the genome is read once.

    Lifted verbatim from ch_aml_library_design.ipynb cell 11. Keys are int for
    1-22 and the strings 'X'/'Y', with no 'chr' prefix -- that is exactly how
    pegg indexes chrom_dict.
    """
    chrom_dict = {}
    for idx in index_list:
        record = records[idx]
        found = re.search(r'chromosome (\w+)', record.description)
        if not found:
            continue                             # scaffolds and the mitochondrion
        name = found.group(1)
        chrom_dict[name if name in ('X', 'Y') else int(name)] = record.seq

    print('chromosomes loaded:', len(chrom_dict))
    return chrom_dict


def maf_from_h2m(df, species="m"):
    """Rename H2M's columns into MAF names, into a FRESH table.

    Lifted from cell 13. Keeps only rows H2M actually modelled (start is not
    NaN) and int-casts the positions, which arrive as float whenever any row
    was unmodelled.
    """
    suffix = "_" + species
    out = df[df["start" + suffix].notna()].copy()
    out = out.rename(columns=dict(
        (src + suffix, dst) for src, dst in H2M_RENAME.items()))
    out["Start_Position"] = out["Start_Position"].astype(int)
    out["End_Position"] = out["End_Position"].astype(int)
    return out


def normalize_chromosome(df):
    """Strip any 'chr' prefix and int-key 1-22, leaving X/Y/MT as strings.

    Lifted from cell 13. chrom_dict is keyed that way, so a 'chr7' or a string
    '7' here becomes a KeyError deep inside the per-mutation design loop --
    minutes after the genome and annotation have already been paid for.
    """
    df = df.copy()
    df['Chromosome'] = (df['Chromosome']
                        .astype(str).str.replace('^chr', '', regex=True))
    df['Chromosome'] = df['Chromosome'].apply(
        lambda c: c if c in ('X', 'Y', 'MT', 'M') else int(c))
    return df


def assert_chroms_present(df, chrom_dict):
    """Fail here, not with a KeyError hours in. Lifted from cell 13."""
    missing = sorted(set(df['Chromosome']) - set(chrom_dict), key=str)
    if missing:
        raise ValueError('not in chrom_dict: %s -- is this the right genome build?'
                         % missing)
    print('chromosomes in the design input:',
          sorted(set(df['Chromosome']), key=str))
    return df


def assert_maf_columns(df):
    missing = [c for c in MAF_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(
            "input is missing MAF columns %s. For 'cBioPortal' input pegg "
            "needs %s. If this came from H2M, run maf_from_h2m() first."
            % (missing, list(MAF_COLUMNS)))
    return df


def annotate_frame(df, db, tx_column=None, gene_column="Hugo_Symbol"):
    """Attach the reading frame each row's bystanders will be designed in.

    Lifted from cell 13. Two routes:

      tx_column given  -> bystander.add_cds_to_variants(df, db, tx_column=...)
          one transcript per row, the one already chosen upstream. Use
          tx_id_m for a mouse design off H2M output. NOTE the pegg default is
          tx_column='tx_id_h' (bystander.py:1340) -- on a mouse-only table
          that is a KeyError, and on an H2M table carrying both columns it
          silently annotates the HUMAN frame. Always pass it explicitly.

      tx_column None   -> bystander.cds_for_variants(df, db), canonical
          transcript per gene_column.

    Then cds_annotation_report() BEFORE designing: it catches a CDS that is not
    a whole number of codons, overlapping blocks, transcripts with no
    annotation, and variants outside their own CDS -- each of which otherwise
    just produces zero bystanders with no error to explain why.
    """
    from pegg import bystander

    if tx_column is not None:
        if tx_column not in df.columns:
            raise SystemExit("tx_column %r is not in the input; columns are %s"
                             % (tx_column, list(df.columns)))
        df, gene_cds = bystander.add_cds_to_variants(df, db, tx_column=tx_column)
        report_key = tx_column
    else:
        df, gene_cds = bystander.cds_for_variants(df, db, gene_column=gene_column)
        report_key = gene_column

    print(bystander.cds_annotation_report(gene_cds, mutations=df,
                                          gene_column=report_key))
    return df, gene_cds, report_key


def prepare(df, chrom_dict, db=None, tx_column=None, gene_column="Hugo_Symbol",
            require_valid_cds=False):
    """The whole cell-13 sequence, in the order the notebook runs it.

    Order matters: frame annotation first (so a bad transcript is reported
    against the names the input actually carries), then chromosome
    normalization, then the presence assert, then reset_index -- pegg's
    mutation_idx IS the input row index, so a non-contiguous index makes every
    downstream merge-back ambiguous.
    """
    assert_maf_columns(df)
    gene_cds, report_key = None, gene_column
    if db is not None:
        df, gene_cds, report_key = annotate_frame(
            df, db, tx_column=tx_column, gene_column=gene_column)
        if require_valid_cds and "cds_valid" in df.columns:
            bad = df[~df["cds_valid"].astype(bool)]
            if len(bad):
                raise SystemExit(
                    "%d of %d rows have cds_valid=False; silent bystanders "
                    "would be designed in the wrong frame or not at all. "
                    "Rows: %s" % (len(bad), len(df),
                                  bad[gene_column].tolist()[:10]))
    df = normalize_chromosome(df)
    assert_chroms_present(df, chrom_dict)
    df = df.reset_index(drop=True)
    return df, gene_cds, report_key


def main():
    ap = argparse.ArgumentParser(
        description="Prepare a variant table for pegg and report what it "
                    "would design against. Writes the prepared table and "
                    "stops -- design.py does the designing.")
    ap.add_argument("variants", help="CSV: MAF columns, or raw H2M output")
    ap.add_argument("-o", "--out", required=True, help="prepared CSV to write")
    ap.add_argument("--build", required=True, choices=sorted(
        set(list(GENOMES) + list(ANNOTATIONS))),
        help="genome build of the variant coordinates")
    ap.add_argument("--from-h2m", metavar="SPECIES", choices=["h", "m"],
                    help="rename H2M columns first ('m' = the mouse orthologue)")
    ap.add_argument("--tx-column", default=None,
                    help="transcript-id column to take the reading frame from "
                         "(e.g. tx_id_m). Omit to use the canonical transcript "
                         "per gene.")
    ap.add_argument("--gene-column", default="Hugo_Symbol")
    ap.add_argument("--no-frame", action="store_true",
                    help="skip CDS annotation (only if you will not use "
                         "silent_bystander=True)")
    ap.add_argument("--require-valid-cds", action="store_true",
                    help="fail unless every row has cds_valid=True")
    args = ap.parse_args()

    df = pd.read_csv(args.variants)
    if args.from_h2m:
        df = maf_from_h2m(df, species=args.from_h2m)
        print("renamed H2M columns -> MAF (%d rows modelled)" % len(df))

    chrom_dict, _ = load_genome(args.build)
    db = None if args.no_frame else load_annotation(args.build)

    df, _gene_cds, _key = prepare(
        df, chrom_dict, db=db, tx_column=args.tx_column,
        gene_column=args.gene_column,
        require_valid_cds=args.require_valid_cds)

    df.to_csv(args.out, index=False)
    print("wrote %s (%d rows, %d cols)" % (args.out, df.shape[0], df.shape[1]))
    print(df["Variant_Type"].value_counts().to_string())


if __name__ == "__main__":
    main()
