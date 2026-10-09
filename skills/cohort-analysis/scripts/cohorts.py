"""Public cohort resolution, loading, disease restriction and variant collapsing.

Lifted verbatim from:
  - early-leukemia-evolution/ch_aml_mutation_list.ipynb cell 5
      the four read_csv calls (sep='\t', comment='#', na_values='Not Applicable')
  - early-leukemia-evolution/ch_aml_mutation_list.ipynb cell 7
      the CANCER_TYPE_DETAILED == 'Acute Myeloid Leukemia' restriction
  - early-leukemia-evolution/ch_aml_mutation_list.ipynb cell 13
      VARIANT_KEY, LOF_CLASSES, filter_variants (with its full docstring)
  - early-leukemia-evolution/ch_aml_mutation_list.ipynb cell 14
      select_variants, generalised to take (primary_df, fallback_df) instead of
      closing over the module-level df_mutation_heme / df_mutation_genie
  - hsc-ch-mutation/.../CH.ipynb cell 28
      t_dict, cosmic_color_dict, get_cosmic_role

Do not edit the lifted bodies here without editing there -- these must not drift.
Added here and NOT in the notebooks: the $LAB_DATA_ROOT resolver, the cohort
registry, patient-level collapsing, the frequency-with-denominator table, and
the Cancer Gene Census provisioning check.
"""

from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# $LAB_DATA_ROOT  (never a hardcoded path; never an f-string with nested
# quotes -- f'{os.environ["X"]}/y' is a SyntaxError before Python 3.11 and
# the pegg env is 3.9/3.10)
# --------------------------------------------------------------------------

def lab_data_root():
    """$LAB_DATA_ROOT, or a message naming the export line."""
    root = os.environ.get("LAB_DATA_ROOT")
    if not root:
        raise SystemExit(
            "LAB_DATA_ROOT is not set -- see docs/data-layout.md\n"
            "    export LAB_DATA_ROOT=~/Documents/GitHub/Database")
    if not os.path.isdir(root):
        raise SystemExit("LAB_DATA_ROOT points at " + root + ", which does not exist")
    return root


# --------------------------------------------------------------------------
# cohort registry
# --------------------------------------------------------------------------
# Only v17.0 is exercised by the source notebooks; the other GENIE releases are
# on disk but their column stability is unverified, so they are listed and
# warned about rather than advertised.
_GENIE_VERIFIED = "v17.0"

COHORTS = {
    "genie": {
        "dir": "PublicDatabase/AACR-GENIE/{version}",
        "mutations": "data_mutations_extended.txt",
        "sample": "data_clinical_sample.txt",
        "patient": "data_clinical_patient.txt",
        "sv": "data_sv.txt",
        "default_version": _GENIE_VERIFIED,
        "versions": ["v13.1", "v14.0", "v15.0", "v16.1", "v17.0", "v18.0"],
        "build": "GRCh37",
        "controlled": True,
    },
    "heme": {
        "dir": "PublicDatabase/MSK-HEME/heme_msk_impact_2022",
        "mutations": "data_mutations.txt",
        "sample": "data_clinical_sample.txt",
        "patient": "data_clinical_patient.txt",
        "sv": None,
        "default_version": None,
        "versions": [None],
        "build": "GRCh37",
        "controlled": False,
    },
    "impact": {
        "dir": "PublicDatabase/MSK-IMPACT/msk_impact_2017",
        "mutations": "data_mutations.txt",
        "sample": "data_clinical_sample.txt",
        "patient": "data_clinical_patient.txt",
        "sv": None,
        "default_version": None,
        "versions": [None],
        "build": "GRCh37",
        "controlled": False,
    },
}


def resolve_cohort(name, version=None):
    """Paths for a registered cohort under $LAB_DATA_ROOT.

    Returns a dict with 'name', 'version', 'build', 'mutations', 'sample',
    'patient', 'sv' (None where the cohort has no SV table).
    """
    key = name.lower()
    if key not in COHORTS:
        raise SystemExit(
            "unknown cohort " + repr(name) + "; known: " + ", ".join(sorted(COHORTS)))
    spec = COHORTS[key]
    ver = version or spec["default_version"]
    if ver is not None and ver not in spec["versions"]:
        raise SystemExit(
            "cohort " + key + " has no version " + repr(ver) +
            "; known: " + ", ".join([str(v) for v in spec["versions"]]))
    if key == "genie" and ver != _GENIE_VERIFIED:
        print("WARNING: only GENIE " + _GENIE_VERIFIED + " is exercised by the "
              "source notebooks; column stability in " + str(ver) + " is unverified")

    sub = spec["dir"].format(version=ver) if ver else spec["dir"]
    base = os.path.join(lab_data_root(), sub)
    if not os.path.isdir(base):
        raise SystemExit(
            "cohort directory not found: " + base + "\n"
            "    expected under $LAB_DATA_ROOT -- see docs/data-layout.md")

    out = {"name": key, "version": ver, "build": spec["build"],
           "controlled": spec["controlled"], "dir": base}
    for field in ("mutations", "sample", "patient", "sv"):
        fname = spec[field]
        out[field] = os.path.join(base, fname) if fname else None
    return out


def _read_cbio(path, usecols=None):
    """cBioPortal / MAF table read, exactly as the notebooks read it.

    comment='#' is what transparently skips cBioPortal's four '#' metadata
    lines in the clinical files, MSK-HEME's two ('#genome_nexus_version',
    '#isoform') in its mutations file, and nothing at all in GENIE's mutations
    file whose real header is line 1. A skiprows= hardcoding handles only one
    of the three.
    """
    return pd.read_csv(path, header=0, sep="\t", comment="#",
                       na_values="Not Applicable", usecols=usecols,
                       low_memory=False)


def load_cohort(name, version=None, genes=None, usecols=None):
    """Load a cohort's mutation and sample tables.

    genes restricts the MAF to those Hugo_Symbol values while reading in
    chunks, which is what keeps GENIE v17's ~3 GB data_mutations_extended.txt
    off the heap. Without genes the whole table is read, as the notebooks do.
    """
    paths = resolve_cohort(name, version)
    if paths["controlled"]:
        print("NOTE: " + paths["name"] + " is controlled-access under a data use "
              "agreement. Read it from $LAB_DATA_ROOT; never copy it into a repo.")

    if genes:
        wanted = set(genes)
        chunks = []
        reader = pd.read_csv(paths["mutations"], header=0, sep="\t", comment="#",
                             na_values="Not Applicable", usecols=usecols,
                             chunksize=500000, low_memory=False)
        for chunk in reader:
            chunks.append(chunk[chunk["Hugo_Symbol"].isin(wanted)])
        mut = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    else:
        mut = _read_cbio(paths["mutations"], usecols=usecols)

    samp = _read_cbio(paths["sample"])
    return mut.reset_index(drop=True), samp, paths


def restrict_to_disease(mut, samp, cancer_type=None, cancer_type_detailed=None):
    """Subset a MAF to one disease, by name, and report the denominator.

    Addressing by column NAME is not optional: GENIE's data_clinical_sample.txt
    starts PATIENT_ID, SAMPLE_ID, ... with CANCER_TYPE at index 6, while
    MSK-HEME's starts SAMPLE_ID, PATIENT_ID, ... with CANCER_TYPE_DETAILED at
    index 10. Index-based code written against GENIE silently returns the wrong
    rows on MSK-HEME.

    CANCER_TYPE and CANCER_TYPE_DETAILED are NOT interchangeable. On GENIE v17,
    CANCER_TYPE == 'B-Lymphoblastic Leukemia/Lymphoma' gives 1459 samples while
    an exact CANCER_TYPE_DETAILED match gives 886, because the detailed field
    splits B-ALL into ~10 subtypes. Which column was used is printed and
    returned, never left implicit.
    """
    if cancer_type is None and cancer_type_detailed is None:
        return mut.reset_index(drop=True), {
            "column": None, "value": None,
            "n_samples": int(samp["SAMPLE_ID"].nunique()),
            "n_samples_in_maf": int(mut["Tumor_Sample_Barcode"].nunique()),
        }
    if cancer_type is not None and cancer_type_detailed is not None:
        raise SystemExit("pass cancer_type OR cancer_type_detailed, not both -- "
                         "they select different sample sets")

    col = "CANCER_TYPE" if cancer_type is not None else "CANCER_TYPE_DETAILED"
    val = cancer_type if cancer_type is not None else cancer_type_detailed
    wanted = [val] if isinstance(val, str) else list(val)
    if col not in samp.columns:
        raise SystemExit(col + " is not a column of this cohort's data_clinical_sample.txt")

    ids = samp[samp[col].isin(wanted)]["SAMPLE_ID"].unique()
    sub = mut[mut["Tumor_Sample_Barcode"].isin(ids)].reset_index(drop=True)
    prov = {"column": col, "value": wanted, "n_samples": int(len(ids)),
            "n_samples_in_maf": int(sub["Tumor_Sample_Barcode"].nunique())}
    print("filtered on %s in %r: %d samples in the cohort, %d of them carry a "
          "mutation in the loaded MAF" % (col, wanted, prov["n_samples"],
                                          prov["n_samples_in_maf"]))
    return sub, prov


# --------------------------------------------------------------------------
# patient vs sample
# --------------------------------------------------------------------------

def patient_id(sample_id):
    """GENIE sample barcode -> patient id.

    GENIE ids are GENIE-[CENTER]-[PATIENT]-[SAMPLE], so the patient is the
    first three dash-separated fields. One patient commonly has a primary and a
    relapse sample, which is why a sample-level count is not a patient-level
    count.
    """
    if not isinstance(sample_id, str):
        return sample_id
    parts = sample_id.split("-")
    if len(parts) >= 4:
        return "-".join(parts[:3])
    return sample_id


def add_patient_col(mut, sample_col="Tumor_Sample_Barcode", out_col="Patient_ID"):
    out = mut.copy()
    out[out_col] = out[sample_col].map(patient_id)
    return out


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- ch_aml_mutation_list.ipynb cell 13
# --------------------------------------------------------------------------

VARIANT_KEY = ['Chromosome', 'Start_Position', 'End_Position',
               'Reference_Allele', 'Tumor_Seq_Allele2']

LOF_CLASSES = ['Nonsense_Mutation', 'Frame_Shift_Del',
               'Frame_Shift_Ins', 'Splice_Site']


def filter_variants(df, classifications=None, protein_change=None,
                    key=None, count_col='count', per_sample=True,
                    sample_col='Tumor_Sample_Barcode'):
    """
    Filters a MAF-style mutation table and collapses recurrent variants into one
    row each, carrying a count of how often the variant was observed.

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        MAF-style table, e.g. one gene's rows pulled from cBioPortal.

    classifications
        *type = list or None*

        Variant_Classification values to keep. Default = None (no filter).

    protein_change
        *type = str, list or None*

        HGVSp_Short value(s) to keep, e.g. 'p.R882H'. Default = None.

    key
        *type = list or None*

        Columns defining variant identity. Default = None (VARIANT_KEY).

    count_col
        *type = str*

        Name of the appended count column. Default = 'count'.

    per_sample
        *type = bool*

        If True, count distinct samples rather than rows, so that a sample
        sequenced twice does not inflate the count. Default = True.

    sample_col
        *type = str*

        Column holding the sample identifier. Default = 'Tumor_Sample_Barcode'.
    """
    out = df

    if classifications is not None:
        out = out[out['Variant_Classification'].isin(classifications)]

    if protein_change is not None:
        wanted = [protein_change] if isinstance(protein_change, str) \
            else protein_change
        out = out[out['HGVSp_Short'].isin(wanted)]

    key = [k for k in (key or VARIANT_KEY) if k in out.columns]
    if len(key) == 0:
        raise ValueError('none of the key columns are present in df')

    if len(out) == 0:
        return out.assign(**{count_col: pd.Series(dtype=int)}).reset_index(drop=True)

    if per_sample and sample_col in out.columns:
        counts = (out.groupby(key, dropna=False)[sample_col]
                  .nunique().rename(count_col).reset_index())
    else:
        counts = (out.groupby(key, dropna=False)
                  .size().rename(count_col).reset_index())

    out = out.drop_duplicates(subset=key).merge(counts, on=key, how='left')

    return out.sort_values(count_col, ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------
# LIFTED from ch_aml_mutation_list.ipynb cell 14, generalised to take the two
# frames as arguments instead of closing over module-level globals.
# --------------------------------------------------------------------------

def select_variants(gene, wanted, primary_df, fallback_df=None,
                    classifications=None, primary_name="primary",
                    fallback_name="GENIE", sample_col="Tumor_Sample_Barcode"):
    """Rows for `wanted` HGVSp_Short in `gene`, from primary, else fallback.

    Some variants asked for by name are absent from the disease cohort even
    though they are common elsewhere (KRAS G12C: 0 MSK-HEME AML, 4992 GENIE;
    TP53 W146C: 0 / 2). Selecting them from the disease cohort alone would
    silently yield an empty table and drop them from the library, so fall back
    to pan-cancer for any the disease cohort does not carry -- and say which
    came from where.

    The fallback is a RECURRENCE CHECK, not a frequency source: a count from
    pan-cancer GENIE is out of a different, far larger denominator. The
    source_cohort column records which one each row came from.
    """
    out = []
    for v in wanted:
        hit = filter_variants(
            primary_df[primary_df['Hugo_Symbol'] == gene],
            classifications=classifications, protein_change=v,
            sample_col=sample_col)
        src = primary_name
        if len(hit) == 0 and fallback_df is not None:
            hit = filter_variants(
                fallback_df[fallback_df['Hugo_Symbol'] == gene],
                classifications=classifications, protein_change=v,
                sample_col=sample_col)
            src = fallback_name
        if len(hit) == 0:
            print('   %-6s %-10s NOT FOUND in either cohort' % (gene, v))
            continue
        row = hit.iloc[[0]].copy()
        row['source_cohort'] = src
        out.append(row)
        print('   %-6s %-10s from %-5s (n=%d)' % (gene, v, src, hit.iloc[0]['count']))
    if not out:
        return pd.DataFrame()
    return pd.concat(out, ignore_index=True)


def normalize_protein_change(values):
    """Accept 'R882H' as well as 'p.R882H'.

    HGVSp_Short matching is exact string matching and the strings are not
    canonical: GENIE writes 'p.R882H' but a user types 'R882H'. Frameshifts
    carry their new stop ('p.G646Wfs*12'), so a prefix like 'p.G646Wfs'
    matches nothing -- this only adds/strips the 'p.' , it cannot repair a
    truncated frameshift string.
    """
    wanted = [values] if isinstance(values, str) else list(values)
    out = []
    for v in wanted:
        v = str(v).strip()
        out.append(v)
        if v.startswith("p."):
            out.append(v[2:])
        else:
            out.append("p." + v)
    return list(dict.fromkeys(out))


# --------------------------------------------------------------------------
# frequency WITH its denominator
# --------------------------------------------------------------------------

def _wilson(k, n, z=1.959963984540054):
    """Wilson score interval -- behaves at k=0 and k=n, unlike Wald."""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = float(k) / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / d
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def variant_frequency(mut, n_cohort, per="patient",
                      sample_col="Tumor_Sample_Barcode",
                      classifications=None, protein_change=None, key=None):
    """Per-variant frequency table carrying its own denominator.

    A count without a denominator is not a frequency: the same variant is
    '3 patients' in MSK-HEME AML (296 samples) and '68 patients' in pan-cancer
    GENIE (~200k). n_cohort, freq and the 95% CI are COLUMNS so they survive
    into whatever table or figure a reviewer reads.

    per='patient' collapses GENIE's primary+relapse samples to one patient and
    is the default; per='sample' reproduces the notebooks' sample-level nunique.
    """
    if per not in ("patient", "sample"):
        raise SystemExit("per must be 'patient' or 'sample'")
    work = mut
    col = sample_col
    if per == "patient":
        work = add_patient_col(mut, sample_col=sample_col)
        col = "Patient_ID"

    out = filter_variants(work, classifications=classifications,
                          protein_change=protein_change, key=key,
                          count_col="n_" + per, per_sample=True, sample_col=col)
    if len(out) == 0:
        return out
    counts = out["n_" + per].astype(float).values
    out["n_cohort"] = int(n_cohort)
    out["freq"] = counts / float(n_cohort) if n_cohort else np.nan
    lo, hi = [], []
    for k in counts:
        a, b = _wilson(k, n_cohort)
        lo.append(a)
        hi.append(b)
    out["ci_low"] = lo
    out["ci_high"] = hi
    out["counted_per"] = per
    return out


def gene_frequency(mut, n_cohort, per="patient",
                   sample_col="Tumor_Sample_Barcode", cosmic=None):
    """Genes ranked by mutated patients, with denominator and COSMIC role."""
    work = mut
    col = sample_col
    if per == "patient":
        work = add_patient_col(mut, sample_col=sample_col)
        col = "Patient_ID"
    out = (work.groupby("Hugo_Symbol")[col].nunique()
           .sort_values(ascending=False).rename("n_" + per).reset_index())
    out["n_cohort"] = int(n_cohort)
    out["freq"] = out["n_" + per] / float(n_cohort) if n_cohort else np.nan
    out["counted_per"] = per
    if cosmic is not None:
        out["cosmic_role"] = out["Hugo_Symbol"].map(
            lambda g: get_cosmic_role(g, cosmic))
    return out


# --------------------------------------------------------------------------
# protein position
# LIFTED VERBATIM -- ch_aml_mutation_list.ipynb cell 22
# --------------------------------------------------------------------------

PROTEIN_AA = {
    'DNMT3A': 912, 'TET2': 2002, 'ASXL1': 1541, 'IDH1': 414, 'IDH2': 452,
    'SF3B1': 1304, 'SRSF2': 221, 'U2AF1': 240, 'JAK2': 1132, 'TP53': 393,
    'PPM1D': 605, 'FLT3': 993, 'NRAS': 189, 'KRAS': 189, 'KIT': 976,
    'PTPN11': 593, 'NPM1': 294, 'RUNX1': 480, 'WT1': 517,
}


def aa_position(hgvsp):
    """First residue number in an HGVSp_Short string, or None.

    Handles 'p.R882H', 'p.R771*', 'p.Q527fs' and the frameshift/indel forms; the
    number is all that matters here, not what the change is.
    """
    if not isinstance(hgvsp, str):
        return None
    found = re.search(r'(\d+)', hgvsp)
    return int(found.group(1)) if found else None


def residue_counts(df, gene, classifications=None, protein_len=None,
                   sample_col='Tumor_Sample_Barcode'):
    """Per-residue mutation counts, counted once per sample.

    Lifted from `_counts` in cell 22; protein_len is a parameter here because
    PROTEIN_AA covers only the 19 CH/AML genes the notebook plotted, and a
    KeyError on gene 20 is a worse failure than an unbounded axis.
    """
    sub = df[df['Hugo_Symbol'] == gene].copy()
    if classifications is not None:
        sub = sub[sub['Variant_Classification'].isin(classifications)]
    sub['aa_pos'] = sub['HGVSp_Short'].map(aa_position)
    length = protein_len or PROTEIN_AA.get(gene)
    if length is None:
        sub = sub[sub['aa_pos'].notna()]
    else:
        sub = sub[sub['aa_pos'].notna() & sub['aa_pos'].between(1, length)]
    # one row per sample per change, so a sample sequenced twice is not counted twice
    if sample_col in sub.columns:
        sub = sub.drop_duplicates(subset=[sample_col, 'HGVSp_Short'])
    return sub.groupby('aa_pos').size()


# --------------------------------------------------------------------------
# COSMIC Cancer Gene Census
# LIFTED VERBATIM (t_dict, cosmic_color_dict, get_cosmic_role) -- CH.ipynb cell 28
# The hardcoded path into B-ALL-in-vivo-base-editing is NOT lifted.
# --------------------------------------------------------------------------

t_dict = {
    'TSG': 'TSG', 'TSG, fusion': 'TSG',
    'Undefined': 'Undefined', 'fusion': 'Undefined',
    'oncogene': 'Oncogene',
    'oncogene, TSG': 'Oncogene/TSG', 'oncogene, TSG, fusion': 'Oncogene/TSG',
    'oncogene, fusion': 'Oncogene',
}

cosmic_color_dict = {'Oncogene': 'tab:red', 'Oncogene/TSG': 'tab:purple',
                     'TSG': 'tab:blue', 'Undefined': 'tab:grey'}

_CENSUS_NAMES = ["Census_all.csv", "cancer_gene_census.csv",
                 "Cancer_Gene_Census.csv"]


def load_census(path=None):
    """Cancer Gene Census, from $LAB_DATA_ROOT/PublicDatabase/COSMIC/.

    This file is genuinely not provisioned yet: that directory currently holds
    only Cosmic_MutantCensus_Tsv_v99_GRCh37. Rather than fail on a missing
    path, say what to download. The census read is .fillna('Undefined') to
    match the notebook.
    """
    if path:
        return pd.read_csv(path).fillna("Undefined")
    cdir = os.path.join(lab_data_root(), "PublicDatabase", "COSMIC")
    for name in _CENSUS_NAMES:
        p = os.path.join(cdir, name)
        if os.path.exists(p):
            return pd.read_csv(p).fillna("Undefined")
    raise SystemExit(
        "Cancer Gene Census not found under " + cdir + "\n"
        "    Download it from cancer.sanger.ac.uk (free academic login),\n"
        "    'Cancer Gene Census' -> Census_all.csv, and save it as\n"
        "    " + os.path.join(cdir, "Census_all.csv") + "\n"
        "    It is a DIFFERENT table from Cosmic_MutantCensus_v99_GRCh37.tsv,\n"
        "    which has UPPER_SNAKE columns and no 'Role in Cancer'.\n"
        "    Or pass --census /path/to/Census_all.csv explicitly.")


def get_cosmic_role(gene, cosmic):
    """COSMIC Role in Cancer, collapsed to the four plotting categories.

    .values[0] means a gene with several census rows gets whichever sorted
    first, and t_dict.get(..., 'Undefined') maps any role string a future
    COSMIC release adds to 'Undefined' -- the seven keys are the complete
    v2024 set, so a new release needs t_dict updated, not the default relied on.
    """
    sub = cosmic[cosmic['Gene Symbol'] == gene]
    if len(sub) == 0:
        return 'Undefined'
    raw = sub['Role in Cancer'].values[0]
    if raw not in t_dict:
        print("WARNING: COSMIC role " + repr(raw) + " for " + gene +
              " is not in t_dict -- mapped to Undefined; update t_dict")
    return t_dict.get(raw, 'Undefined')
