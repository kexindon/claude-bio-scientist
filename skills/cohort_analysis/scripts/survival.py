"""Patient-level outcome analysis on a local cohort: Cox PH, log-rank, KM.

Lifted verbatim from hsc-ch-mutation/20260301-committee-meeting-clinical-analysis/
CH.ipynb:
  - cell 4:  NUMERIC_COLUMNS
  - cell 6:  normalize_flt3, unique_sample_gene, build_pair_survival_dataset,
             fit_binary_cox, km_plot
  - cell 35: lr, stars  (log-rank helpers)
  - cell 1:  CH_FUNCTIONAL_GROUPS, CH, CH_UNIQUE, GENE_TO_GROUP
  - cell 8:  the de novo subset, Time_to_OS_years, base_surv, gene_freq,
             top20_non_ch
  - cell 9:  top_mutation_per_gene
Do not edit these bodies here without editing there -- the two must not drift.

Added here and NOT in the notebook: load_cohort_table (the pyreadr/CSV reader
with select-by-name), cohort_qc, the min_per_arm parameter in place of the
buried n>=5 constant, the 'double_vs_rest' and 'four_arm' contrasts from
Figure_2.R, the covariates passthrough, convergence reporting, and the
pairs-skipped ledger.

This module never loads, copies or embeds any cohort file. Patient-level
cohort data stays outside the repo.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- CH.ipynb cell 1
# --------------------------------------------------------------------------

CH_FUNCTIONAL_GROUPS = {
    "Epigenetic":            ["DNMT3A", "TET2", "ASXL1", "IDH1", "IDH2"],
    "Splicing":              ["SF3B1", "SRSF2", "U2AF1"],
    "Kinase signaling":      ["JAK2"],
    "DNA damage response":   ["TP53"],
}

# Flatten to the working CH list, preserving group order
CH = [g for genes in CH_FUNCTIONAL_GROUPS.values() for g in genes]
CH_UNIQUE = list(dict.fromkeys(CH))

# Reverse lookup: gene -> functional group
GENE_TO_GROUP = {
    gene: group
    for group, genes in CH_FUNCTIONAL_GROUPS.items()
    for gene in genes
}

# len(CH_UNIQUE) == 10, not 9. CH.ipynb's own markdown says "9 unique CH genes,
# so the volcano panel is a full 3x3 layout" and every grid in the code is
# plt.subplots(3, 4). Anything built from that prose makes the wrong grid.
assert len(CH_UNIQUE) == 10, "CH_UNIQUE changed; the 3x4 grid assumption moves"


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- CH.ipynb cell 4
# --------------------------------------------------------------------------

NUMERIC_COLUMNS = [
    "VAF", "Age", "Hemoglobin", "Platelet", "WBC", "LDH",
    "BM_blast_percent", "PB_blast_percent", "PB_wbc_percent",
    "Time_to_OS", "Censor", "mut_freq_gene", "mut_freq_pt", "mut_freq_bin",
]


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- CH.ipynb cell 6
# --------------------------------------------------------------------------

def normalize_flt3(df):
    df = df.copy()
    mask_itd = (df["Gene"] == "FLT3") & (df["variant_type"].isin({"ITD", "INDEL"}))
    mask_tkd = (df["Gene"] == "FLT3") & (df["variant_type"].isin({"SNV", "Deletion", "other"}))
    df.loc[mask_itd, "Gene"] = "FLT3-ITD"
    df.loc[mask_tkd, "Gene"] = "FLT3-TKD"
    return df


def unique_sample_gene(df):
    return df[["Sample", "Gene"]].dropna().drop_duplicates().copy()


def build_pair_survival_dataset(base_surv, sample_gene_df, ch_gene, driver_gene):
    ch_samples = set(sample_gene_df.loc[sample_gene_df["Gene"] == ch_gene, "Sample"])
    driver_samples = set(sample_gene_df.loc[sample_gene_df["Gene"] == driver_gene, "Sample"])

    double_samples = ch_samples & driver_samples
    driver_only_samples = driver_samples - ch_samples

    temp = base_surv.copy()
    temp["group"] = np.where(
        temp["Sample"].isin(double_samples),
        "Double",
        np.where(temp["Sample"].isin(driver_only_samples), "Driver_only", np.nan),
    )
    temp = temp.dropna(subset=["group"]).copy()
    temp["group_binary"] = (temp["group"] == "Double").astype(int)
    return temp


def fit_binary_cox(df, duration_col="Time_to_OS_years", event_col="Censor", group_col="group_binary"):
    temp = df[[duration_col, event_col, group_col]].dropna().copy()
    if temp.empty or temp[group_col].nunique() < 2:
        return None
    cph = CoxPHFitter()
    cph.fit(temp, duration_col=duration_col, event_col=event_col)
    row = cph.summary.iloc[0]
    return {
        "HR": float(np.exp(row["coef"])),
        "lower_95": float(np.exp(row["coef lower 95%"])),
        "upper_95": float(np.exp(row["coef upper 95%"])),
        "p_value": float(row["p"]),
        "model": cph,
    }


def km_plot(ax, df, title):
    kmf = KaplanMeierFitter()
    palette = {"Driver_only": "#e08214", "Double": "#B24745"}
    for group, grp in df.groupby("group"):
        kmf.fit(grp["Time_to_OS_years"], grp["Censor"], label=group)
        kmf.plot(ax=ax, ci_show=False, color=palette.get(group, "grey"))
    ax.set_xlabel("Years")
    ax.set_ylabel("Survival probability")
    ax.set_title(title)
    return ax


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- CH.ipynb cell 35
# --------------------------------------------------------------------------

def lr(a, b):
    r = logrank_test(a["Time_to_OS_years"], b["Time_to_OS_years"],
                     event_observed_A=a["Censor"], event_observed_B=b["Censor"])
    return r.p_value


def stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- CH.ipynb cell 9
# --------------------------------------------------------------------------

def top_mutation_per_gene(gene_list, df):
    rows = []
    for gene in gene_list:
        sub = df[df["Gene"] == gene]
        if sub.empty:
            continue
        freq = sub.groupby("amino_acid_change")["Sample"].nunique().sort_values(ascending=False)
        rows.append({
            "Gene": gene,
            "amino_acid_change": freq.index[0],
            "n_patients": int(freq.iloc[0]),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# loading -- wrapper around the notebook's cells 3/4/8
# --------------------------------------------------------------------------

def load_cohort_table(path, object_name=None, subset=None):
    """Load a cohort to canonical long format from .RData/.rds or CSV/TSV.

    The file stays where it is. Nothing is copied into the repo and no path is
    hardcoded -- the caller passes one, and it should live outside any git
    working tree.

    pyreadr.read_r returns an OrderedDict of every object in the .RData, hence
    the notebook's next(iter(result.values())). If the file ever gains a second
    object that picks an arbitrary one, so select by name when you can.
    """
    low = path.lower()
    if low.endswith(".rdata") or low.endswith(".rda") or low.endswith(".rds"):
        try:
            import pyreadr
        except ImportError as exc:
            raise SystemExit(
                "pyreadr is required to read " + path + "\n"
                "    conda install -c conda-forge pyreadr   (or pip install pyreadr)\n"
                "    note pyreadr reads data frames only, not arbitrary R objects"
            ) from exc
        result = pyreadr.read_r(path)
        if object_name is not None:
            if object_name not in result:
                raise SystemExit(
                    "no object " + repr(object_name) + " in " + path +
                    "; found: " + ", ".join([str(k) for k in result.keys()]))
            df = result[object_name]
        elif len(result) == 1:
            df = next(iter(result.values()))
        else:
            names = ", ".join([str(k) for k in result.keys()])
            raise SystemExit(path + " holds several objects (" + names +
                             "); pass --rdata-object to choose one")
    elif low.endswith(".tsv") or low.endswith(".txt"):
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)

    df = df.copy()
    for col in [c for c in NUMERIC_COLUMNS if c in df.columns]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    required = ["Sample", "Gene"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit("cohort table is missing required column(s): " +
                         ", ".join(missing))

    if subset is not None:
        if "Subset" not in df.columns:
            raise SystemExit("--subset given but the table has no Subset column")
        before = len(df)
        df = df.loc[df["Subset"] == subset].copy()
        print("Subset == %r: %d of %d rows, %d unique samples"
              % (subset, len(df), before, df["Sample"].nunique()))

    if "variant_type" in df.columns and (df["Gene"] == "FLT3").any():
        df = normalize_flt3(df)
        leftover = df.loc[df["Gene"] == "FLT3"]
        if len(leftover):
            kinds = sorted(set(leftover["variant_type"].dropna().astype(str)))
            print("WARNING: %d FLT3 rows did not split into ITD/TKD "
                  "(variant_type %s) and keep the bare symbol FLT3. "
                  "Three FLT3 entities now compete in the frequency ranking, "
                  "and the unsplit one can displace a real driver out of a "
                  "top-N list." % (len(leftover), kinds))

    if "Time_to_OS" in df.columns:
        df["Time_to_OS_years"] = pd.to_numeric(df["Time_to_OS"], errors="coerce") / 365.0
    return df.reset_index(drop=True)


def base_survival(df):
    """One row per patient: Sample, Time_to_OS_years, Censor."""
    for col in ("Time_to_OS_years", "Censor"):
        if col not in df.columns:
            raise SystemExit("missing " + col + "; Time_to_OS (days) and Censor "
                             "are required for any survival analysis")
    return df[["Sample", "Time_to_OS_years", "Censor"]].drop_duplicates().copy()


def cohort_qc(df):
    """One-page sanity report to read BEFORE trusting any HR or KM curve.

    The two errors this catches are the ones that reach a figure silently:
    Time_to_OS is in DAYS and nothing in the column name says so, and Censor's
    polarity is a convention. A Cox model fits happily on days and returns an
    HR per day that looks fine and is comparable to nothing published.
    """
    rep = {}
    rep["n_rows"] = int(len(df))
    rep["n_patients"] = int(df["Sample"].nunique())
    rep["n_genes"] = int(df["Gene"].nunique()) if "Gene" in df.columns else None

    if "Subset" in df.columns:
        rep["subset_levels"] = df["Subset"].value_counts().to_dict()

    if "Time_to_OS" in df.columns:
        t = pd.to_numeric(df["Time_to_OS"], errors="coerce")
        rep["Time_to_OS_max"] = float(np.nanmax(t.values)) if t.notna().any() else None
        if rep["Time_to_OS_max"] is not None:
            if rep["Time_to_OS_max"] < 30:
                rep["time_unit_verdict"] = (
                    "max is %.1f -- this looks like YEARS already. Dividing by "
                    "365 again would be wrong." % rep["Time_to_OS_max"])
            elif rep["Time_to_OS_max"] > 400:
                rep["time_unit_verdict"] = (
                    "max is %.0f -- consistent with DAYS (MetaAML max is 6780); "
                    "Time_to_OS_years = Time_to_OS / 365." % rep["Time_to_OS_max"])
            else:
                rep["time_unit_verdict"] = (
                    "max is %.0f -- ambiguous between days and months. Check the "
                    "source before fitting anything." % rep["Time_to_OS_max"])

    if "Censor" in df.columns:
        c = pd.to_numeric(df["Censor"], errors="coerce")
        levels = sorted(set(c.dropna().unique().tolist()))
        rep["Censor_levels"] = levels
        if not set(levels).issubset({0.0, 1.0}):
            rep["censor_verdict"] = ("Censor holds %r -- lifelines wants 1 = event "
                                     "(death) observed, 0 = censored." % (levels,))
        else:
            per_pt = df[["Sample", "Censor"]].drop_duplicates()
            n_event = int(pd.to_numeric(per_pt["Censor"], errors="coerce").sum())
            rep["n_events"] = n_event
            rep["n_censored"] = int(per_pt["Sample"].nunique() - n_event)
            rep["censor_verdict"] = (
                "%d events / %d censored, i.e. 1 is being read as DEATH. If that "
                "is backwards for this cohort every HR is inverted."
                % (n_event, rep["n_censored"]))

    if "Time_to_OS_years" in df.columns:
        per_pt = df[["Sample", "Time_to_OS_years"]].drop_duplicates()
        vals = pd.to_numeric(per_pt["Time_to_OS_years"], errors="coerce").dropna()
        if len(vals):
            rep["median_followup_years"] = float(np.median(vals.values))
    return rep


def print_qc(rep):
    print("--- cohort QC ---")
    for key in ("n_rows", "n_patients", "n_genes", "n_events", "n_censored",
                "median_followup_years", "Time_to_OS_max", "Censor_levels",
                "time_unit_verdict", "censor_verdict", "subset_levels"):
        if key in rep and rep[key] is not None:
            print("  %-22s %s" % (key, rep[key]))


# --------------------------------------------------------------------------
# contrasts
# --------------------------------------------------------------------------
# CH.ipynb and the original Figure_2.R do NOT fit the same model. The notebook
# keeps Double vs Driver_only and drops wild-type and focal-only patients; the
# R code (Figure_2.R lines 574-609) builds a four-level strata by rowSums of
# STRATA1/2/3 -- 0 = wild-type, 1 = focal-only, 2 = comparator-only, 6 =
# double -- then fits coxph(Surv(Time_to_OS, Censor) ~ STRATA_all) with
# STRATA_all = ifelse(STRATA == 6, 6, 0), i.e. DOUBLE vs EVERYONE ELSE, and
# draws the 4-arm KM. Same data, different question, different HR. The R guard
# is n_pts_both >= 10; the notebook's is n >= 5 per arm.

CONTRASTS = ("double_vs_driver_only", "double_vs_rest", "four_arm")

_CONTRAST_DOC = {
    "double_vs_driver_only":
        "among comparator-mutant patients, does adding the focal mutation "
        "change outcome (wild-type and focal-only patients are DROPPED)",
    "double_vs_rest":
        "double-mutant vs every other patient in the cohort "
        "(Figure_2.R STRATA_all == 6 vs 0)",
    "four_arm":
        "wild-type / focal-only / comparator-only / double, as four arms "
        "(Figure_2.R survfit; KM only, no single HR)",
}


def build_contrast(base_surv, sample_gene_df, focal_gene, comparator_gene,
                   contrast="double_vs_driver_only"):
    """Per-patient frame with an `arm` column, for one of the three contrasts.

    `arm` is always one of Wildtype / Focal_only / Comparator_only / Double.
    group_binary is the 0/1 the Cox model fits, and which patients it compares
    depends on the contrast -- which is why the contrast name travels with
    every result row and belongs in every KM legend.
    """
    if contrast not in CONTRASTS:
        raise SystemExit("contrast must be one of " + ", ".join(CONTRASTS))

    focal = set(sample_gene_df.loc[sample_gene_df["Gene"] == focal_gene, "Sample"])
    comp = set(sample_gene_df.loc[sample_gene_df["Gene"] == comparator_gene, "Sample"])

    temp = base_surv.copy()
    in_f = temp["Sample"].isin(focal)
    in_c = temp["Sample"].isin(comp)
    temp["arm"] = np.select(
        [in_f & in_c, in_f & ~in_c, ~in_f & in_c],
        ["Double", "Focal_only", "Comparator_only"],
        default="Wildtype")

    if contrast == "double_vs_driver_only":
        temp = temp.loc[temp["arm"].isin(["Double", "Comparator_only"])].copy()
        # the notebook's own group labels, so km_plot and the published
        # figures keep working unchanged
        temp["group"] = np.where(temp["arm"] == "Double", "Double", "Driver_only")
    else:
        temp["group"] = temp["arm"]

    temp["group_binary"] = (temp["arm"] == "Double").astype(int)
    temp["contrast"] = contrast
    return temp.reset_index(drop=True)


def fit_cox(df, covariates=None, duration_col="Time_to_OS_years",
            event_col="Censor", group_col="group_binary"):
    """Cox PH on group_binary, optionally adjusted for covariates.

    The notebook's model is univariate on group_binary alone. For AML outcome
    that is not a safe default: NPM1 and FLT3-ITD are themselves among the
    usual comparators AND are the dominant prognostic confounders, and
    Age/Sex/Cytogenetics/Risk/FAB sit unused in the table. Pass covariates=
    to adjust, and report the adjusted HR beside the univariate one.

    Convergence is reported, not swallowed. CH.ipynb runs under a blanket
    warnings.filterwarnings('ignore') twice over, so a Cox fit that barely
    converged returns an HR and a CI with no visible complaint.
    """
    cols = [duration_col, event_col, group_col] + list(covariates or [])
    missing = [c for c in cols if c not in df.columns]
    if missing:
        return {"ok": False, "reason": "missing column(s) " + ", ".join(missing)}

    temp = df[cols].dropna().copy()
    for c in covariates or []:
        if temp[c].dtype == object:
            temp = pd.get_dummies(temp, columns=[c], drop_first=True)
    if temp.empty or temp[group_col].nunique() < 2:
        return {"ok": False, "reason": "fewer than 2 groups after dropna"}

    cph = CoxPHFitter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cph.fit(temp, duration_col=duration_col, event_col=event_col)
    notes = sorted(set([str(w.message).split("\n")[0][:120] for w in caught]))

    row = cph.summary.loc[group_col]
    return {
        "ok": True,
        "HR": float(np.exp(row["coef"])),
        "lower_95": float(np.exp(row["coef lower 95%"])),
        "upper_95": float(np.exp(row["coef upper 95%"])),
        "p_value": float(row["p"]),
        "n": int(len(temp)),
        "covariates": list(covariates or []),
        "fit_warnings": "; ".join(notes),
        "model": cph,
    }


def pair_survival_grid(base_surv, sample_gene_df, focal_genes, comparator_genes,
                       contrast="double_vs_driver_only", min_per_arm=5,
                       covariates=None):
    """Cox over every focal x comparator pair, plus the skipped-pair ledger.

    Returns (results, skipped). min_per_arm replaces CH.ipynb's buried
    `if n_double < 5 or n_driver_only < 5: continue`; Figure_2.R's equivalent
    guard is n_pts_both >= 10, so the two disagree on which pairs exist at all.

    A pair below the threshold is MISSING from `results`, not non-significant.
    `skipped` names why, per pair, so a later .pivot cannot pass NaN off as a
    null result.
    """
    rows = []
    skipped = []
    plot_data = {}

    for focal in focal_genes:
        for comp in comparator_genes:
            if focal == comp:
                continue
            pair_key = "%s + %s" % (focal, comp)
            temp = build_contrast(base_surv, sample_gene_df, focal, comp,
                                  contrast=contrast)
            n_double = int((temp["arm"] == "Double").sum())
            n_other = int((temp["arm"] != "Double").sum())

            if contrast == "four_arm":
                plot_data[pair_key] = temp
                arms = temp["arm"].value_counts().to_dict()
                if n_double < min_per_arm:
                    skipped.append({"focal_gene": focal, "comparator_gene": comp,
                                    "pair": pair_key, "n_double": n_double,
                                    "n_other": n_other,
                                    "reason": "n_double < min_per_arm=%d" % min_per_arm})
                    continue
                rows.append({"focal_gene": focal, "comparator_gene": comp,
                             "pair": pair_key, "contrast": contrast,
                             "n_double": n_double, "n_other": n_other,
                             "arms": arms})
                continue

            if n_double < min_per_arm or n_other < min_per_arm:
                skipped.append({"focal_gene": focal, "comparator_gene": comp,
                                "pair": pair_key, "n_double": n_double,
                                "n_other": n_other,
                                "reason": "arm below min_per_arm=%d" % min_per_arm})
                continue

            fit = fit_cox(temp, covariates=covariates)
            if not fit["ok"]:
                skipped.append({"focal_gene": focal, "comparator_gene": comp,
                                "pair": pair_key, "n_double": n_double,
                                "n_other": n_other, "reason": fit["reason"]})
                continue

            plot_data[pair_key] = temp
            rows.append({
                "focal_gene": focal, "comparator_gene": comp, "pair": pair_key,
                "contrast": contrast, "contrast_meaning": _CONTRAST_DOC[contrast],
                "HR": fit["HR"], "lower_95": fit["lower_95"],
                "upper_95": fit["upper_95"], "p_value": fit["p_value"],
                "n_double": n_double, "n_other": n_other, "n_total": fit["n"],
                "covariates": ",".join(fit["covariates"]),
                "fit_warnings": fit["fit_warnings"],
            })

    res = pd.DataFrame(rows)
    skip = pd.DataFrame(skipped)
    n_tested = len(res)
    print("%d pairs tested, %d skipped (of %d considered), contrast=%s, "
          "min_per_arm=%d" % (n_tested, len(skip), n_tested + len(skip),
                              contrast, min_per_arm))

    if len(res) and "p_value" in res.columns:
        res["q_value_bonf"] = multipletests(res["p_value"], method="bonferroni")[1]
        res["q_value_fdr"] = multipletests(res["p_value"], method="fdr_bh")[1]
        res["survival_status"] = np.select(
            [(res["p_value"] <= 0.05) & (res["HR"] < 1),
             (res["p_value"] <= 0.05) & (res["HR"] > 1)],
            ["Favorable", "Unfavorable"], default="NS")
        # CH.ipynb's own thresholds, kept explicit: survival_status is RAW p,
        # while its label / sig_pairs selection used q_value_bonf <= 0.15.
        res["significance_basis"] = "survival_status: raw p <= 0.05"
    return res.reset_index(drop=True), skip, plot_data


def km_frame(plot_df, arms=None):
    """Per-arm survival input for a KM plot: the mechanics, not the style.

    Returns a list of (label, times, events). Colours, fonts, panel size, the
    risk table and PDF export come from the sci-plot skill.
    """
    out = []
    order = arms or ["Wildtype", "Focal_only", "Comparator_only", "Double",
                     "Driver_only"]
    present = [a for a in order if a in set(plot_df["group"])]
    for a in present:
        sub = plot_df.loc[plot_df["group"] == a].dropna(
            subset=["Time_to_OS_years", "Censor"])
        out.append((a + " (n=%d)" % len(sub),
                    sub["Time_to_OS_years"].values, sub["Censor"].values))
    return out


def logrank_matrix(plot_df):
    """Pairwise log-rank p between every pair of arms present."""
    arms = sorted(set(plot_df["group"]))
    rows = []
    for i in range(len(arms)):
        for j in range(i + 1, len(arms)):
            a = plot_df.loc[plot_df["group"] == arms[i]].dropna(
                subset=["Time_to_OS_years", "Censor"])
            b = plot_df.loc[plot_df["group"] == arms[j]].dropna(
                subset=["Time_to_OS_years", "Censor"])
            if len(a) == 0 or len(b) == 0:
                continue
            p = lr(a, b)
            rows.append({"arm_a": arms[i], "arm_b": arms[j], "n_a": len(a),
                         "n_b": len(b), "logrank_p": p, "stars": stars(p)})
    return pd.DataFrame(rows)
