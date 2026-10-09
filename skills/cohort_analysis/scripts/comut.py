"""Co-occurrence / mutual-exclusivity testing: Fisher's exact + BH FDR.

Lifted verbatim from hsc-ch-mutation/20260301-committee-meeting-clinical-analysis/
co-mutation-leukemia.ipynb:
  - cell 11: build_binary_matrix, fisher_or_p_from_binary, pairwise_or_fisher,
             classify_significance
  - cell 12: p_to_stars
  - cell 13: safe_log_or
Do not edit these bodies here without editing there -- the two must not drift.

Also lifted, as a documented COMPATIBILITY SHIM only, from CH.ipynb cell 6:
  - fisher_for_pair
It is kept so published numbers can be reproduced bit-for-bit, and it is not
what the CLI calls. It substitutes the magic literal 0.06982398 for a true odds
ratio of 0; see its docstring.

Added here and NOT in the notebooks: pairwise_grid (the focal x comparator loop
with BOTH within-focal and global BH), the or_floored flag, and the
no-test-performed accounting.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection, multipletests


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- co-mutation-leukemia.ipynb cell 11
# --------------------------------------------------------------------------

def build_binary_matrix(df_long, patient_col="PatientID", gene_col="Gene"):
    """
    Build Patient x Gene binary matrix from long mutation table.
    """
    # Keep unique patient-gene pairs to avoid double counting
    df_u = df_long[[patient_col, gene_col]].drop_duplicates()
    # Assign 1 for presence
    df_u["value"] = 1
    # Pivot to wide; fill missing with 0
    mat = df_u.pivot_table(index=patient_col, columns=gene_col, values="value", fill_value=0, aggfunc="max")
    # Ensure integer 0/1
    return mat.astype(np.int8)


def fisher_or_p_from_binary(x, y):
    """
    Compute OR and Fisher p-value from two binary vectors x,y (0/1).
    Returns: (or_value, p_value, table_counts)
    table layout:
        y=1  y=0
    x=1   a    b
    x=0   c    d
    """
    # a: both 1, b: x=1,y=0, c: x=0,y=1, d: both 0
    a = int(np.sum((x == 1) & (y == 1)))
    b = int(np.sum((x == 1) & (y == 0)))
    c = int(np.sum((x == 0) & (y == 1)))
    d = int(np.sum((x == 0) & (y == 0)))

    table = np.array([[a, b], [c, d]], dtype=np.int64)

    # Fisher exact; scipy returns OR and p
    # NOTE: if b*c==0 OR can be inf or 0; that's okay but handle later for log.
    or_value, p_value = fisher_exact(table, alternative="two-sided")
    return or_value, p_value, (a, b, c, d)


def pairwise_or_fisher(
    mat, gene1, genes2, patient_subset=None, fdr_method="bh"
):
    """
    Compute OR/p/q and counts for gene1 vs each gene in genes2
    within an optional patient subset.
    """
    if patient_subset is not None:
        mat_sub = mat.loc[mat.index.intersection(patient_subset)]
    else:
        mat_sub = mat

    # Make sure gene exists
    if gene1 not in mat_sub.columns:
        raise ValueError("%s not in matrix columns" % gene1)

    x = mat_sub[gene1].values
    results = []

    for g2 in genes2:
        if g2 not in mat_sub.columns:
            continue
        y = mat_sub[g2].values

        or_value, p_value, (a, b, c, d) = fisher_or_p_from_binary(x, y)

        # s1: gene1 mutated count; s2: gene2 mutated count; sintersect: co-mut
        s1 = int(np.sum(x == 1))
        s2 = int(np.sum(y == 1))
        sintersect = a

        results.append({
            "gene1": gene1,
            "gene2": g2,
            "n": int(len(x)),
            "s1": s1,
            "s2": s2,
            "sintersect": sintersect,
            "a": a, "b": b, "c": c, "d": d,
            "odds_ratio": or_value,
            "p": p_value,
        })

    df_res = pd.DataFrame(results)

    # FDR correction
    if len(df_res) > 0:
        q = fdrcorrection(df_res["p"].values)[1]
        df_res["q"] = q
    else:
        df_res["q"] = np.array([])

    return df_res


def classify_significance(df_res, q_thresh=0.05):
    """
    Add 'significant' label like your R code:
    - q < thresh and OR > 1 -> Co-occuring
    - q < thresh and OR < 1 -> Mutually exclusive
    - else NS
    """
    df = df_res.copy()
    df["significant"] = "NS"
    mask = df["q"].notna() & (df["q"] < q_thresh) & df["odds_ratio"].notna()
    df.loc[mask & (df["odds_ratio"] > 1), "significant"] = "Co-occuring"
    df.loc[mask & (df["odds_ratio"] < 1), "significant"] = "Mutually exclusive"
    return df


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- co-mutation-leukemia.ipynb cell 13
# --------------------------------------------------------------------------

def safe_log_or(or_values, eps=1e-12, cap=50):
    """
    Compute log(OR) safely:
    - Replace 0 with eps
    - Cap inf to cap
    """
    or_arr = np.asarray(or_values, dtype=float)
    or_arr = np.where(or_arr == 0, eps, or_arr)
    log_or = np.log(or_arr)
    log_or = np.where(np.isinf(log_or) & (log_or > 0), cap, log_or)
    log_or = np.where(np.isinf(log_or) & (log_or < 0), -cap, log_or)
    return log_or


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- co-mutation-leukemia.ipynb cell 12
# --------------------------------------------------------------------------

def p_to_stars(p, levels=(0.001, 0.01, 0.1)):
    """
    Convert p/q value to significance stars.
    levels: (***, **, *) thresholds; last one can be dot.
    """
    if p is None or np.isnan(p):
        return ""
    if p < levels[0]:
        return "***"
    if p < levels[1]:
        return "**"
    if p < levels[2]:
        return "*"
    return ""


# --------------------------------------------------------------------------
# LIFTED VERBATIM -- CH.ipynb cell 6. COMPATIBILITY SHIM, not the engine.
# --------------------------------------------------------------------------

def fisher_for_pair(sample_gene_df, gene_a, gene_b):
    """Reproduces CH.ipynb's published odds ratios exactly, bugs included.

    DO NOT use this for a new analysis. It rewrites a true odds ratio of 0 to
    the literal 0.06982398 and a NaN odds ratio to 1.0, so that log(OR) plots
    without an infinity. On the real de novo MetaAML cohort that fires for 15
    of 200 focal x comparator cells, and those 15 land at log(OR) = -2.66 on
    every volcano, bubble and heatmap as if that were a measured effect size,
    when it actually means "zero double-mutant patients". The p-value is
    unaffected, so such a cell reads as a strong, significant mutual
    exclusivity at a fabricated effect size.

    It is kept only so a published figure can be regenerated byte-identically.
    pairwise_grid() uses safe_log_or and an explicit or_floored flag instead.
    """
    samples = sorted(sample_gene_df["Sample"].unique())
    a_set = set(sample_gene_df.loc[sample_gene_df["Gene"] == gene_a, "Sample"])
    b_set = set(sample_gene_df.loc[sample_gene_df["Gene"] == gene_b, "Sample"])

    a = pd.Series([s in a_set for s in samples], index=samples, dtype=int)
    b = pd.Series([s in b_set for s in samples], index=samples, dtype=int)
    table = pd.crosstab(a, b).reindex(index=[0, 1], columns=[0, 1], fill_value=0)

    odds_ratio, p_value = fisher_exact(table.values)
    if np.isnan(odds_ratio):
        odds_ratio = 1.0
    if odds_ratio == 0:
        odds_ratio = 0.06982398

    n_double = int(table.loc[1, 1])
    return {
        "odds_ratio": float(odds_ratio),
        "p_value": float(p_value),
        "n_double": n_double,
        "table": table,
    }


# --------------------------------------------------------------------------
# the focal x comparator grid -- wrapper, not lifted
# --------------------------------------------------------------------------

def pairwise_grid(mat, focal_genes, comparator_genes, patient_subset=None,
                  q_thresh=0.05, eps=1e-12, cap=50):
    """focal x comparator Fisher grid with BOTH FDR families reported.

    Returns one row per tested pair with the full a/b/c/d contingency table,
    odds_ratio, p, log_or, or_floored, and TWO q columns:

      q_within_focal  BH applied within each focal gene (one family per
                      volcano panel; what CH.ipynb cell 11 does, 20 tests)
      q_global        BH applied across the whole grid (what the original
                      Figure_2.R line 451 does, 200 tests)

    These give different q values on the same data -- the within-focal family
    is systematically more permissive. Both are emitted so a figure can name
    which one it used, and so the discrepancy with the R original is visible
    rather than discovered by a reviewer.

    or_floored marks the cells where the odds ratio was exactly 0 or NaN, i.e.
    where log_or is a floor/cap and not a measured effect size.
    """
    missing = [g for g in list(focal_genes) + list(comparator_genes)
               if g not in mat.columns]
    if missing:
        print("not in the matrix (0 mutated patients), skipped: " +
              ", ".join(sorted(set(missing))))

    frames = []
    for focal in focal_genes:
        if focal not in mat.columns:
            continue
        others = [g for g in comparator_genes if g != focal]
        res = pairwise_or_fisher(mat, focal, others, patient_subset=patient_subset)
        if len(res) == 0:
            continue
        res = res.rename(columns={"q": "q_within_focal"})
        frames.append(res)

    if not frames:
        return pd.DataFrame()

    out = pd.concat(frames, ignore_index=True)
    out["q_global"] = multipletests(out["p"].values, method="fdr_bh")[1]
    out["log_or"] = safe_log_or(out["odds_ratio"].values, eps=eps, cap=cap)
    ors = out["odds_ratio"].astype(float)
    out["or_floored"] = (ors == 0) | ors.isna() | np.isinf(ors)

    out = classify_significance(
        out.rename(columns={"q_within_focal": "q"}), q_thresh=q_thresh)
    out = out.rename(columns={"q": "q_within_focal"})
    out["q_thresh"] = q_thresh
    out["fdr_family"] = "within_focal"
    return out.reset_index(drop=True)


def comut_volcano_frame(res, correction="q_within_focal"):
    """Plot-ready frame for a co-occurrence volcano.

    This is the DATA only -- x, y, point size, category and whether the point
    is a floored odds ratio. Fonts, palette, panel size and PDF export come
    from the sci-plot skill. Keeping the mechanics here and the style there is
    deliberate: nothing in this skill declares rcParams.

    `correction` names which q column drives y and the category, and is
    returned in the frame so the figure can state it. CH.ipynb's volcano and
    its bubble plot disagreed about this -- the volcano plotted -log10(q) while
    the bubble was called with q_col='p_value' -- so the same data showed
    different points as significant in two figures.
    """
    if correction not in res.columns:
        raise SystemExit("no column " + repr(correction) + " in the results; "
                         "expected q_within_focal, q_global or p")
    out = res.copy()
    out["x_log_or"] = out["log_or"]
    out["y_neglog10"] = -np.log10(np.clip(out[correction].values, 1e-300, 1.0))
    out["size_n_double"] = out["sintersect"]
    out["correction_used"] = correction
    return out
