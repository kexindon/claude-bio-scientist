"""Guide-level LFC / empirical-FDR engine for pooled sensor screens.

Lifted verbatim from
  fsrlab-sensor-analysis/NGS-scripts/STEP9_LFC_FDR_calculation/LFC_FDR_calculation.ipynb
  cell  5  -> rpm
  cell  7  -> LFC_table_generator
  cell  9  -> attach_editing_info
  cell 11  -> bootstrap_null, plot_null_distributions
  cell 13  -> empirical_p_value, two_sided_FDR
  cell 15  -> run_screen_pipeline
  cell 19  -> _mle_from_compact  (here: mle_from_compact, BASE_crispr as an argument)
Do not edit the bodies here without editing there -- these two must not drift.
The published LFC-FDR/*.csv tables were produced by that notebook, so any change
to the math here silently de-calibrates every figure that reads those tables.

The one piece of NEW code is `empirical_p_matrix`, a np.searchsorted
reimplementation of the `empirical_p_value` double loop. It is verified to return
values *identical* to the loop (np.array_equal on both tails, including NaN, +-inf
and exact ties against the null) and turns ~2 min per tail per condition into
milliseconds. `empirical_p_value` is kept exported so the equivalence can be
re-checked at any time -- see `two_sided_FDR(..., legacy_loop=True)`.

Imports are deliberately minimal. The notebook's import cell also pulls
upsetplot / adjustText / matplotlib_venn / sklearn / seaborn, none of which this
pipeline uses; importing them here would make the FDR table unobtainable on a
clean env for no reason. matplotlib is imported lazily inside
`plot_null_distributions`, so the statistics run on a headless box with no
plotting stack. The notebook's `warnings.filterwarnings('ignore')` is also NOT
carried over -- the warnings it hides (log2 divide-by-zero, nanmean of an empty
slice, merge dtype changes) are exactly the ones that matter here.

Figure styling is NOT this module's business. `plot_null_distributions` draws the
QC panel you glance at once; anything going into a paper (volcano, waterfall,
LFC-vs-editing) is handed to the `sci-plot` skill.
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests

CONTROL_CLASSES = ['non-targeting control', 'safe-targeting control']

# Editing columns attach_editing_info pulls out of each MLE CSV.
MLE_COLS = ['Guide_ID', 'Reads_aligned_all_amplicons', 'corr_perc',
            'target_base_edit_perc', 'byproduct_INDEL_perc', 'byproduct_sub_perc']


# ---------------------------------------------------------------------------
# Step 1 -- depth normalization                                  (notebook c.5)
# ---------------------------------------------------------------------------

def rpm(df):
    df_rpm = df.copy()
    pseudocount = 1
    samples = list(df.columns[2:])
    for k in samples:
        col = np.asarray(df[k]) + pseudocount
        sum_col = np.sum(col)
        rpm = 1e6 * (col / sum_col)
        df_rpm[k] = rpm
    return df_rpm


# ---------------------------------------------------------------------------
# Step 2 -- LFC, replicate summaries, z-scores                   (notebook c.7)
# ---------------------------------------------------------------------------

def LFC_table_generator(RPM_df, samp_dict, t0_comparison=['input_rep1', 'input_rep2', 'input_rep3'],
                        ambiguity_threshold=2/3,
                        baseline='t0',
                        RPM_df_bc=None, samp_dict_bc=None):
    """
    Function for generating LFC values
    RPM_df = dataframe containing RPM values (generated with pseudocount of +1)
    samp_dict = dictionary of samples and corresponding replicates
    t0_comparison = list containing base mean values for LFC calculation (used when baseline='t0')
    ambiguity_threshold = minimum fraction of replicates agreeing on direction for concordant mean (default: 2/3)
    baseline = 't0' (default) or 'barcoding'.
        't0': single global baseline = median of t0_comparison columns in RPM_df.
        'barcoding': per-condition baseline = median of matching barcoding replicates in RPM_df_bc.
            Requires RPM_df_bc (barcoding RPM dataframe) and samp_dict_bc (dict mapping
            condition -> list of barcoding replicate column names).
    """

    # global baseline (used only in 't0' mode)
    if baseline == 't0':
        t0_median = np.median(RPM_df[t0_comparison], axis=1)
    elif baseline == 'barcoding':
        if RPM_df_bc is None or samp_dict_bc is None:
            raise ValueError("baseline='barcoding' requires RPM_df_bc and samp_dict_bc")
    else:
        raise ValueError("baseline must be 't0' or 'barcoding', got " + repr(baseline))

    df_holder = []
    for condition_of_interest in list(samp_dict.keys()):

        samp_of_interest = samp_dict[condition_of_interest]

        # determine per-condition baseline
        if baseline == 't0':
            cond_baseline = t0_median
        else:  # barcoding
            samp_of_interest_bc = samp_dict_bc[condition_of_interest]
            cond_baseline = np.median(RPM_df_bc[samp_of_interest_bc], axis=1)

        # calculate log2fc
        lfc_list = []
        for k in samp_of_interest:
            col = RPM_df[k]
            lfc = np.log2(col / cond_baseline)
            lfc_list.append(lfc)

        lfc_table = pd.DataFrame(dict(zip(samp_of_interest, lfc_list)))
        lfc_table['LFC_avg_' + condition_of_interest] = np.average(lfc_table[samp_of_interest], axis=1)
        lfc_table['LFC_median_' + condition_of_interest] = np.median(lfc_table[samp_of_interest], axis=1)

        # concordant mean: majority direction determined by sign vote across replicates
        rep_lfcs = lfc_table[samp_of_interest].values  # (n_guides, n_reps)
        n_reps = len(samp_of_interest)
        signs = np.sign(rep_lfcs)  # -1, 0, or 1

        n_positive = np.sum(signs > 0, axis=1)  # number of replicates with positive LFC
        n_negative = np.sum(signs < 0, axis=1)  # number of replicates with negative LFC
        majority_positive = n_positive >= n_negative  # True = majority direction is positive

        agrees = np.where(majority_positive[:, np.newaxis], signs >= 0, signs <= 0)

        frac_concordant = np.sum(agrees, axis=1) / n_reps
        is_directional = frac_concordant >= ambiguity_threshold  # True = concordant, False = ambiguous

        rep_lfcs_concordant = np.where(agrees, rep_lfcs, np.nan)
        concordant_mean = np.nanmean(rep_lfcs_concordant, axis=1)
        concordant_mean[~is_directional] = np.nan

        lfc_table['LFC_concordant_mean_' + condition_of_interest] = concordant_mean
        lfc_table['LFC_directional_' + condition_of_interest] = is_directional

        # z-scores for each summary metric: z = (S - mean(S)) / std(S)
        for summary in ['avg', 'median', 'concordant_mean']:
            col = lfc_table['LFC_' + summary + '_' + condition_of_interest]
            avg = np.average(col[~np.isnan(col)])
            std = np.std(col[~np.isnan(col)])
            lfc_table['z_score_' + summary + '_' + condition_of_interest] = (col - avg) / std

        lfc_table['Input_median'] = cond_baseline
        ids = RPM_df[['sgRNA', 'gene']]

        LFC_df = pd.concat((ids, lfc_table), axis=1)
        df_holder.append(LFC_df)

    # and combine different tissues into same table
    combined_LFC = pd.concat(df_holder, axis=1)
    # axis=1 concat repeats sgRNA/gene/Input_median once per condition -- this
    # dedup is load-bearing, not cosmetic. Dropping it leaves duplicate column
    # labels and every later df['sgRNA'] returns a DataFrame, not a Series.
    combined_LFC = combined_LFC.loc[:, ~combined_LFC.columns.duplicated()].copy()
    combined_LFC = combined_LFC.rename(columns={'sgRNA': 'gRNA_id', 'gene': 'Gene'})

    return combined_LFC


# ---------------------------------------------------------------------------
# Step 3 -- sensor editing measurements                          (notebook c.9)
# ---------------------------------------------------------------------------

def attach_editing_info(LFC_df, fp_editing, mle_names, mle_affix, valid_grnas, verbose=True):
    """
    Merge sensor editing measurements into an LFC table.

    LFC_df = LFC dataframe with 'gRNA_id' column
    fp_editing = path to folder containing MLE editing CSVs
    mle_names = list of condition names (e.g. ['input', 'd5', 'd15', 'bm', 'spleen', 'men'])
    mle_affix = file name suffix (e.g. 'ABE_MLE'); use '' for no suffix
    valid_grnas = array-like of valid guide IDs to keep (e.g. df_counts['sgRNA'])
    """
    LFC = LFC_df.copy()
    for condition in mle_names:
        fname = condition + '_' + mle_affix + '.csv' if mle_affix else condition + '.csv'
        editing = pd.read_csv(os.path.join(fp_editing, fname))
        editing = editing[MLE_COLS].rename(
            columns={'Guide_ID': 'gRNA_id',
                     'Reads_aligned_all_amplicons': 'sensor_reads_' + condition,
                     'corr_perc': 'corr_perc_' + condition,
                     'target_base_edit_perc': 'target_base_edit_perc_' + condition,
                     'byproduct_INDEL_perc': 'byproduct_INDEL_perc_' + condition,
                     'byproduct_sub_perc': 'byproduct_sub_perc_' + condition}
        )
        editing = editing[editing['gRNA_id'].isin(valid_grnas)]
        before = len(LFC)
        # NOTE: default inner join, no size assertion in the notebook. A guide
        # missing from one MLE CSV is dropped from the WHOLE table, and the
        # merges cascade, so one bad condition file shrinks everything. We do not
        # change the join (that would change published output) -- we report it.
        LFC = pd.merge(LFC, editing, on='gRNA_id')
        if verbose and len(LFC) != before:
            print('  [warn] MLE merge "' + condition + '" dropped '
                  + str(before - len(LFC)) + ' guides (' + str(before)
                  + ' -> ' + str(len(LFC)) + ')')
    return LFC


# ---------------------------------------------------------------------------
# Step 4 -- empirical null from control guides                   (notebook c.11)
# ---------------------------------------------------------------------------

def bootstrap_null(LFC_df, conditions, samp_dict, control_classes=CONTROL_CLASSES,
                   n_bootstrap=10000, seed=0):
    """
    Build empirical null distributions by bootstrap-resampling control guide LFCs.

    Returns (null_dict {condition: array of n_bootstrap LFCs}, controls DataFrame).
    seed=0 and n_bootstrap=10000 are what the published tables used -- changing
    either moves p-values near the 1/(n_bootstrap+1) floor and changes hit counts.
    """
    np.random.seed(seed)
    controls = LFC_df[LFC_df['classification'].isin(control_classes)].reset_index(drop=True)
    if len(controls) == 0:
        raise SystemExit(
            "No control guides found. 'classification' must contain the exact "
            "strings " + repr(control_classes) + "; present values were "
            + repr(sorted(LFC_df['classification'].dropna().unique())[:12]))
    null_dict = {}
    for condition in conditions:
        samp_of_interest = samp_dict[condition]
        nt_gRNA_values = controls[samp_of_interest].to_numpy().flatten()
        n_guides = len(nt_gRNA_values) // len(samp_of_interest)
        print(condition + ': ' + str(len(nt_gRNA_values)) + ' control guide LFC values ('
              + str(n_guides) + ' guides x ' + str(len(samp_of_interest)) + ' replicates)')
        if len(nt_gRNA_values) < 100:
            print('  [warn] fewer than 100 pooled control LFCs -- the bootstrap '
                  'cannot manufacture resolution the controls do not have. '
                  'Check the null QC panel before trusting the FDRs.')
        null_dict[condition] = np.random.choice(nt_gRNA_values, size=n_bootstrap, replace=True)
    return null_dict, controls


def plot_null_distributions(null_dict, controls, samp_dict, bins=np.linspace(-11, 8, 80),
                            name_dict=None, color='tab:purple', save_path=None, show=True):
    """QC panel: bootstrapped null overlaid with each replicate's own control LFCs.

    This is the glance-at-it-once diagnostic, deliberately plain. Paper figures
    (volcano, waterfall, LFC vs editing) go through the `sci-plot` skill.
    """
    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['pdf.fonttype'] = 42  # keep PDF text editable

    conditions = list(null_dict.keys())
    n = len(conditions)
    if name_dict is None:
        name_dict = dict((c, c) for c in conditions)

    fig, ax = plt.subplots(1, n, figsize=(3.5 * n, 3.5), sharex=True, sharey=False)
    if n == 1:
        ax = [ax]

    for idx, condition in enumerate(conditions):
        ax[idx].hist(null_dict[condition], density=True, bins=bins, alpha=.7,
                     color=color, label='Bootstrapped Null')
        for idx2, k in enumerate(samp_dict[condition]):
            label = 'Replicates' if idx2 == 0 else None
            ax[idx].hist(controls[k], density=True, bins=bins, alpha=.2, linewidth=.5,
                         color='black', histtype='step', label=label)
        ax[idx].set_title(name_dict[condition])
        ax[idx].spines[['top', 'right']].set_visible(False)
        ax[idx].legend(fontsize=8, loc='upper left')

    ax[0].set_ylabel('Density')
    ax[n // 2].set_xlabel('LFC')
    fig.tight_layout()

    if save_path is not None:
        d = os.path.dirname(save_path)
        if d:
            os.makedirs(d, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print('Figure saved to ' + save_path)
    if show:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Step 5-6 -- empirical p, combine across replicates, BH         (notebook c.13)
# ---------------------------------------------------------------------------

def empirical_p_value(observed, null_distribution):
    """Scalar reference implementation, kept for equivalence checking."""
    S_high = sum(null_distribution >= observed)
    S_low = sum(null_distribution <= observed)
    N = len(null_distribution)

    # Compute two-sided p-value
    p_high = (S_high + 1) / (N + 1)
    p_low = (S_low + 1) / (N + 1)
    return p_high, p_low


def empirical_p_matrix(observed, null_distribution):
    """Vectorized `empirical_p_value` over a 2-D (guides x replicates) array.

    NEW code, not in the notebook. Verified np.array_equal to the scalar loop on
    both tails for random data, exact ties against null values, +inf, -inf and
    NaN -- see the NaN note below, which is the one case naive searchsorted gets
    wrong.

      p_high = (1 + #(null >= obs)) / (1 + N),  #(null >= obs) = N - searchsorted(left)
      p_low  = (1 + #(null <= obs)) / (1 + N),  #(null <= obs) = searchsorted(right)
    """
    obs = np.asarray(observed, dtype=float)
    s = np.sort(np.asarray(null_distribution, dtype=float))
    N = len(s)

    # np.searchsorted sorts NaN to the end, which would give p_low = 1.0; the
    # scalar loop's comparisons are all False for NaN, so it returns the 1/(N+1)
    # floor on BOTH tails. Substitute a dummy, then restore the floor.
    is_nan = np.isnan(obs)
    safe = np.where(is_nan, 0.0, obs)

    p_high = (N - np.searchsorted(s, safe, side='left') + 1) / (N + 1)
    p_low = (np.searchsorted(s, safe, side='right') + 1) / (N + 1)

    floor = 1.0 / (N + 1)
    p_high = np.where(is_nan, floor, p_high)
    p_low = np.where(is_nan, floor, p_low)
    return p_high, p_low


def _combine(p_matrix, method):
    """Combine each row of p-values. combine_pvalues returns a tuple in old
    scipy and a SignificanceResult in new scipy; [1] is the p-value in both."""
    return np.array([combine_pvalues(p_vals, method=method)[1] for p_vals in p_matrix])


def two_sided_FDR(condition_of_interest, samp_of_interest, LFC_df, null_dict,
                  legacy_loop=False):
    '''
    condition_of_interest = name of condition (e.g. "spleen")
    samp_of_interest = list of all replicates ["spleen1", "spleen2", "spleen3"]
    LFC_df = dataframe with LFC values
    null_dict = dictionary containing bootstrapped null distributions for each tissue & editor combination

    Mutates LFC_df IN PLACE and returns the same object. run_screen_pipeline
    relies on that: its df_holder ends up holding N aliases of one DataFrame, and
    the pipeline is only correct because the concat is followed by the
    .duplicated() dedup. Do not "fix" the aliasing without rereading the caller.

    legacy_loop=True uses the original pure-Python double loop (minutes per tail
    on a 13k-guide screen) instead of the vectorized form, for equivalence checks.
    '''

    # Load the null distribution
    nt_gRNA_values = null_dict[condition_of_interest]

    # Targeting gRNAs (including NT for completeness)
    gRNA_values = np.asarray(LFC_df[samp_of_interest])

    # Compute p-values for each replicate separately (two-tailed)
    if legacy_loop:
        p_values_high = np.array([
            [empirical_p_value(rep, nt_gRNA_values)[0] for rep in gRNA]
            for gRNA in gRNA_values
        ])
        p_values_low = np.array([
            [empirical_p_value(rep, nt_gRNA_values)[1] for rep in gRNA]
            for gRNA in gRNA_values
        ])
    else:
        p_values_high, p_values_low = empirical_p_matrix(gRNA_values, nt_gRNA_values)

    # Combine p-values across replicates using Fisher's and Stouffer's methods
    combined_p_values_high = _combine(p_values_high, 'fisher')
    combined_p_values_low = _combine(p_values_low, 'fisher')
    combined_p_values_high_stouffers = _combine(p_values_high, 'stouffer')
    combined_p_values_low_stouffers = _combine(p_values_low, 'stouffer')

    # Adjust for multiple testing (Benjamini-Hochberg FDR). NOTE: this corrects
    # across EVERY row of the table, control guides included -- published
    # behavior, and the reason a control guide can come out with FDR < 0.1.
    adjusted_p_values_high = multipletests(combined_p_values_high, method='fdr_bh')[1]
    adjusted_p_values_low = multipletests(combined_p_values_low, method='fdr_bh')[1]
    adjusted_p_values_high_stouffers = multipletests(combined_p_values_high_stouffers, method='fdr_bh')[1]
    adjusted_p_values_low_stouffers = multipletests(combined_p_values_low_stouffers, method='fdr_bh')[1]

    # Store results in the dataframe
    c = condition_of_interest
    LFC_df['p_high_unadjusted_' + c + '_fishers'] = combined_p_values_high
    LFC_df['p_low_unadjusted_' + c + '_fishers'] = combined_p_values_low
    LFC_df['FDR_high_' + c + '_fishers'] = adjusted_p_values_high
    LFC_df['FDR_low_' + c + '_fishers'] = adjusted_p_values_low

    LFC_df['p_high_unadjusted_' + c + '_stouffers'] = combined_p_values_high_stouffers
    LFC_df['p_low_unadjusted_' + c + '_stouffers'] = combined_p_values_low_stouffers
    LFC_df['FDR_high_' + c + '_stouffers'] = adjusted_p_values_high_stouffers
    LFC_df['FDR_low_' + c + '_stouffers'] = adjusted_p_values_low_stouffers

    # Combine low and high FDRs for ease of manipulation and plotting.
    # The notebook does this with an iterrows() loop over min(FDR_low, FDR_high);
    # these two arrays are the same two columns it reads, elementwise.
    min_FDR = []
    for _, val in LFC_df.iterrows():
        FDR_low = val['FDR_low_' + c + '_fishers']
        FDR_high = val['FDR_high_' + c + '_fishers']
        min_FDR.append(min(FDR_low, FDR_high))

    LFC_df['FDR_' + c + '_fishers'] = min_FDR

    min_FDR = []
    for _, val in LFC_df.iterrows():
        FDR_low = val['FDR_low_' + c + '_stouffers']
        FDR_high = val['FDR_high_' + c + '_stouffers']
        min_FDR.append(min(FDR_low, FDR_high))

    LFC_df['FDR_' + c + '_stouffers'] = min_FDR

    return LFC_df


# ---------------------------------------------------------------------------
# MAIN                                                          (notebook c.15)
# ---------------------------------------------------------------------------

def run_screen_pipeline(screen_name, df_counts, library, samps_to_merge, new_names,
                        conditions_null,
                        fp_editing=None, mle_names=None, mle_affix='',
                        to_comparison_samps=None,
                        df_counts_bc=None, samps_to_merge_bc=None, new_names_bc=None,
                        control_classes=CONTROL_CLASSES,
                        ambiguity_threshold=2/3, n_bootstrap=10000, seed=0,
                        null_color='tab:purple',
                        save_path=None, show=True, legacy_loop=False):
    """
    Run the full LFC -> null -> FDR pipeline for a single screen.

    Returns (combined_LFC_FDR, null_dict, controls).

    Baseline mode is auto-detected: df_counts_bc given -> 'barcoding',
    otherwise 't0' from to_comparison_samps. Passing both is rejected here
    (the notebook silently ignored to_comparison_samps in that case, which made
    a barcoding run read like a T0 run in the config).

    Set fp_editing=None to skip the sensor-editing merge entirely (screens with
    no sensor data).
    """
    print('===== ' + screen_name + ' =====')

    # Determine baseline mode
    if df_counts_bc is not None:
        if samps_to_merge_bc is None or new_names_bc is None:
            raise ValueError("df_counts_bc provided but samps_to_merge_bc or new_names_bc is missing")
        if to_comparison_samps is not None:
            raise ValueError(
                "both df_counts_bc and to_comparison_samps given -- the notebook "
                "silently ignored the T0 baseline here. Pick one.")
        baseline = 'barcoding'
    else:
        if to_comparison_samps is None:
            raise ValueError("to_comparison_samps is required when not using barcoding baseline")
        baseline = 't0'

    if len(samps_to_merge) != len(new_names):
        raise ValueError('samps_to_merge has ' + str(len(samps_to_merge)) + ' groups but '
                         + 'new_names has ' + str(len(new_names)) + ' labels')
    missing = [c for c in conditions_null if c not in new_names]
    if missing:
        raise ValueError('conditions_null entries not in new_names: ' + repr(missing))

    # Step 1: RPM normalization
    df_rpm = rpm(df_counts)
    samp_dict = dict(zip(new_names, samps_to_merge))

    # Step 2: LFC table generation
    if baseline == 'barcoding':
        df_rpm_bc = rpm(df_counts_bc)
        samp_dict_bc = dict(zip(new_names_bc, samps_to_merge_bc))
        LFC = LFC_table_generator(df_rpm, samp_dict,
                                  ambiguity_threshold=ambiguity_threshold,
                                  baseline='barcoding',
                                  RPM_df_bc=df_rpm_bc, samp_dict_bc=samp_dict_bc)
    else:
        LFC = LFC_table_generator(df_rpm, samp_dict, t0_comparison=to_comparison_samps,
                                  ambiguity_threshold=ambiguity_threshold)

    # Step 3: Attach library and editing information
    n_before = len(LFC)
    LFC = pd.merge(LFC, library, on='gRNA_id')
    if len(LFC) != n_before:
        print('  [warn] library merge dropped ' + str(n_before - len(LFC))
              + ' guides (' + str(n_before) + ' -> ' + str(len(LFC))
              + ') -- guides absent from the library annotation')
    if fp_editing is not None:
        print('  MLE folder: ' + fp_editing + ' (conditions: ' + ','.join(mle_names or []) + ')')
        LFC = attach_editing_info(LFC, fp_editing=fp_editing,
                                  mle_names=mle_names, mle_affix=mle_affix,
                                  valid_grnas=df_counts['sgRNA'])
    else:
        print('  MLE folder: none (sensor editing merge skipped)')

    # Step 4: Bootstrap null distribution
    null_dict, controls = bootstrap_null(LFC, conditions_null, samp_dict,
                                         control_classes=control_classes,
                                         n_bootstrap=n_bootstrap, seed=seed)

    # Visualize null distributions
    plot_null_distributions(null_dict, controls, samp_dict, color=null_color,
                            save_path=save_path, show=show)

    # Steps 5-6: Empirical p-values, combined p-values, and FDR
    # two_sided_FDR mutates LFC in place, so df_holder is N aliases of one frame;
    # the dedup below is what makes that correct.
    df_holder = []
    for condition in conditions_null:
        samp_of_interest = samp_dict[condition]
        LFC_df = two_sided_FDR(condition, samp_of_interest, LFC, null_dict,
                               legacy_loop=legacy_loop)
        df_holder.append(LFC_df)

    combined_LFC_FDR = pd.concat(df_holder, axis=1)
    combined_LFC_FDR = combined_LFC_FDR.loc[:, ~combined_LFC_FDR.columns.duplicated()].copy()

    print_hit_summary(combined_LFC_FDR, conditions_null)
    return combined_LFC_FDR, null_dict, controls


def print_hit_summary(combined_LFC_FDR, conditions_null, fdr_cut=0.1, lfc_cut=1):
    """Per-condition hit counts, lifted from the tail of run_screen_pipeline."""
    n_total = len(combined_LFC_FDR)
    for condition in conditions_null:
        fisher_mask = combined_LFC_FDR['FDR_' + condition + '_fishers'] < fdr_cut
        stouffer_mask = combined_LFC_FDR['FDR_' + condition + '_stouffers'] < fdr_cut
        overlap_mask = fisher_mask & stouffer_mask

        print('  ' + condition + ' (n=' + str(n_total) + '): '
              + str(fisher_mask.sum()) + ' Fisher hits, '
              + str(stouffer_mask.sum()) + ' Stouffer hits, '
              + str(overlap_mask.sum()) + ' overlap (FDR<' + str(fdr_cut) + ')')

        for method, mask in [('Fisher', fisher_mask), ('Stouffer', stouffer_mask)]:
            hits = combined_LFC_FDR[mask]
            for summary in ['avg', 'median', 'concordant_mean']:
                lfc_col = 'LFC_' + summary + '_' + condition
                if lfc_col not in hits.columns:
                    continue
                n_dep = (hits[lfc_col] < -lfc_cut).sum()
                n_enr = (hits[lfc_col] > lfc_cut).sum()
                print('    [' + method + '] |' + lfc_col + '|>' + str(lfc_cut) + ': '
                      + str(n_dep) + ' depleted, ' + str(n_enr) + ' enriched')

        if 'classification' in combined_LFC_FDR.columns:
            ctrl_hits = combined_LFC_FDR[fisher_mask & combined_LFC_FDR[
                'classification'].isin(CONTROL_CLASSES)]
            if len(ctrl_hits):
                print('    [note] ' + str(len(ctrl_hits)) + ' of the Fisher hits are CONTROL '
                      'guides -- BH correction includes controls. Exclude them before '
                      'reporting biology.')


# ---------------------------------------------------------------------------
# MLE rebuild from CRISPResso compact files                     (notebook c.19)
# ---------------------------------------------------------------------------
# Kept here, not in a separate skill, only so that the replicate-exclusion list
# used for the counts and the one used for the editing table come from ONE config.
# Generating the compact files themselves is `ngs_screen` territory.

def mle_from_compact(base_crispr, screen_dir, groups, condition_names, out_dir):
    """
    base_crispr    : root of the crispresso output tree (e.g. 'crispresso/focused/')
    groups         : list of lists of per-sample stems (compact file name without suffix)
    condition_names: list of condition labels, one per group
    out_dir        : destination folder for MLE CSVs
    """
    cols_to_add = [
        'Reads_in_input', 'Reads_lowqual', 'Reads_aligned_all_amplicons',
        'Reads_aligned_WT', 'Reads_aligned_HDR', 'Reads_unaligned', 'WT',
        'correct_edit', 'target_base_edit', 'byproduct_all',
        'byproduct_INDEL', 'byproduct_sub', 'byproduct_ambiguous',
    ]
    os.makedirs(out_dir, exist_ok=True)
    for group, cond in zip(groups, condition_names):
        dfs = [
            pd.read_csv(os.path.join(base_crispr, screen_dir, s + '_compact_unfiltered.csv'))
            for s in group
        ]
        merged = dfs[0][cols_to_add].copy()
        for df in dfs[1:]:
            merged = merged + df[cols_to_add].values
        merged.insert(0, 'Guide_ID', dfs[0]['Guide_ID'])
        r = merged['Reads_aligned_all_amplicons']
        merged['corr_perc'] = 100 * merged['correct_edit'] / r
        merged['target_base_edit_perc'] = 100 * merged['target_base_edit'] / r
        merged['WT_perc'] = 100 * merged['WT'] / r
        merged['byproduct_all_perc'] = 100 * merged['byproduct_all'] / r
        merged['byproduct_INDEL_perc'] = 100 * merged['byproduct_INDEL'] / r
        merged['byproduct_sub_perc'] = 100 * merged['byproduct_sub'] / r
        merged['byproduct_ambiguous_perc'] = 100 * merged['byproduct_ambiguous'] / r
        merged = merged.fillna(0)
        merged.to_csv(os.path.join(out_dir, cond + '.csv'), index=False)
    print('Saved ' + str(len(groups)) + ' MLE files -> ' + out_dir)


# ---------------------------------------------------------------------------
# shared library annotation (optional, via $LAB_DATA_ROOT)
# ---------------------------------------------------------------------------

LIB_COLS = ['gRNA_id', 'gene_name_m_corrected', 'gene_name_h',
            'HGVSp_m', 'HGVSp_h', 'classification', 'classification_m']


def load_library(path, columns=LIB_COLS):
    """Read a library annotation table.

    `path` may be an explicit file path, or `lab:<name>` to resolve
    $LAB_DATA_ROOT/libraries/<name>.csv.
    """
    if path.startswith('lab:'):
        root = os.environ.get('LAB_DATA_ROOT')
        if not root:
            raise SystemExit('LAB_DATA_ROOT is not set -- see docs/data-layout.md')
        path = os.path.join(root, 'libraries', path[4:] + '.csv')
    if not os.path.exists(path):
        raise SystemExit('library annotation not found: ' + path)
    lib = pd.read_csv(path)
    keep = [c for c in columns if c in lib.columns]
    for required in ['gRNA_id', 'classification']:
        if required not in keep:
            raise SystemExit('library ' + path + " is missing required column '"
                             + required + "'")
    return lib[keep]
