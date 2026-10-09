# Output column dictionary and the statistics behind each column

Source of the math: `fsrlab-sensor-analysis/docs/source/analysis.md`, section
"STEP 9 — Empirical LFC-FDR calculation". Source of the code:
`NGS-scripts/STEP9_LFC_FDR_calculation/LFC_FDR_calculation.ipynb`.

`<c>` below is a condition label from `new_names` / the `samps` keys —
`bm`, `spleen`, `men`, `d5`, `d15`, `lib`, `input`.

## The six steps

**1. RPM, pseudocount 1.** Per sample column *j*:

```
RPM_ij = 1e6 * (count_ij + 1) / sum_k (count_kj + 1)
```

The +1 is added before the depth sum, so the denominator is
`total_reads + n_guides`, not `total_reads`. It is what keeps zero-count guides
out of `log2(0)`; drop it and `-inf`/NaN LFCs propagate into `combine_pvalues`
and then into `multipletests`.

**2. LFC per replicate, against a baseline.**

```
LFC_ij = log2( RPM_ij,sample / median RPM_i,baseline )
```

Two baseline modes, selected by what the config supplies:

| mode | baseline | config |
|---|---|---|
| `t0` | median RPM across the T0/input columns, one global vector | `baseline_samps: [input]` |
| `barcoding` | median RPM across the *matching condition's* barcoding replicates | `counts_bc` + `samps_bc` |

In barcoding mode the baseline is taken **positionally** (`np.median(..., axis=1)`
on the barcoding RPM frame, then divided row-wise) — there is no merge on
`sgRNA`. The two count tables must therefore be the same guides in the same
order. `check_inputs.py` asserts this.

**3. Replicate summaries.** Per condition:

| column | meaning |
|---|---|
| `<replicate>` | the raw per-replicate LFC, under the replicate's own column name |
| `LFC_avg_<c>` | `np.average` across replicates |
| `LFC_median_<c>` | `np.median` across replicates — **the house default for reporting** |
| `LFC_concordant_mean_<c>` | mean over only the replicates agreeing with the majority sign; NaN when under 2/3 agree |
| `LFC_directional_<c>` | bool: did ≥ 2/3 of replicates agree on direction |
| `Input_median` | the baseline vector used (one column, last condition's wins) |

Concordant mean, exactly: sign-vote across replicates; ties break **positive**
(`n_positive >= n_negative`); replicates whose sign matches the majority
(including zeros on the majority side) are averaged; if the agreeing fraction is
below `ambiguity_threshold = 2/3` the value is set to NaN.

**4. Z-scores.** For each of the three summaries, across all guides in that
condition, computed on non-NaN values only:

```
z = (S - mean(S[~isnan])) / std(S[~isnan])
```

`np.std`, so population SD (ddof=0). Columns: `z_score_avg_<c>`,
`z_score_median_<c>`, `z_score_concordant_mean_<c>`.

**5. Empirical null, and the two tail p-values.** Pool the LFCs of every
`safe-targeting control` + `non-targeting control` guide across that condition's
replicates (a flat array of `n_control_guides × n_replicates` values), then draw
`n_bootstrap = 10000` values from it with replacement at `seed = 0`. Per guide ×
replicate:

```
p_high = (1 + #(null >= LFC_obs)) / (1 + N_null)
p_low  = (1 + #(null <= LFC_obs)) / (1 + N_null)
```

so the smallest attainable p is `1/10001 ≈ 1.0e-4` regardless of how many real
control observations stand behind the null. Tails are kept separate so direction
is assessed independently.

**6. Combine across replicates, then BH.** `scipy.stats.combine_pvalues` with
both `method='fisher'` and `method='stouffer'`, then
`statsmodels.stats.multitest.multipletests(method='fdr_bh')` on each, then

```
FDR_<c> = min(FDR_high_<c>, FDR_low_<c>)
```

Fisher is sensitive to one extreme replicate; Stouffer needs more consistent
evidence across replicates. Both assume the tests are independent, which
biological replicates drawn from one shared initial guide pool do not satisfy —
so FDR < 0.1 here is a ranking convention, not a calibrated error rate.

## Significance columns, per condition in `conditions_null`

| column | |
|---|---|
| `p_high_unadjusted_<c>_fishers` | Fisher-combined enrichment p |
| `p_low_unadjusted_<c>_fishers` | Fisher-combined depletion p |
| `FDR_high_<c>_fishers` | BH over the enrichment p's |
| `FDR_low_<c>_fishers` | BH over the depletion p's |
| `FDR_<c>_fishers` | `min` of the two — **what to filter on** |
| `*_stouffers` | the same five with Stouffer's method |

BH runs over **every row of the table, control guides included**. That is
published behavior, and it is why a control guide can land at FDR < 0.1.

## Sensor editing columns, per condition in `mle_names`

`sensor_reads_<c>` (= `Reads_aligned_all_amplicons`), `corr_perc_<c>`,
`target_base_edit_perc_<c>`, `byproduct_INDEL_perc_<c>`, `byproduct_sub_perc_<c>`.
Merged in from STEP 7's per-condition MLE CSVs on `Guide_ID` → `gRNA_id`. These
make LFC-vs-editing calibration possible from the single output table.

## Library annotation columns

Carried over from the library CSV: `gRNA_id`, `gene_name_m_corrected`,
`gene_name_h`, `HGVSp_m`, `HGVSp_h`, `classification`, `classification_m`, plus
`Gene` (renamed from the count table's `gene`).

`classification` is matched on **exact strings**. `targeting guide`,
`safe-targeting control`, `non-targeting control`, and in the focused library
`essential truncation guide` — which is a positive control for depletion and is
deliberately **not** part of the null.

## Why this instead of MAGeCK

From `analysis.md`: MAGeCK was run first and consistently over-called, flagging
large fractions of guides as enriched or depleted including in noisy tissues
(CBE meninges). MAGeCK is built for CRISPR-nuclease screens where several guides
per gene collapse into one gene-level score; this analysis is about individual
guide behavior and wants replicate-level variability preserved rather than
collapsed early. The pipeline here is the conservative non-parametric
alternative: the null comes from the library's own control guides, per tissue and
per screen, with no parametric model imposed.

So "why does MAGeCK call so many hits" is answered by *this skill*, but running
MAGeCK itself is STEP 8 and belongs to `ngs-screen`.

## Reference shape

`ABE_OG_LFC_FDR.csv`: 1070 guides × 152 columns, from 5 conditions in
`conditions_null` plus `lib`, with 6 MLE conditions merged. 1034 targeting
guides, 26 safe-targeting and 10 non-targeting controls.
