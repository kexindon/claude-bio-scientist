---
name: screen-stats
description: >-
  Call hits in a pooled base- or prime-editing sensor screen: from a guide x
  sample count matrix to per-guide log2 fold-change and an empirical FDR built
  from the library's own control guides. Use whenever asked to "call hits",
  "compute LFC and FDR", "which guides are enriched or depleted in bone marrow
  / spleen / meninges", "is this guide significant", "why are all my p-values
  1e-4", "why not just use MAGeCK", "re-run without spleen2", "make the
  LFC-FDR table". Starts at the count matrix, not FASTQs (ngs-screen); plots
  the table through sci-plot, not here.
---

# screen-stats — counts → LFC → empirical FDR → hit list

Three commands. Preflight, run, and (when a replicate is in question) diff.

```bash
S=~/.claude/skills/screen-stats/scripts

python $S/check_inputs.py --config screens.yaml     # before a run, not after
python $S/run_screen.py   --config screens.yaml --out LFC-FDR
python $S/compare_hits.py LFC-FDR/A_LFC_FDR.csv LFC-FDR/B_LFC_FDR.csv
```

Per screen you get:

```
LFC-FDR/<screen>_LFC_FDR.csv              one wide row per guide (152 cols for ABE_OG)
LFC-FDR/<screen>_null_dist.pdf            null QC panel — look at it
LFC-FDR/<screen>_LFC_FDR.provenance.json  seed, n_bootstrap, replicate lists, MLE
                                          folder, control counts, counts file hash
```

## Where the boundaries are

Starts at the count matrix. FASTQs, CRISPResso, sensor-amplicon extraction,
building the count matrix and running MAGeCK are **`ngs-screen`**. Designing the
library and its control guides is **`library-design`**. Every figure built from
the output table — volcano, waterfall, LFC-vs-editing, tissue heatmap — goes
through **`sci-plot`**; nothing in this skill defines a figure style. The one
plot it does draw, the null-distribution QC panel, is deliberately plain and is
for glancing at once.

Why MAGeCK was dropped *is* this skill's question; running MAGeCK is not.

## Reporting default

**`LFC_median_<cond>` + `FDR_<cond>_fishers` < 0.1, |LFC| > 1, controls
excluded.** The pipeline also computes `LFC_avg`, `LFC_concordant_mean`, three
z-scores and the whole Stouffer family — all of which belong in supplementary
data. The published figures are median + Fisher. Plot anything else and the
figure will not match the paper.

## Inputs

**1. Count matrix** — tab-separated, the first two columns *exactly* `sgRNA` and
`gene`, then one integer column per sequenced sample. Column order beyond those
two does not matter, but every name must match the config exactly.

**2. Library annotation** — CSV with `gRNA_id` and `classification`, matched on
the **exact strings** `targeting guide`, `safe-targeting control`,
`non-targeting control`. `gene_name_m_corrected`, `gene_name_h`, `HGVSp_m`,
`HGVSp_h`, `classification_m` are carried through when present. Either a path,
or `lab:<name>` to resolve `$LAB_DATA_ROOT/libraries/<name>.csv`.

**3. MLE editing CSVs** (optional) — one per condition, keyed on `Guide_ID`, from
STEP 7. Omit `mle_dir` and the merge is skipped entirely.

**4. The config** — one YAML (or JSON, same keys) that owns the replicate
exclusions, so the counts and the editing table can never disagree about which
mice are in:

```yaml
screens:
  ABE_FC_EPO_T0:
    counts: counts/focused/FC_ABE_EPO_COUNTS.txt
    library: MBESv2_focused.csv
    samps:
      bm:     [bm1, bm3, bm4, bm5]            # bm2 out: 90/10 skew
      spleen: [spleen1, spleen3, spleen4, spleen5]   # spleen2 out: same
      d15:    [d15-rep1, d15-rep2, d15-rep3, d15-rep4, d15-rep5]
      d5:     [d5]
      input:  [input]
      lib:    [lib]
    baseline_samps: [input]                   # T0 mode
    conditions_null: [d5, d15, spleen, bm]    # ORDER MATTERS — see Gotchas
    mle_dir: MLE/focused/ABE_EPO
    mle_names: [input, d5, d15, bm, spleen]
    mle_affix: ''
```

For the barcoding baseline, give `counts_bc` + `samps_bc` **instead of**
`baseline_samps`, never both. `python run_screen.py --print-config-doc` prints
the full format. A single screen can also be driven entirely from flags
(`--counts --library --samps 'bm=bm1,bm3' --baseline-samps input
--conditions-null bm`).

*Which* replicate to exclude is decided upstream, by the skew/correlation QC in
the project's `global_analysis.ipynb`, not here. This skill records the decision
and lets you measure what it cost.

## The method, in one pass

RPM with a **+1 pseudocount** → per-replicate **log2 FC** against a T0-input or
per-condition barcoding median → replicate **mean / median / concordant-mean**
(2/3 directional agreement) and z-scores → a **bootstrap null** of 10,000 values
resampled from the pooled `safe-targeting` + `non-targeting` control LFCs, built
**separately per tissue and per screen** → **two-sided empirical p** per guide
per replicate, enrichment and depletion tails kept apart → combined across
replicates with **both Fisher and Stouffer** → **BH-adjusted** → reported as
`FDR = min(FDR_high, FDR_low)`.

It is the conservative non-parametric alternative to MAGeCK, which over-called
badly on these screens because it is built to collapse several guides per gene
into one gene-level score, while this analysis is about individual guides.

Full formulas, every output column, and the MAGeCK argument a PI will ask about:
`references/columns.md`.

**The honest caveat, from the lab's own docs:** Fisher and Stouffer both assume
the combined tests are independent. Biological replicates here share one initial
guide pool and have correlated variability, so that assumption is only
approximate and FDR < 0.1 is a ranking convention, not a calibrated error rate.
Fisher is the reported default because it is conventional, not because it is
better behaved — it can call a hit off one extreme replicate, which is exactly
why the run prints the Fisher/Stouffer overlap. A large Fisher-only count is a
signal to go and look at the replicates.

## Gotchas

**The order of `conditions_null` changes the answer.** `np.random.seed(seed)` is
called once, then each condition draws from the shared RNG stream, so a
condition's null depends on its *position* in the list. Verified on real data:
`ABE_OG` bone marrow gives **168** Fisher hits at `[d5, d15, spleen, bm, men]`
and **177** at `[bm, d5, d15, spleen, men]` — same counts, same seed, same LFCs.
This is a bug in the source, not a feature, and it is not fixed here because
fixing it would stop reproducing the published tables. So: **keep
`conditions_null` in the order the published config used**, and treat the
provenance JSON's `conditions_null` as part of the result. If you only want one
tissue re-tested, re-run the whole list and read that tissue's columns, rather
than running the tissue alone.

**All p-values bottom out at 1.0e-4.** The floor is `1/(n_bootstrap+1)`, and it
is there whether 36 control guides or 3,600 stand behind the null. "My p-values
are all 0.0001" means the guide is past the end of the null, not that it is
10,000-fold significant. `ABE_OG` has **36 control guides** (26 safe-targeting +
10 non-targeting), so the 10,000-value bootstrap for bone marrow is resampling
180 underlying numbers. Bootstrapping cannot manufacture information. Look at the
null QC panel.

**Control guides are in the multiple-testing correction, and can be called
hits.** BH runs over every row of the table. In the published
`ABE_OG_LFC_FDR.csv`, 2 of the 168 bone-marrow Fisher hits are control guides.
That is intended behavior — do not "fix" it — but filter
`~classification.isin(['safe-targeting control','non-targeting control'])`
before you report a hit list, or you will report a control guide as biology.
`compare_hits.py` excludes controls by default; the run summary prints how many
of the hits are controls.

**Every merge is a silent inner join.** The library merge and the six cascading
MLE merges have no size assertion in the source. One guide missing from one MLE
CSV disappears from the *whole* table, and because the editing merges chain, one
bad condition file shrinks everything. On the `ABE_OG` fixture nothing drops
(1070 in, 1070 out), which is precisely why this goes unnoticed. The scripts here
print before/after counts and `check_inputs.py` predicts the loss per file.

**A stray index column corrupts every RPM value.** `rpm()` treats
`df.columns[2:]` as sample columns, full stop. A `to_csv`/`read_csv` round trip
that keeps the index turns `sgRNA` into an annotation column and `gene` into a
"sample". Nothing errors; every number is wrong. Preflighted.

**The barcoding baseline is positional, not merged.** In barcoding mode the
baseline is `np.median` over the barcoding table's rows and is then divided
row-wise — there is no join on `sgRNA`. Two tables with the same guides in a
different order will silently pair the wrong guides. Preflighted.

**`LFC_concordant_mean_*` is NaN for ambiguous guides, by design.** A comparison
(`> 1`) treats NaN as False, which is fine for counting, but a `.mean()` or a
volcano over that column drops ambiguous guides without saying so. The sign vote
also breaks ties *positive* (`n_positive >= n_negative`), so a 2-up/2-down guide
is called "positive" with 0.5 agreement, fails the 2/3 threshold, and becomes
NaN. Correct, but not symmetric.

**`seed=0` and `n_bootstrap=10000` are what the published tables used.** Changing
either moves p-values near the floor and changes hit counts. They are recorded in
the provenance sidecar; if "the numbers changed since last time", diff that file
first.

**Do not put single-replicate conditions in `conditions_null`.** `lib`,
`plasmidlib`, `d5` when it is one column: Fisher and Stouffer pass a lone p-value
through unchanged, so the "combined" p carries no replicate-consistency evidence
at all. The existing configs correctly leave them out; `check_inputs.py` warns.

**`*_LFC_FDR_FSR.csv` is not reproducible and must not be mixed in.** The
`LFC-FDR/` folder of the B-ALL project holds six `*_LFC_FDR_FSR.csv` tables read
by ten figure notebooks, and **nothing in any repo generates them** — there is no
`to_csv(..._FSR)` anywhere. They are narrower (101 vs 132 columns) and they *add
back* the replicates the main tables exclude: `ABE_FC_EPO_T0_LFC_FDR_FSR.csv`
carries `bm2` and `spleen2`, the two abnormally skewed samples. They disagree
substantially — for `ABE_FC_EPO_T0` bone marrow, 30 enriched hits in the main
table versus 10 in the FSR table, only 9 shared. Treat them as an unreproducible
legacy artifact: never put the two families in one figure, and if a figure
disagrees with a freshly computed table, check which family it loaded. Run
`compare_hits.py main.csv main_FSR.csv` to see the damage for your own screen.

**Importing the notebook's import cell would break a clean environment.** It
pulls `upsetplot`, `adjustText`, `matplotlib_venn` and `sklearn`, none of which
the pipeline uses. `screenstats.py` imports only `os`, numpy, pandas, scipy and
statsmodels at module level, and matplotlib lazily inside the one plotting
function — so the FDR table is obtainable on a headless box with no plotting
stack configured. Seaborn is not needed either. The notebook's blanket
`warnings.filterwarnings('ignore')` is also not carried over — the warnings it
hides (log2 divide-by-zero, nanmean of an empty slice, merge dtype changes) are
the ones that matter here.

## Environment

Plain scientific Python. **No `pegg`, no scikit-learn, no pickled models** — this
skill is *not* under the Python 3.9/3.10 + `scikit-learn==1.1.1` + `numpy<2` pin
that `library-design` and `h2m-variant` live with. Run it in whatever modern env
you like.

Verified end to end on Python 3.9.6 with numpy 1.26.4, pandas 2.2.2, scipy
1.13.0, statsmodels 0.14.4. Floors, not pins: `pandas>=1.5`, `scipy>=1.9`
(`combine_pvalues` changed its return type to a `SignificanceResult`; the code
indexes `[1]`, which works either way), `statsmodels>=0.13`. numpy 2.x is fine.
`PyYAML` only if the config is YAML rather than JSON. No `$LAB_DATA_ROOT` is
needed unless the library is given as `lab:<name>`.

Laptop-scale, ~10 s for a 1070-guide × 5-condition screen.

## Verifying

`scripts/screenstats.py` is lifted verbatim from
`fsrlab-sensor-analysis/NGS-scripts/STEP9_LFC_FDR_calculation/LFC_FDR_calculation.ipynb`
(cells 5, 7, 9, 11, 13, 15, 19), so the two cannot drift. Regenerating `ABE_OG`
through `run_screen.py` reproduces the published `LFC-FDR/ABE_OG_LFC_FDR.csv`
**exactly**: 1070 × 152, identical column order, all 144 numeric columns equal at
`atol=0`.

The one piece of new code is `empirical_p_matrix`, an `np.searchsorted`
replacement for the original pure-Python double loop. On the same screen the loop
takes 82 s and the vectorized form 10 s, with **bit-identical** output — verified
both as a whole table and elementwise against ties, `+inf`, `-inf` and NaN. (Note
for anyone re-deriving it: naive `searchsorted` is *wrong* for NaN, which it sorts
to the end, giving `p_low = 1.0` where the loop's all-False comparisons give the
`1/(N+1)` floor on both tails. `empirical_p_matrix` masks NaN and restores the
floor.) To re-check at any time:

```bash
python $S/run_screen.py --config screens.yaml --out /tmp/vec
python $S/run_screen.py --config screens.yaml --out /tmp/loop --legacy-loop
```

After a run, worth asking:

- Does the null QC panel show the bootstrapped null sitting on top of the
  per-replicate control histograms? A null narrower than the replicates means
  too few controls; one much wider means a replicate is misbehaving.
- Is the Fisher/Stouffer overlap close to the Fisher count? A big gap means hits
  resting on single replicates.
- Does the provenance JSON's `n_control_lfcs_per_condition` clear ~100?
- Do the figure notebooks still find their columns?
  `python $S/check_inputs.py --check-downstream LFC-FDR/X_LFC_FDR.csv --notebooks *.ipynb`
  (advisory only — f-string-built column names are invisible to it).

`new_names` / the `samps` keys are a **stable public key**. Downstream figure
code selects columns by `f'LFC_median_{cond}'`, so renaming `bm` to `bonemarrow`
breaks every figure notebook silently.
