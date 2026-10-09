---
name: sci-plot
description: Kexin Dong's publication figure style for Python (matplotlib + seaborn) — fonts, open axes, black-edged marks, sequential-palette colors with fixed meanings (tissues/sample types, ABE vs CBE, COSMIC roles, enriched vs depleted, controls), PDF export with editable text, and tested recipes for screen plots (volcano, waterfall, bar + replicate dots, box + exact-P brackets, correlation scatter, z-score heatmap, correlation matrix, sample PCA, h2m-style bars). Defaults to Nature-style panels at final print size (mm figure widths, 7 pt text, bars with SD caps and dots, legends above the plot, italic gene names, group and fold-change/exact-P brackets, dot + mean plots, log hexbin density, per-exon heatmaps, protein domain schematics, aligned panel letters). Learned from the B-ALL in vivo base-editing and H2M paper figure notebooks and two published Nature-family prime-editing figures. Use this skill whenever writing, fixing, or restyling any Python plot, figure panel, or figure notebook for Kexin's research — "plot the LFCs", "make figure 3c", "volcano for bone marrow", "add an extended figure panel", "clean this up for the paper", "make it look like a Nature figure", "final size for submission", "make it match my other figures" — even when style isn't mentioned. Prefer it over generic chart-styling guidance for static matplotlib/seaborn figures; not for interactive dashboards or HTML pages.
---

# Scientific plotting in Kexin's style

These conventions come from the published figure notebooks in
`~/Documents/GitHub/B-ALL-in-vivo-base-editing/*_CLEANED.ipynb` and
`~/Documents/GitHub/h2m-public/figure_generation/*.ipynb`, plus the panel
conventions of two published Nature-family prime-editing figures (Nat Biotech
2024 PEGG/sensor paper, `~/Downloads/s41587-024-02172-9.pdf`; Nature 2024 La
paper, `~/Downloads/s41586-024-07259-6.pdf`). The goal is panels that
drop straight into Illustrator and look like they belong to the same paper:
open axes, black-edged marks, colors that mean the same thing in every panel,
and PDFs whose text stays editable.

## Bundled files

- `scripts/kdplot.py` — tested implementation of everything below: palettes,
  semantic color maps, `setup()`, axis helpers, `save()`, and one function per
  recurring plot type. Read the docstring/body of the function you need; it is
  the reference implementation.
- `references/gallery.md` — one-off chart types from the papers (pies/donuts,
  Venns, stacked counts, histogram overlays, error-bar validation scatter,
  pairwise hexbin grid, protein-position scatter, log–log threshold scatter,
  cumulative curves, broken axes, ridge plots), each with its source notebook.
  Read it when the request isn't covered by a kdplot function.

## Workflow

1. **Decide where the code lives.**
   - *Figure notebook in a repo* (anything headed for a paper or shared with
     co-authors): keep it self-contained so it reruns for anyone who clones the
     repo. Put the header cell below at the top (or reuse the one already there),
     and write plotting code inline by adapting the matching kdplot function.
     Don't import from `~/.claude` in these notebooks — collaborators don't have it.
     If several panels need the same helper, offer to copy `kdplot.py` into the
     repo instead of duplicating code.
   - *Quick look, exploration, scratch script*: import kdplot directly —
     ```python
     import sys, os
     sys.path.insert(0, os.path.expanduser('~/.claude/skills/sci-plot/scripts'))
     import kdplot as kp
     kp.setup()                  # final print size by default — see "Two sizes"
     ```
   - **Defaults for everything new**: final print size (journal mode) and the
     Nature panel conventions below. The one exception: when adding a panel to
     an existing B-ALL / h2m figure drawn in the older draw-big style, match that
     figure (`kp.setup(mode='screen')`, stars, `$R_P$ | $R_S$`) so panels in one
     figure stay consistent — and say you did.
2. **Pick the recipe** from the table at the bottom (or the gallery).
3. **Apply the style rules** below.
4. **Save** as PDF (`bbox_inches='tight'`).
5. **Look at the result before calling it done.** Run the code, write a PNG
   (`dpi=150` is enough to inspect) and view it. adjust_text collisions,
   legends sitting on data, clipped tick labels, and colors that don't match the
   rest of the figure only show up in the render. Fix them — widen limits, drop
   or shorten labels, move the legend — then re-check.

## Header cell

For new notebooks use this (print size, Nature look); in an existing notebook
keep the header it already has:

```python
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import scipy.stats
from adjustText import adjust_text
import warnings
import os
warnings.filterwarnings('ignore')

# ── Global style: final print size, Nature look ───────────────────────────
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['text.usetex']  = False
plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']
mpl.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Helvetica',   # italic gene names,
                     'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})  # P = 6.7 × 10$^{-3}$
mpl.rcParams.update({'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
                     'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
                     'axes.linewidth': .6, 'lines.linewidth': .75, 'patch.linewidth': .5,
                     'xtick.major.width': .6, 'ytick.major.width': .6,
                     'xtick.major.size': 3, 'ytick.major.size': 3,
                     'errorbar.capsize': 2, 'legend.frameon': False, 'figure.dpi': 200})
MM = 1 / 25.4   # figsize=(183 * MM, 110 * MM): 89 mm single column, 183 mm double

blues   = sns.color_palette('Blues').as_hex()
greens  = sns.color_palette('Greens').as_hex()
reds    = sns.color_palette('Reds').as_hex()
yellows = sns.color_palette('Oranges').as_hex()
purples = sns.color_palette('Purples').as_hex()
greys   = sns.color_palette('Greys').as_hex()
teals   = sns.color_palette('BuPu').as_hex()
BASE = '/Users/kexindong/Documents/GitHub/<repo>/'
```

`pdf.fonttype = 42` embeds TrueType so text is editable in Illustrator — never
drop it. h2m notebooks use `['Arial']` instead of Helvetica.

## Style rules

### Axes
- Open frame on every Cartesian plot: `ax.spines[['top', 'right']].set_visible(False)`. No gridlines (exception: the h2m bar style below).
- Ticks outward: at print size `ax.tick_params(axis='both', which='major', labelsize=6, direction='out', length=3, width=.6)` (screen mode: 16 / 6 / 1).
- Use `fig, ax = plt.subplots(...)` and explicit `ax.` methods; `fig.tight_layout()` before saving.
- Blank the categorical axis label (`ax.set_ylabel('')`) — the tick labels already say it. Break long tick labels with `\n` (`'SNP in ≥ 3\nSamples'`) rather than rotating.
- Put scale in the label and rescale ticks for big counts: ticks `[0, 5e3, …]` labelled `[0, 5, …]`, label `'No. of Mutations (x1,000)'`.
- Shared axes for panel rows (`sharey=True`); label the outer axes only.

### Font sizes

| Element | **Print size (default)** | Screen: single panel (~4–6 in) | Screen: panel row (12–20 in) |
|---|---|---|---|
| Axis labels | 7 | 16 | 12–14 |
| Tick labels | 6 | 16 | 11–13 |
| Titles | 7, plain | 16 | 14 |
| Legend | 6–7 (`frameon=False`) | 12 | 10 |
| Point labels | 5–6 | 10–12 | 7–9 |
| Heatmap annotations | 5–6 | 11–13 | 10 |
| Panel letters | 8, bold | 22, bold | 22, bold |

Nothing below 5 pt at print size. h2m figures used bigger screen canvases with
proportionally bigger fonts (20–40); only go that large when matching an
existing h2m panel.

### Two sizes: screen vs journal

- **Journal** (default, `kp.setup()`): draw at final print size, as the
  Nature / Nat Biotech panels are — 7 pt labels, 6 pt ticks, 5 pt minimum,
  0.6 pt axes, 3 pt ticks, 8 pt bold panel letters. Size the figure in mm:
  `plt.subplots(figsize=kp.mm(183, 110))` (89 mm single column, 120 mm 1.5
  column, 183 mm double) and save the PDF at that size — no rescaling in
  Illustrator.
- **Screen** (`kp.setup(mode='screen')`): draw big and scale down in
  Illustrator — how the B-ALL and h2m notebooks were made. Only for adding a
  panel to one of those existing figures, or when the user asks for it.

Every kdplot recipe scales its fonts, lines and markers through `kp.fs()`,
`kp.lw()`, `kp.ms()`, `kp.area()`, so the same call works in both modes; use
those helpers for any size you hard-code yourself. In a self-contained notebook
the header cell above already carries the print-size rcParams — write explicit
sizes as 7 (labels) / 6 (ticks) pt.

### Colors

Colors come from seaborn sequential palettes as 6-step hex lists
(`blues[0]` lightest … `blues[5]` darkest; `[3]` is the usual fill, `[2]`
highlights, `[4]/[5]` dark accents) or, in h2m, colormap shades such as
`plt.cm.Oranges(0.8)`. Colors carry meaning, so reuse these maps rather than
inventing new ones — a reader should be able to recognize "bone marrow" or
"enriched" by color across the whole paper:

**Categories with no fixed meaning** (conditions, groups, arbitrary levels) come
from a categorical panel — `kp.categorical(n, style)`, n = 1..8. Both canonical
sets are built only from the sequential palettes above, so a categorical panel
never introduces a hue the paper does not already use.

| Style | Shades | Use | n=4 |
|---|---|---|---|
| `normal` **(default)** | alternating [3] / [4] | most panels; differs in value as well as hue, so it survives greyscale | `#539ecd` `#e05206` `#2c944c` `#8d89c0` |
| `soft` | [2], same hue order | dense panels, many overlapping marks, anywhere `normal` reads too loud | `#89bedc` `#fda057` `#8ed08b` `#afaed4` |
| `nature` | sampled from published Nature figures | muted and slightly desaturated — nothing fully saturated, which is what keeps those panels calm at eight categories; best separation here (dE 27) | `#3c6c9c` `#d86c3c` `#9cc090` `#9c3054` |
| `nature_soft` | the pale twin | panel backgrounds, light fills | `#a8c0d8` `#d8b484` `#9cc090` `#6c84b4` |
| `nature_paired` | light/dark twins per hue | **nested conditions** where the pairing is the point (D16 / D16-NUT, control / treated within a genotype); use `[0]+[1]`, `[2]+[3]`, … | `#a8c0d8` `#3c6c9c` `#9cc090` `#488454` |

Full order for both: blue, orange, green, purple, red, pink, grey, then violet
(`normal`) or amber (`soft`). `soft`'s eighth is an amber rather than a second
periwinkle — `teals[2]` and `purples[2]` are only dE 8 apart, not separable.

`PALETTES` also carries alternatives kept for reference, should a figure need
them: `deep` (all [4], for thin lines), `optimised` (highest separation, dE 40),
`warm` (orange/pink led), `okabe` (Okabe-Ito) and `tol_muted` (Paul Tol).

Don't hand-mix a new categorical set from `[3]` of every sequential palette:
`Oranges[3]` ≈ `YlOrBr[3]` and `Purples[3]` ≈ `BuPu[3]` (dE 6.5 — the same colour
to the eye).

**Categories that DO carry meaning** use the fixed maps below, so a reader
recognizes "bone marrow" or "enriched" by color across the whole paper:

| Meaning | Colors (`kdplot` constant) |
|---|---|
| Sample type, in this order | Plasmid `black`, Input `yellows[3]`, D5 `purples[2]`, D15 `purples[5]`, Spleen `blues[3]`, Bone Marrow `greens[3]`, Meninges `reds[3]` (`SAMPLE_ORDER`, `SAMPLE_COLORS`) |
| Editor | ABE `tab:purple`, CBE `tab:blue` (`EDITOR_COLORS`) |
| COSMIC role | Oncogene `tab:red`, Oncogene/TSG `tab:purple`, TSG `tab:blue`, Undefined `tab:grey` (`COSMIC_COLORS`) |
| Hit direction | enriched `reds[2]`, depleted `blues[2]`, not significant `lightgrey` (`DIRECTION_COLORS`) |
| Guide class | targeting `blues[3]`, controls `#e41a1c` in scatters / `grey` in histograms (`GUIDE_COLORS`) |
| Editing outcome | target w/ bystanders `blues[4]`, pure correct `teals[2]` (`EDIT_TYPE_COLORS`) |
| Method overlap | Both `#e63946`, A only `#f4a261`, B only `#457b9d`, Neither `#cccccc` (`OVERLAP_COLORS`) |
| Pathogenicity (h2m) | pathogenic `darkred`, benign `darkblue`, ambiguous `grey`; pair colors = `blend()` of both calls |
| Edit outcome, Nature look | intended edit `#e38683` (salmon), indels `#b9b9b9` (`OUTCOME_COLORS`) |
| Protein schematics | domains `#adc2dc` / `#7189b7`, backbone `#9f9f9f`, tag `#b79a72`, mutation `#ba494b` (`SCHEMATIC_COLORS`) |
| One dot color per target site | `#e58782` `#7598cb` `#813b7a` `#473d89` (`categorical(4, 'nature_sites')`) |

The sample-type colors follow Figure 2 (b, c). A few older notebooks (Fig 4,
Ext Fig 10/13) swapped spleen and bone marrow — don't copy those. Grey always
means background / not significant / control. For categories with no mapping,
pick explicitly from the `tab:` colors or a palette shade; never rely on the
default color cycle.

### Marks
- **Grouped data → bar + scatter, always.** Whenever replicates exist this is the
  default, not an embellishment: the bar carries the mean, the points carry n and
  spread, and a bar alone hides both. Default recipe: `kp.bar_dots()` — mean bar,
  black SD error bar with caps, one dot per replicate, bar positions returned for
  brackets. `kp.bar_strip()` (seaborn, settings below) remains for horizontal
  bars, log axes and the older B-ALL panels.
  - bar: `edgecolor='black', linewidth=1`, `zorder=1`
  - points: every replicate, jittered, **`alpha=.55` by default** so overlapping
    replicates stay readable — drop to .3–.4 when n is large; thin black edge
    (`linewidth=.5`), `zorder=3` so they sit above the bar
  - error bar: SD or 95% CI, black, `lw=1`, `capsize=3`
  - `point_color='black'` matches the published Nature panels; leave it `None` to
    colour points by category

  Plain count bars (no replicates) use `alpha=.8` and no points.
- **Scatter**: dense clouds small and soft (`s=10–30`, `alpha=.3–.5`, `linewidths=0`); points that matter bigger with a black edge (`s≈50`, `edgecolors='black', linewidths=.5`, higher `zorder`). Size can encode significance (waterfall FDR bins → `s` 10/40/80/160).
- **Box plots**: `showfliers=False` plus every point on top (the points already show outliers).
- **Reference lines**: thin dashed — zero line and y = x in black/grey `lw=.8–1`, thresholds (FDR 0.1, |LFC| 1) grey `'--'`; h2m score cutoffs use `'-.'`.
- **Heatmaps**: LFC/z-scores → `cmap='seismic'`, symmetric `vmin=-5, vmax=5`, `annot=True`, `linewidth=1, linecolor='white'`. Correlations → `cmap='Reds'`, `vmin=0, vmax=1`, annotated, no clustering.
- **Pies/donuts**: no text on the wedges; `edgecolor='white'`, `explode` to separate classes, colormap shades.

### Traps worth knowing

Each of these produced a wrong or unreadable panel in practice.

- **Box-plot whiskers are not the maximum.** With `showfliers=False` the visible
  top of a box is Q3 + 1.5·IQR. A right-skewed group can therefore look *lower*
  than a tight one whose true maximum is smaller — the opposite of the truth. If
  a reader might compare maxima, mark the true max (`marker='v'`, black, above
  the box) or show the points.
- **One bin per integer breaks at scale.** `np.arange(0, max+2) - 0.5` is fine
  for counts in the tens; at 10⁴–10⁵ it puts tens of thousands of sub-pixel bins
  on a 3-inch axis — the bars disappear and the PDF balloons (33 MB in one case,
  0.07 MB after the fix). When the range spans orders of magnitude use
  `np.logspace(log10(.5), log10(hi*1.5), 25)` with `set_xscale('log')`, and
  `.clip(lower=.5)` so zero-count categories stay visible on the log axis.
- **Mouse vs human gene symbols.** Mouse symbols are title case, and *Trp53 ≠
  TP53*. Upper-casing alone leaves Trp53 unmatched, so the gene is silently
  dropped from every panel that groups by gene — no error, just a smaller count.
  Map with `kp.human_symbol()` before matching, and assert nothing lands in
  "other".
- **Rank plots need their sort key stated.** When ranking per mutation or gene,
  say what the order is in the axis label ("ordered by best bystander score"),
  or the reader assumes it is the plotted quantity.
- **Look at the render.** adjust_text collisions, legends sitting on data,
  clipped rotated tick labels and near-identical colours only appear in the
  image. Write a PNG at `dpi=150`, open it, fix, repeat.

### Text & annotation
- Point labels via `adjust_text` with thin black leader lines (`arrowprops=dict(arrowstyle='-', color='black', lw=.5)`), each label colored like its point. Mutation labels read '*Kras* G12D' — gene in italics; fall back to the human change with `' (H)'` (`kdplot.mutation_label`, `italic_gene=False` for plain).
- Correlations: always report both — `'ρ = 0.83 and R = 0.84'` in the title (`kdplot.corr_text`; `style='kd'` gives the older `'$R_P$ = 0.84 | $R_S$ = 0.83'`).
- Comparisons: exact P on a bracket, not stars — statannotations with `test='t-test_ind'`, `comparisons_correction='bonferroni'`, `loc='outside'`, labels from `kp.format_p` (`kp.box_strip(..., pairs=…)` does this; `stars=True` for the old look). Between two bars, fold change + P via `kp.compare_bracket`. Name the test in the figure legend.
- Legends: above the plot, one row of flat squares, no frame (`kp.legend_top`). Hide seaborn's auto legend (`legend=False`). Black-edged swatches (`kp.legend_patches` or `legend_top(edge=True)`) only for outlined bars.
- Log axes: `LogLocator(base=10)` + `LogFormatterMathtext()`; powers as math text (`'FDR < $10^{-2}$'`).
- Percent axes: tick labels like `['20', '40', '60', '80', '100%']` (unit on the last tick only).

### Nature panel conventions (default)

What makes the two reference figures read as Nature panels. These are the
defaults for every new figure; drop them only to match an existing draw-big
B-ALL / h2m figure or when the user asks for something else:

- **Legend above the plot**, one row of flat square swatches, no frame
  (`kp.legend_top`). Nested conditions: black-edged swatches stacked top-left.
- **Gene names in italics** inside labels and legends: `f"{kp.italic('DNMT1')} +5 G-to-T"`.
- **Bars**: mean bar, black SD error bar with caps, every replicate as a dot
  (`kp.bar_dots`). Flat fill + black dots, or black outline + dots in the bar
  color with a black rim. Axis starts at 0.
- **Groups of bars get a bracket label** above them ('Parental cells') instead of
  a legend or a facet title (`kp.group_bracket`); a second tier of labels under
  the ticks names higher-level groups (`kp.axis_groups`).
- **Comparisons are written out, not starred**: fold change and exact P on a
  bracket joining the two bars — `0.36×, P = 6.7 × 10⁻³` (`kp.compare_bracket`
  + `kp.fold_p_text`). P values: italic *P*, two decimals down to 0.01, then
  `× 10ⁿ` (`kp.format_p`). State the test in the legend text.
- **Long category names rotate 45°**, right-anchored (`rotation_mode='anchor'`);
  short ones stay horizontal.
- **Summary across sites**: one dot per site in its own color, a short black
  line at the mean (`kp.dot_mean`).
- **Thousands of points**: hexbin with log-scaled counts in `inferno`, slim
  color bar labelled on its right; title gives `ρ = … and R = …`
  (`kp.hexbin_density`).
- **Position along a gene**: per-exon small-multiple heatmaps, + and − strand
  rows, widths proportional to exon length; counts in Blues, efficiency in Reds
  (`kp.exon_heatmap`). **Constructs**: grey backbone, blue domain blocks, red
  mutation bars (`kp.domain_diagram`).
- **Panel letters**: bold lowercase, 8 pt in print, aligned per row at the top-left
  of each panel's full extent (`kp.panel_labels(axs)` after `tight_layout`).
- Titles are plain (not bold) and centered; axis labels in sentence case with
  units in parentheses — 'Average correct editing (%)'.

### Saving
- `fig.savefig(BASE + 'figures/figure4/<name>.pdf', bbox_inches='tight')`; `os.makedirs(..., exist_ok=True)` first. PDF is the default (vector, editable text); PNG only for previews/slides at `dpi=300`; `transparent=True` only when asked.
- Name files after the panel (`2e.pdf`, `ext_4c.pdf`) or its content (`volcano_leukemia_ABE_bm.pdf`).
- In notebooks, `plt.show()` after saving.

### Notebook layout
- First markdown cell: `# Figure N`, authors, `Last Updated: <Month D, YYYY>`, one-line figure title.
- Each panel: a markdown cell `## 2b` with a one-sentence legend of what is shown and how it's colored, then the code cell.
- Define cutoffs as named constants at the top of the cell (`min_edit = 20`, `FDR_cutoff = 0.1`, `LFC_cutoff = 1.0`) and use `# ── Section ───` banners inside long cells.

## Recipes

| Plot | kdplot function | Original panel |
|---|---|---|
| **Grouped bars (default)**: mean + SD caps + replicate dots | `bar_dots` | Nature 2024 La Fig b, d; Nat Biotech 2024 PEGG Fig b |
| Seaborn bar + replicate dots: horizontal, log axis, older B-ALL look | `bar_strip` | B-ALL Fig 2b, 2d; Fig 4 validation |
| Box + all points + exact-P brackets (`stars=True` for stars) | `box_strip` | h2m Ext Fig 4c/4d; B-ALL Fig 4 |
| Volcano (median LFC vs −log10 FDR) | `volcano` | B-ALL Fig 2f |
| Waterfall (guides ranked by LFC, sized by FDR) | `waterfall` | B-ALL Fig 2e |
| Paired correlation scatter (screens, replicates) | `corr_scatter` | B-ALL Ext Fig 13c |
| LFC / z-score heatmap | `z_heatmap` | B-ALL Fig 3c/3d, Fig 4 barcoding |
| Sample correlation matrix | `corr_matrix` | B-ALL Ext Fig 8, 14 |
| Sample PCA from counts | `pca_samples` | B-ALL Fig 2c |
| h2m pale horizontal bars with grid | `pale_barh` (+ `economist_ax`) | h2m Fig 1h, 2g |
| Group bracket / comparison bracket / 2nd-tier x labels | `group_bracket`, `compare_bracket` (+ `fold_p_text`, `format_p`), `axis_groups` | Nature 2024 La Fig b, d |
| Horizontal bars + grey replicate dots, zero line | `bar_strip` (`point_color='grey'`) + `ax.axvline(0)` | Nat Biotech 2024 PEGG Fig a |
| Dot per site + mean line | `dot_mean` | Nature 2024 La Fig e |
| Log hexbin density + colour bar | `hexbin_density` | Nat Biotech 2024 PEGG Fig c |
| Per-exon strand heatmaps | `exon_heatmap` (stack two on `fig.subfigures(2, 1)`) | Nat Biotech 2024 PEGG Fig d |
| Protein / construct schematic | `domain_diagram` | Nature 2024 La Fig a |
| Panel letters, figure sized in mm | `panel_labels`, `mm` | both Nature figures |
| Anything else | `references/gallery.md` | — |

**h2m flavor** — when working in `h2m-public` or asked for the H2M look: Arial,
colors as colormap shades (`plt.cm.Oranges(0.8/0.6/0.4)`, pale fill
`cm(0.2)` with dark outline `cm(0.9)`), and for bars/histograms the
`economist_ax` treatment: `#A8BAC4` gridlines behind the data, only a heavy
butt-capped left spine, value labels on top for horizontal bars.
