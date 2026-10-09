# Gallery — one-off chart types from the papers

Snippets assume the header cell from SKILL.md (palettes `blues`, `greens`, …
defined) and, where noted, `import kdplot as kp`. Paths are relative to
`~/Documents/GitHub/`. Open the source notebook when you need the full context.

The font sizes, line widths and marker sizes below are copied from the original
draw-big notebooks (screen mode). At the default print size, pass them through
`kp.fs()`, `kp.lw()`, `kp.ms()`, `kp.area()` — or use 7 pt labels / 6 pt ticks —
and apply the Nature conventions from SKILL.md (legend on top, exact P, italic
gene names, ρ and R).

## Contents
1. Pie and nested donut
2. Venn diagram
3. Stacked single bar with counts
4. Count bars with rescaled axis
5. Overlaid histograms (ABE vs CBE; targeting vs controls)
6. Screen vs validation scatter with error bars
7. Pairwise correlation grid (hexbin, lower triangle)
8. Values along a protein, with domain boundaries
9. Paired score scatter (log–log thresholds; 0–1 scores with zoom box)
10. Cumulative "modeled %" curve
11. Broken-axis histogram
12. Ridge plot

---

## 1. Pie and nested donut
Source: `h2m-public/figure_generation/figure1.ipynb` (1e, 1f); `B-ALL-in-vivo-base-editing/figure3_CLEANED.ipynb` (3a).
No wedge labels — numbers go in the legend text or Illustrator.

```python
fig, ax = plt.subplots(figsize=(5, 5))
ax.axis('equal')
wedges, _ = ax.pie([n_multi, n_one, n_none], radius=1.4,
                   colors=[plt.cm.Purples(0.8), plt.cm.Purples(0.6), 'grey'],
                   explode=(0, 0.1, 0.2), startangle=0)
plt.setp(wedges, edgecolor='white')

# nested donut: outer ring = modeled vs not, inner ring = class breakdown
outer, _ = ax.pie([modeled, total - modeled], radius=1.4, colors=[plt.cm.Greens(0.7), 'white'])
plt.setp(outer, width=0.4, edgecolor='white')
inner, _ = ax.pie([n2, n0, n1, n3, rest], radius=1,
                  colors=[plt.cm.Oranges(0.8), plt.cm.Oranges(0.6), plt.cm.Oranges(0.4), 'grey', 'white'])
plt.setp(inner, width=0.4, edgecolor='white')
```

## 2. Venn diagram
Source: `B-ALL-in-vivo-base-editing/figure3_CLEANED.ipynb` (3b), `figure1_CLEANED.ipynb` (1d).

```python
from matplotlib_venn import venn2
v = venn2([og_hits, focused_hits], set_labels=('OG', 'Focused'),
          set_colors=('#f4a26188', '#457b9d88'), ax=ax)
for t in v.set_labels:
    if t: t.set_fontsize(14)
for t in v.subset_labels:
    if t: t.set_fontsize(16)
ax.set_title('Spleen', fontsize=16)
```

## 3. Stacked single bar with counts
Source: `B-ALL-in-vivo-base-editing/figure3_CLEANED.ipynb` (3a-2).

```python
palette_cat = {'In vitro only': yellows[3], 'Spleen/BM only': blues[3], 'Meninges only': reds[3],
               'In vivo only (multi)': greens[3], 'In vitro + In vivo': purples[5]}
fig, axes = plt.subplots(1, 2, figsize=(6, 5), sharey=True)
for ax, counts, title in zip(axes, [enr_counts, dep_counts], ['Enrichers', 'Depletors']):
    bottom = 0
    for cat in cat_order:
        val = counts.get(cat, 0)
        if val == 0:
            continue
        ax.bar(0, val, bottom=bottom, color=palette_cat[cat], edgecolor='none', width=0.5)
        if val > 5:
            ax.text(0, bottom + val / 2, str(val), ha='center', va='center',
                    fontsize=13, color='white', fontweight='bold')
        bottom += val
    ax.set_xlim(-0.5, 0.5); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=16)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
handles = [Patch(color=palette_cat[c], label=c) for c in cat_order]
fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.55, 0.5),
           fontsize=12, frameon=False, title='Hit category', title_fontsize=12)
```

## 4. Count bars with rescaled axis
Source: `B-ALL-in-vivo-base-editing/figure1_CLEANED.ipynb` (1b).

```python
fig, ax = plt.subplots(figsize=(8, 6))
ax.bar(range(5), heights, edgecolor='black', linewidth=1, alpha=.8,
       color=['tab:purple', 'tab:purple', 'tab:blue', 'tab:blue', 'tab:green'])
ax.set_yticks([0, 5e3, 10e3, 15e3, 20e3, 25e3]); ax.set_yticklabels([0, 5, 10, 15, 20, 25])
ax.set_ylabel('No. of Mutations (x1,000)', fontsize=16)
ax.set_xticks(range(5))
ax.set_xticklabels(['SNP in ≥ 3\nSamples', 'A>G\nor C>T', 'H2M\nMapped', 'Flank\nSize ≥ 2', 'Optimal NG\nPAM present'])
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=16)
```

## 5. Overlaid histograms
ABE vs CBE — source: `B-ALL-in-vivo-base-editing/figure1_CLEANED.ipynb` (1c):
```python
bins = np.linspace(0, 280, 29)
ax.hist(n_cbe, bins, edgecolor='black', linewidth=1, alpha=.7, color='tab:blue', label='CBE')
ax.hist(n_abe, bins, edgecolor='black', linewidth=1, alpha=.7, color='tab:purple', label='ABE')
ax.set_yscale('symlog')
ax.legend(fontsize=15)
```
Targeting vs controls with dynamic range — source: `B-ALL-in-vivo-base-editing/extended_figure13_CLEANED.ipynb` (13b):
```python
bins = np.linspace(-5, 5, 100)
ax.hist(ctrl, density=True, bins=bins, alpha=0.6, color='grey', label='Controls')
ax.hist(targ, density=True, bins=bins, alpha=0.7, color=tissue_color, label='Targeting Guides')
dr = np.percentile(targ, 95) - np.percentile(targ, 5)
ax.text(0.97, 0.97, f'DR = {dr:.2f}', transform=ax.transAxes, fontsize=10, ha='right', va='top', color=tissue_color)
ax.spines[['top', 'right', 'left']].set_visible(False)
ax.set_yticks([])
```

## 6. Screen vs validation scatter with error bars
Source: `B-ALL-in-vivo-base-editing/figure4_CLEANED.ipynb` (validation correlation).

```python
ax.errorbar(x, y, xerr=x_sd, yerr=y_sd, fmt='o', color=tissue_color, capsize=5, alpha=0.8,
            elinewidth=1, markersize=7, markeredgecolor='black', markeredgewidth=0.8)
for xi, yi, name in zip(x, y, names):
    ax.annotate(name, (xi, yi), fontsize=7, ha='left', va='bottom', xytext=(3, 3), textcoords='offset points')
ax.axhline(0, color='black', linestyle='dashed', linewidth=0.8)
ax.axvline(0, color='black', linestyle='dashed', linewidth=0.8)
ax.axline((0, 0), slope=1, color='black', linestyle='dashed', linewidth=0.8)
r, p = scipy.stats.pearsonr(x, y)
ax.set_title(f'{tissue}\n$r={r:.2f}$, $p={p:.3f}$', fontsize=12)
```

## 7. Pairwise correlation grid (hexbin, lower triangle)
Source: `B-ALL-in-vivo-base-editing/extended_figures.ipynb` (KO vs BE correlations).

```python
n = len(cols)
fig, ax = plt.subplots(n, n, figsize=(10, 10), sharex=True, sharey=True)
for i, a in enumerate(cols):
    for j, b in enumerate(cols):
        if j >= i:
            fig.delaxes(ax[i][j]); continue
        ax[i][j].hexbin(df[b], df[a], bins='log', gridsize=30, linewidth=0, extent=(-8, 8, -8, 8))
        ax[i][j].plot([-8, 6], [-8, 6], linestyle='dashed', color='black', linewidth=1)
        ax[i][j].axhline(0, linestyle='dashed', color='black', linewidth=1)
        ax[i][j].axvline(0, linestyle='dashed', color='black', linewidth=1)
        ax[i][j].set_title(kp.corr_text(df[b], df[a]), fontsize=13)
        ax[i][j].spines[['top', 'right']].set_visible(False)
        ax[i][j].tick_params(axis='both', which='major', labelsize=13)
    ax[i][0].set_ylabel(labels[a], fontsize=13)
    ax[n - 1][i].set_xlabel(labels[cols[i]], fontsize=13)
```

## 8. Values along a protein, with domain boundaries
Source: `h2m-public/figure_generation/figure2.ipynb` (2d), `extended_figure6.ipynb` (6c).
Class colors used there: `{0: 'tab:orange', 1: 'tab:green', 4: 'tab:purple'}` (H2M class); domains IG `tab:green`, TM `tab:orange`, Kinase `tab:purple`.

```python
fig, ax = plt.subplots(figsize=(12, 6))
ax.scatter(data=df, x='loc_h', y='count', c='color', alpha=.3, s=120, edgecolor='white', linewidths=.5)
for b in domain_boundaries:
    ax.plot([b, b], [-10, ymax], linewidth=1, linestyle='dashed', color='black')
ax.plot([0, protein_len], [5, 5], linewidth=2, linestyle='-.', color='black')      # recurrence threshold
ax.xaxis.set_ticks(domain_starts)
ax.set_xlabel('Codon', fontsize=20); ax.set_ylabel('# of patients', fontsize=20)
kp.label_points(ax, hot['loc_h'], hot['count'], hot['HGVSp_h'], colors=hot['color'], fontsize=18)
ax.spines[['top', 'right']].set_visible(False)
```

## 9. Paired score scatter
Log–log %rank with binding thresholds — source: `h2m-public/figure_generation/figure2.ipynb` (2g/2h), `extended_figure5.ipynb`.
Class colors: strong `darkred`, both `tab:purple`, human-only `tab:orange`, mouse-only `tab:blue`, neither `grey`.
```python
from matplotlib.ticker import FixedLocator, FuncFormatter
fig, ax = plt.subplots(figsize=(9, 9))
ax.scatter(data=diff, x='rank_m', y='rank_h', c='color', alpha=.3, s=80, edgecolor='black', linewidths=.5, marker='^')
ax.scatter(data=same, x='rank_m', y='rank_h', c='color', alpha=.5, s=120, edgecolor='black', linewidths=.5)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_aspect('equal')
ax.plot([.01, 100], [.01, 100], linewidth=2, linestyle='--', color='grey')
for t, c in [(.5, 'darkred'), (2, 'tab:purple')]:
    ax.plot([.01, t], [t, t], linewidth=2, linestyle='--', color=c)
    ax.plot([t, t], [.01, t], linewidth=2, linestyle='--', color=c)
ticks, labels = [.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100%']
for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, pos: labels[pos]))
kp.label_points(ax, sel['rank_m'], sel['rank_h'], sel['label'], colors=sel['color'], fontsize=15, box=True)
```
0–1 scores (AlphaMissense/SIFT) with zoom box — source: `h2m-public/figure_generation/figure2.ipynb` (2i, 2k).
Cutoffs `0.34`/`0.56` as grey `'-.'` lines; point color = `kp.blend(color_1, color_2)` of the two calls;
title `f'{n} Paralogous Mutation Pairs\nPearson Corr={r_p:.2f}|Spearman Corr={r_s:.2f}'`.
Mark the zoom region with a black `lw=2` rectangle; the zoom panel has all four spines `lw=3` and no ticks.

## 10. Cumulative "modeled %" curve
Source: `h2m-public/figure_generation/figure2.ipynb` (2e).

```python
from sksurv.nonparametric import kaplan_meier_estimator
for region, color in region_colors.items():
    m = df['region'] == region
    t, p, _ = kaplan_meier_estimator(df['status'][m], df['flank_size'][m], conf_type='log-log')
    ax.step(t, p, where='post', color=color, linewidth=3, label=f'{region} (n = {m.sum()})')
ax.set_ylim(0, 1)
ax.yaxis.set_ticks([.2, .4, .6, .8, 1]); ax.yaxis.set_ticklabels(['20', '40', '60', '80', '100%'])
ax.legend(loc='best', fontsize=14, frameon=False)
```

## 11. Broken-axis histogram
Source: `h2m-public/figure_generation/figure1.ipynb` (1g). Long-tailed counts: a wide main panel plus a narrow ">max" panel.

```python
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9, 6), sharey=True, gridspec_kw={'width_ratios': [20, 1]})
fig.subplots_adjust(wspace=0.05)
for a in (ax, ax2):
    a.hist(values, bins=1000, color=plt.cm.Oranges(0.8))
    kp.economist_ax(a, grid='y')
ax2.spines['left'].set_visible(False)
ax.set_xlim(0, 79); ax2.set_xlim(tail - 1, tail)
ax2.xaxis.set_ticks([tail]); ax2.xaxis.set_ticklabels([f'>{tail}'])
ax.yaxis.set_ticks([3e4, 6e4, 9e4]); ax.yaxis.set_ticklabels(['30', '60', '90K'])   # abbreviated counts
ax.plot([median, median], [0, ymax], linewidth=2, color='black')                    # median marker
d, dd = .015, .02                                                                   # break marks
kw = dict(transform=ax.transAxes, color='k', clip_on=False, linewidth=3)
ax.plot((1 - d, 1 + d), (-dd, dd), **kw); ax.plot((1 - d, 1 + d), (1 - dd, 1 + dd), **kw)
kw.update(transform=ax2.transAxes)
ax2.plot((-d * 20, d * 20), (1 - dd, 1 + dd), **kw); ax2.plot((-d * 20, d * 20), (-dd, dd), **kw)
```

## 12. Ridge plot
Source: `h2m-public/figure_generation/extended_figure4.ipynb` (4a). Stacked KDEs, one row per group.

```python
sns.set_theme(style='white', rc={'axes.facecolor': (0, 0, 0, 0)})
g = sns.FacetGrid(df, row='F', hue='F', aspect=8, height=3, palette='viridis_r', hue_order=[0, 4, 10, 20])
g.map(sns.kdeplot, 'am_pathogenicity', bw_adjust=1, clip_on=False, fill=True, alpha=1, linewidth=1.5)
g.map(sns.kdeplot, 'am_pathogenicity', bw_adjust=1, clip_on=False, color='w', lw=8)
g.map(plt.axhline, y=0, lw=8, clip_on=False)
for i, ax in enumerate(g.axes.flat):
    ax.text(-.25, .1, f'Flank = {[0, 4, 10, 20][i]}', fontsize=40, color=ax.lines[-1].get_color())
    ax.set_xlabel(''); ax.set_ylabel('')
g.figure.subplots_adjust(hspace=-.5)
g.set_titles(''); g.set(yticks=[]); g.despine(bottom=True, left=True)
```
`sns.set_theme` resets rcParams — re-run the header's rcParams afterwards so later panels keep Helvetica and fonttype 42.
