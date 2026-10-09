"""kdplot — Kexin Dong's publication figure style (matplotlib + seaborn).

Distilled from the published figure notebooks in
    B-ALL-in-vivo-base-editing/*_CLEANED.ipynb   (base-editing screens; Helvetica)
    h2m-public/figure_generation/*.ipynb          (H2M paper; Arial)
and from two published Nature-family prime-editing figures
    Nat Biotech 2024 (PEGG / sensor library, s41587-024-02172-9)
    Nature 2024 (La promotes prime editing, s41586-024-07259-6)

Two ways to use it:
  * quick looks / scratch scripts — import it:
        import sys, os
        sys.path.insert(0, os.path.expanduser('~/.claude/skills/sci-plot/scripts'))
        import kdplot as kp
        kp.setup()                   # print size (journal) by default; mode='screen' to draw big
  * repo figure notebooks — keep the notebook self-contained: paste the header
    cell from SKILL.md and copy the body of the function you need from here.
"""
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import seaborn as sns
import scipy.stats

# ── Palettes ─────────────────────────────────────────────────────────────────
# Six-step seaborn sequential palettes as hex, index 0 = lightest … 5 = darkest.
# Workhorse shades: [2] highlights, [3] category fills, [4]/[5] dark accents.
blues   = sns.color_palette('Blues').as_hex()
greens  = sns.color_palette('Greens').as_hex()
reds    = sns.color_palette('Reds').as_hex()
yellows = sns.color_palette('Oranges').as_hex()
purples = sns.color_palette('Purples').as_hex()
greys   = sns.color_palette('Greys').as_hex()
teals   = sns.color_palette('BuPu').as_hex()
pinks   = sns.color_palette('RdPu').as_hex()
browns  = sns.color_palette('YlOrBr').as_hex()

# Sample / tissue types in plotting order (Figure 2b/2c convention).
SAMPLE_ORDER = ['Plasmid', 'Input', 'D5', 'D15', 'Spleen', 'Bone Marrow', 'Meninges']
SAMPLE_COLORS = {
    'Plasmid':     'black',
    'Input':       yellows[3],
    'D5':          purples[2],
    'D15':         purples[5],
    'Spleen':      blues[3],
    'Bone Marrow': greens[3],
    'Meninges':    reds[3],
}
# raw count-table column prefixes → display label (longest prefix wins)
_SAMPLE_PREFIX = {
    'plasmidlib': 'Plasmid', 'plasmid': 'Plasmid', 'lib': 'Plasmid',
    'input': 'Input',
    'd15': 'D15', 'd5': 'D5',
    'spleen': 'Spleen',
    'bonemarrow': 'Bone Marrow', 'bm': 'Bone Marrow',
    'meninges': 'Meninges', 'men': 'Meninges',
}

# ── Categorical panels ───────────────────────────────────────────────────────
# For categories that carry no fixed meaning (conditions, arbitrary groups). Both
# canonical sets are built ONLY from the colours above — the sequential [1]-[4]
# working shades — so a categorical panel never introduces a hue the paper does
# not already use.
#
# NORMAL is the default. It alternates [3] and [4] shades, so neighbouring
# categories differ in value as well as hue and stay separable in greyscale.
PALETTE_NORMAL = ['#539ecd',   # blues[3]
                  '#e05206',   # yellows[4]  (orange)
                  '#2c944c',   # greens[4]
                  '#8d89c0',   # purples[3]
                  '#d52221',   # reds[4]
                  '#e84a9b',   # pinks[3]
                  '#5c5c5c',   # greys[4]
                  '#894da3']   # teals[4]    (violet)
NORMAL_NAMES   = ['blue', 'orange', 'green', 'purple', 'red', 'pink', 'grey', 'violet']

# SOFT is the low-contrast twin: the same hue order one or two steps lighter.
# For dense panels, many overlapping marks, and anywhere NORMAL reads too loud.
# The eighth colour is a soft amber rather than a second periwinkle — teals[2]
# and purples[2] are only dE 8 apart, which is not separable.
PALETTE_SOFT = ['#89bedc',   # blues[2]
                '#fda057',   # yellows[2]  (orange)
                '#8ed08b',   # greens[2]
                '#afaed4',   # purples[2]
                '#fc8161',   # reds[2]
                '#f988ad',   # pinks[2]
                '#adadad',   # greys[2]
                '#feda7e']   # browns[1]   (amber)
SOFT_NAMES   = ['blue', 'orange', 'green', 'purple', 'red', 'pink', 'grey', 'amber']

# Kept for reference — alternatives considered, in case a figure needs them.
# 'deep'      all [4] shades: saturated, for thin lines and small markers
# 'optimised' greedy max-min over every colour in use: highest separation (dE 40)
#             but an order that does not follow the usual blue-orange-green
# 'warm'      orange/pink led, blues held back
# 'okabe'     Okabe & Ito, Color Universal Design (Nature Methods' recommendation)
# 'tol_muted' Paul Tol's muted scheme, the largest colour-blind-safe set
PALETTE_DEEP      = ['#2b7bba', '#e05206', '#2c944c', '#705eaa',
                     '#d52221', '#bc1085', '#5c5c5c', '#894da3']
PALETTE_OPTIMISED = ['#539ecd', '#e05206', '#bc1085', '#56b567',
                     '#572c92', '#feb23f', '#f988ad', '#5c5c5c']
PALETTE_WARM      = ['#f67824', '#539ecd', '#e84a9b', '#2c944c',
                     '#572c92', '#feb23f', '#828282', '#d52221']
PALETTE_OKABE     = ['#0072B2', '#F0E442', '#009E73', '#D55E00',
                     '#CC79A7', '#999999', '#E69F00', '#56B4E9']
PALETTE_TOL_MUTED = ['#332288', '#999933', '#88CCEE', '#CC6677',
                     '#117733', '#AA4499', '#44AA99', '#882255']

# NATURE — sampled from published Nature prime-editing figures. Muted and
# slightly desaturated: nothing is fully saturated, which is what makes those
# panels read as calm even at eight categories. Best separation of the sets here
# (dE 27 at n=8).
PALETTE_NATURE = ['#3c6c9c',   # steel blue
                  '#d86c3c',   # rust
                  '#9cc090',   # moss
                  '#9c3054',   # wine
                  '#e4a8a8',   # blush
                  '#a8c0d8',   # mist
                  '#302454',   # indigo
                  '#828282']   # grey
NATURE_NAMES   = ['steel', 'rust', 'moss', 'wine', 'blush', 'mist', 'indigo', 'grey']

# NATURE-SOFT — the pale twin, as used for panel backgrounds in those figures.
PALETTE_NATURE_SOFT = ['#a8c0d8', '#d8b484', '#9cc090', '#6c84b4',
                       '#f0c0c0', '#adadad', '#d89c90', '#e4a8a8']

# NATURE-PAIRED — light/dark twins of one hue, for nested conditions where the
# pairing is the point (D16 / D16-NUT, control / treated within a genotype).
# Use in pairs: [0]+[1] are one hue, [2]+[3] the next, and so on.
PALETTE_NATURE_PAIRED = ['#a8c0d8', '#3c6c9c',   # blue  light / dark
                         '#9cc090', '#488454',   # green light / dark
                         '#e4a8a8', '#b43c3c',   # red   light / dark
                         '#d8b484', '#d86c3c']   # amber light / dark

# NATURE-SITES — one dot colour per target site / gene in a dot plot
# (Nature 2024 La Fig e: DNMT1, CXCR4, VEGFA, RUNX1).
PALETTE_NATURE_SITES = ['#e58782', '#7598cb', '#813b7a', '#473d89']

PALETTES = {'normal': PALETTE_NORMAL, 'soft': PALETTE_SOFT,
            'nature': PALETTE_NATURE, 'nature_soft': PALETTE_NATURE_SOFT,
            'nature_paired': PALETTE_NATURE_PAIRED, 'nature_sites': PALETTE_NATURE_SITES,
            'deep': PALETTE_DEEP, 'optimised': PALETTE_OPTIMISED,
            'warm': PALETTE_WARM, 'okabe': PALETTE_OKABE,
            'tol_muted': PALETTE_TOL_MUTED}


#Mouse orthologue symbols are title case (Trp53, Kras) where human are upper
#(TP53, KRAS), and a few differ outright -- mouse Trp53 vs human TP53. Comparing
#gene lists case-insensitively is not enough; Trp53 silently falls through and
#the gene is dropped from every panel that groups by it.
MOUSE_TO_HUMAN = {'TRP53': 'TP53'}


def human_symbol(sym):
    """Human gene symbol for a mouse or human symbol.

    Use before matching against a gene list: upper-casing alone leaves Trp53
    unmatched against TP53, which quietly drops the gene from grouped panels."""
    s = str(sym).upper()
    return MOUSE_TO_HUMAN.get(s, s)


def categorical(n, style='normal'):
    """First `n` colours of a categorical panel (n = 1..8).

    style='normal' (the default) or 'soft' (its low-contrast twin). The other
    entries in PALETTES are alternatives kept for reference.

    For categories that DO carry meaning — tissues, editors, COSMIC roles — use
    the semantic maps below instead, so a colour keeps its meaning across the
    paper."""
    try:
        pal = PALETTES[style]
    except KeyError:
        raise ValueError('style must be one of %s' % sorted(PALETTES))
    if not 1 <= n <= len(pal):
        raise ValueError('n must be 1..%d' % len(pal))
    return pal[:n]


EDITOR_COLORS        = {'ABE': 'tab:purple', 'CBE': 'tab:blue'}
COSMIC_COLORS        = {'Oncogene': 'tab:red', 'Oncogene/TSG': 'tab:purple',
                        'TSG': 'tab:blue', 'Undefined': 'tab:grey'}
DIRECTION_COLORS     = {'enriched': reds[2], 'depleted': blues[2], 'ns': 'lightgrey'}
GUIDE_COLORS         = {'targeting': blues[3], 'control': '#e41a1c'}
EDIT_TYPE_COLORS     = {'Target Editing (w/ Bystanders)': blues[4],
                        'Pure Correct Editing': teals[2]}
OVERLAP_COLORS       = {'Both': '#e63946', 'A only': '#f4a261',
                        'B only': '#457b9d', 'Neither': '#cccccc'}
PATHOGENICITY_COLORS = {'pathogenic': 'darkred', 'benign': 'darkblue', 'ambiguous': 'grey'}
GRID_COLOR           = '#A8BAC4'   # h2m bar-chart gridlines
# Nature 2024 La Fig b/d: intended edit vs indels, side by side in every group.
OUTCOME_COLORS       = {'edit': '#e38683', 'indels': '#b9b9b9'}
# Nature 2024 La Fig a: protein / construct schematics.
SCHEMATIC_COLORS     = {'light': '#adc2dc', 'dark': '#7189b7', 'backbone': '#9f9f9f',
                        'tag': '#b79a72', 'mutation': '#ba494b'}

# Waterfall: significant guides sized by FDR bin.
FDR_BINS  = [1e-1, 1e-2, 1e-3, 1e-4]
FDR_SIZES = [10, 40, 80, 160]


# ── Size modes ───────────────────────────────────────────────────────────────
# 'journal' (DEFAULT): draw at final print size, like the Nature / Nat Biotech
#   panels — 5-7 pt text, 0.6 pt axes, 89 mm (single) / 183 mm (double) widths.
# 'screen': draw big — 16 pt text on a 4-6 in panel — and scale down in
#   Illustrator. How the B-ALL and h2m notebooks were made; use it only to add
#   a panel to one of those existing figures.
# Every size in the recipes goes through fs() / lw() / ms() / area(), so one
# setup(mode=...) call switches all of them. Call setup() before plotting.
_STYLE = {'mode': 'journal', 'text': 7 / 16, 'line': 0.6}
MM = 1 / 25.4
COLUMN_MM = {'single': 89, 'onehalf': 120, 'double': 183}
_JOURNAL_RC = {
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'axes.linewidth': .6, 'lines.linewidth': .75, 'patch.linewidth': .5,
    'xtick.major.width': .6, 'ytick.major.width': .6, 'xtick.minor.width': .4, 'ytick.minor.width': .4,
    'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
    'xtick.major.pad': 2, 'ytick.major.pad': 2, 'axes.labelpad': 2, 'axes.titlepad': 4,
    'errorbar.capsize': 2, 'legend.frameon': False, 'legend.handlelength': 1,
    'figure.dpi': 200,
}


def fs(pt):
    """Font size: `pt` in screen mode. In journal mode the screen sizes snap to
    the three print sizes so every panel matches: 14 pt and up (axis labels,
    titles) → 7, 11-13 (ticks, legends, heatmap values) → 6, smaller
    (point labels, dense annotations) → 5."""
    if _STYLE['mode'] != 'journal':
        return pt
    return 7 if pt >= 14 else 6 if pt >= 11 else 5


def lw(w):
    """Line width in the current mode."""
    return w * _STYLE['line']


def ms(d):
    """Marker diameter / tick length / small offset (points) in the current mode."""
    return d * _STYLE['text']


def area(s):
    """Scatter marker area (`s`, points²) in the current mode."""
    return s * _STYLE['text'] ** 2


def _figsize(w, h):
    t = _STYLE['text']
    return (w * t, h * t)


def mm(width, height):
    """figsize from millimetres: plt.subplots(figsize=kp.mm(89, 60)).
    Nature widths: 89 mm single column, 120 mm 1.5 column, 183 mm double."""
    return (width * MM, height * MM)


# ── Setup & colour helpers ───────────────────────────────────────────────────
def setup(font='Helvetica', mode='journal'):
    """Global rcParams from every figure notebook: TrueType (editable) text in
    PDF/PS and a sans-serif font. font='Helvetica' (B-ALL) or 'Arial' (h2m).

    mode='journal' (default) draws at final print size; mode='screen' draws big
    for scaling down in Illustrator (see Size modes). Math text — italic
    gene names via italic(), '× 10$^{-3}$' in P values — is set in Arial
    italic/bold: matplotlib can only read the regular face of macOS's
    Helvetica.ttc, and Arial's letterforms match it at print size."""
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype']  = 42
    mpl.rcParams['text.usetex']  = False
    plt.rcParams['font.family']     = 'sans-serif'
    plt.rcParams['font.sans-serif'] = list(dict.fromkeys([font, 'Helvetica', 'Arial', 'DejaVu Sans']))
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = font
    mpl.rcParams['mathtext.it'] = 'Arial:italic'
    mpl.rcParams['mathtext.bf'] = 'Arial:bold'
    if mode == 'journal':
        _STYLE.update(mode='journal', text=7 / 16, line=0.6)
        mpl.rcParams.update(_JOURNAL_RC)
    elif mode == 'screen':
        _STYLE.update(mode='screen', text=1.0, line=1.0)
        mpl.rcParams.update({k: mpl.rcParamsDefault[k] for k in _JOURNAL_RC})
    else:
        raise ValueError("mode must be 'screen' or 'journal'")


def sample_label(name):
    """Map a raw sample column ('bonemarrow3', 'd15_rep1', 'bm4', 'plasmidlib')
    to its display label ('Bone Marrow', 'D15', …). Unknown names pass through."""
    if name in SAMPLE_COLORS:
        return name
    low = str(name).lower()
    for prefix in sorted(_SAMPLE_PREFIX, key=len, reverse=True):
        if low.startswith(prefix):
            return _SAMPLE_PREFIX[prefix]
    return name


def sample_color(name):
    return SAMPLE_COLORS.get(sample_label(name), '#aaaaaa')


def shade(cmap, level):
    """One shade of a matplotlib colormap — the h2m way of choosing colours,
    e.g. shade('Oranges', .8) strong fill, shade('Blues', .2) pale bar fill."""
    return plt.get_cmap(cmap)(level)


def blend(c1, c2):
    """Average two colours (h2m: colour a human/mouse pair by both calls)."""
    return mcolors.to_hex((np.array(mcolors.to_rgb(c1)) + np.array(mcolors.to_rgb(c2))) / 2)


# ── Axes helpers ─────────────────────────────────────────────────────────────
def style_ax(ax, labelsize=16, hide=('top', 'right')):
    """Default look: open frame (no top/right spine), outward ticks. In journal
    mode tick labels are always 6 pt so they match across panels."""
    ax.spines[list(hide)].set_visible(False)
    ax.tick_params(axis='both', which='major',
                   labelsize=6 if _STYLE['mode'] == 'journal' else labelsize,
                   direction='out', length=ms(6), width=lw(1))
    return ax


def economist_ax(ax, grid='y'):
    """h2m bar/histogram look: light gridlines behind the data, no frame except a
    heavy butt-capped left spine. grid='x' (horizontal bars) also moves the value
    labels to the top axis."""
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=GRID_COLOR, lw=lw(1.2))
    ax.spines[['right', 'top', 'bottom']].set_visible(False)
    ax.spines['left'].set_lw(lw(1.5))
    ax.spines['left'].set_capstyle('butt')
    if grid == 'x':
        ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, length=0)
    return ax


def log_axis(ax, which='x'):
    """Log scale with 10^n tick labels and short minor ticks (Fig 2b)."""
    getattr(ax, f'set_{which}scale')('log')
    axis = ax.xaxis if which == 'x' else ax.yaxis
    axis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=10))
    axis.set_major_formatter(mpl.ticker.LogFormatterMathtext())
    ax.tick_params(axis=which, which='minor', length=ms(3), width=lw(0.8), direction='out')
    return ax


def legend_patches(ax, colors, labels=None, **kw):
    """Legend of black-edged colour swatches; `colors` is {label: colour}.
    Use instead of seaborn's auto legend (pass legend=False to seaborn)."""
    labels = list(colors) if labels is None else labels
    handles = [Patch(facecolor=colors[l], edgecolor='black', linewidth=lw(.8), label=l) for l in labels]
    kw = {'frameon': False, 'fontsize': fs(11), 'loc': 'best', **kw}
    return ax.legend(handles=handles, **kw)


def legend_top(ax, colors, labels=None, ncol=None, edge=False, marker=None, **kw):
    """Nature-style legend: one row of flat square swatches just above the plot,
    no frame (Nature 2024 La Fig b/d). Nat Biotech 2024 PEGG Fig b instead stacks
    black-edged swatches in the top-left corner — edge=True, ncol=1,
    loc='upper left', bbox_to_anchor=(0, 1).
    `colors` is {label: colour}; labels may hold math text, e.g. italic('DNMT1').
    marker='o' gives dot handles (dot plots)."""
    labels = list(colors) if labels is None else labels
    if marker:
        handles = [Line2D([], [], ls='', marker=marker, markersize=ms(9), markerfacecolor=colors[l],
                          markeredgewidth=0, label=l) for l in labels]
    else:
        handles = [Patch(facecolor=colors[l], edgecolor='black' if edge else 'none',
                         linewidth=lw(.8), label=l) for l in labels]
    kw = {'loc': 'lower left', 'bbox_to_anchor': (0, 1.01), 'ncol': ncol or len(labels),
          'frameon': False, 'fontsize': fs(14), 'handlelength': 1, 'handleheight': 1,
          'handletextpad': .4, 'columnspacing': 1.2, 'borderaxespad': 0, **kw}
    return ax.legend(handles=handles, **kw)


def panel_label(ax, letter, dx=0, dy=0):
    """Bold lowercase panel letter at the top-left of the panel's full extent
    (tick labels and axis label included), as in Nature figures. Call after
    fig.tight_layout(); dx/dy nudge it in figure fraction."""
    fig = ax.figure
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    return fig.text(bb.x0 + dx, bb.y1 + dy, letter, ha='right', va='bottom',
                    fontsize=8 if _STYLE['mode'] == 'journal' else 22,
                    fontfamily='Arial', fontweight='bold')


def panel_labels(rows, letters='abcdefghijklmnopqrstuvwxyz'):
    """Letter a whole figure, Nature style: letters run left-to-right, top-to-bottom
    and share one baseline per row (the highest panel top in that row), so a
    legend above one panel doesn't lift only its letter. `rows` is a 2-D list of
    axes, e.g. axs from plt.subplots or [[ax_a, ax_b], [ax_c, None]] (None skips).
    Call after fig.tight_layout()."""
    rows = [list(np.atleast_1d(r)) for r in np.atleast_1d(rows)] if isinstance(rows, np.ndarray) else rows
    first = next(a for r in rows for a in r if a is not None)
    fig = first.figure
    fig.canvas.draw()
    inv, rend = fig.transFigure.inverted(), fig.canvas.get_renderer()
    out, k = [], 0
    for r in rows:
        boxes = [(a, a.get_tightbbox(rend).transformed(inv)) for a in r if a is not None]
        top = max(bb.y1 for _, bb in boxes)
        for a, bb in boxes:
            out.append(fig.text(bb.x0, top, letters[k], ha='right', va='bottom',
                                fontsize=8 if _STYLE['mode'] == 'journal' else 22,
                                fontfamily='Arial', fontweight='bold'))
            k += 1
    return out


# ── Text helpers ─────────────────────────────────────────────────────────────
def italic(s):
    """Math-text italics for gene / protein names inside a label:
    f"{kp.italic('DNMT1')} +5 G-to-T"."""
    return r'$\mathit{' + str(s).replace(' ', r'\ ').replace('-', r'\text{-}') + '}$'


def corr_stats(x, y):
    """(Pearson r, Spearman rho)."""
    return scipy.stats.pearsonr(x, y)[0], scipy.stats.spearmanr(x, y)[0]


def corr_text(x, y, prefix='', style='nature'):
    """Both coefficients, always. style='nature' (default): 'ρ = 0.58 and R = 0.62'
    (Nat Biotech 2024 PEGG Fig c). style='kd': '$R_P$ = 0.62 | $R_S$ = 0.58'
    (the older B-ALL / h2m titles)."""
    r_p, r_s = corr_stats(x, y)
    if style == 'nature':
        return f'{prefix}ρ = {r_s:.2f} and $R$ = {r_p:.2f}'
    return f'{prefix}$R_P$ = {r_p:.2f} | $R_S$ = {r_s:.2f}'


def format_p(p):
    """Nature style: italic P, decimals down to 0.01, then '× 10^n':
    'P = 0.53', 'P = 6.7 × 10$^{-3}$'."""
    if p >= 0.01:
        return f'$P$ = {p:.2f}'
    mant, exp = f'{p:.1e}'.split('e')
    return f'$P$ = {mant} × 10$^{{{int(exp)}}}$'


def fold_p_text(fold, p=None):
    """'0.36×,\\nP = 6.7 × 10^-3' — the comparison label of Nature 2024 La Fig d."""
    return f'{fold:.2f}×' + (f',\n{format_p(p)}' if p is not None else '')


def mutation_label(row, gene='gene_name_m_corrected', mut='HGVSp_m', mut_h='HGVSp_h',
                   italic_gene=True):
    """'*Kras* G12D'-style label for a guide row — gene name in italics (Nature
    convention; italic_gene=False for plain text); falls back to the human
    protein change marked '(H)' when the mouse one is missing or legacy."""
    m = row.get(mut)
    human = pd.isna(m) or m in ('', 'Non-coding', 'Not calculated (legacy)')
    if human:
        m = row.get(mut_h)
    m = str(m)
    if m.startswith('p.'):
        m = m[2:]
    g = italic(row.get(gene)) if italic_gene else row.get(gene)
    return f"{g} {m}{' (H)' if human else ''}"


def label_points(ax, x, y, labels, colors='black', fontsize=10, halo=False, box=False,
                 avoid=None, **adjust_kw):
    """Label points with adjust_text: text coloured like its point, thin black
    leader lines. halo=True adds a white outline (busy backgrounds); box=True a
    translucent rounded box in the label colour (h2m 2g). `avoid` = (xs, ys) of
    points the labels must not cover — pass the highlighted points, not every
    background point (slow, and background points can be covered)."""
    from adjustText import adjust_text
    labels = list(labels)
    if isinstance(colors, str) or (isinstance(colors, tuple) and mcolors.is_color_like(colors)):
        colors = [colors] * len(labels)
    effects = [PathEffects.withStroke(linewidth=lw(2), foreground='w')] if halo else None
    texts = []
    for xi, yi, t, c in zip(x, y, labels, colors):
        kw = dict(fontsize=fs(fontsize), color=c, zorder=5, path_effects=effects)
        if box:
            kw['bbox'] = dict(facecolor='white', edgecolor=c, boxstyle='round,pad=0.3', alpha=.3)
        texts.append(ax.text(xi, yi, t, **kw))
    if texts:
        if avoid is not None:
            adjust_kw.setdefault('x', np.asarray(avoid[0], float))
            adjust_kw.setdefault('y', np.asarray(avoid[1], float))
        adjust_kw.setdefault('expand', (1.2, 1.4))
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='black', lw=lw(0.5)), **adjust_kw)
    return texts


def save(fig, path, dpi=300, transparent=False):
    """Save with house settings (tight bbox). PDF by default: vector, editable text."""
    folder = os.path.dirname(path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    fig.savefig(path, bbox_inches='tight', dpi=dpi, transparent=transparent)
    return path


# ── Brackets & grouped axes (Nature panels) ──────────────────────────────────
def group_bracket(ax, x0, x1, y, text, fontsize=14, tick=None):
    """Label a run of bars (Nature 2024 La Fig b, 'Parental K562 PEmax cells'):
    a line from x0 to x1 at height y (data units) with short down-turned ends,
    text centred above. Draw after the y-limits are final."""
    lo, hi = ax.get_ylim()
    tick = (hi - lo) * .03 if tick is None else tick
    ax.plot([x0, x0, x1, x1], [y - tick, y, y, y - tick], color='black', lw=lw(1), clip_on=False)
    return ax.text((x0 + x1) / 2, y + tick * .6, text, ha='center', va='bottom', fontsize=fs(fontsize))


def compare_bracket(ax, x1, x2, y, text, y1=None, y2=None, fontsize=13, gap=None):
    """Pairwise comparison bracket (Nature 2024 La Fig d): horizontal line at `y`
    joining two bars, legs dropping to just above each bar top (`y1`, `y2` — the
    tops bar_dots() returns), `text` centred above, e.g. fold_p_text(fc, p).
    Draw after the y-limits are final."""
    lo, hi = ax.get_ylim()
    gap = (hi - lo) * .03 if gap is None else gap
    y1 = y - gap * 2 if y1 is None else y1 + gap
    y2 = y - gap * 2 if y2 is None else y2 + gap
    ax.plot([x1, x1, x2, x2], [y1, y, y, y2], color='black', lw=lw(1), clip_on=False)
    return ax.text((x1 + x2) / 2, y + gap * .5, text, ha='center', va='bottom',
                   fontsize=fs(fontsize), linespacing=1.1)


def axis_groups(ax, groups, y=-0.35, fontsize=14, rotation=0, title=None):
    """Second tier of x labels under the tick labels (Nature 2024 La Fig d):
    for each (label, x0, x1) a line spanning x0..x1 and the label beneath it.
    `y` is in axes fraction (negative = below the axis; go lower when the tick
    labels are rotated). `title` prints left of the tier ('3′ end of\\npegRNA:')."""
    trans = mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for label, x0, x1 in groups:
        ax.plot([x0, x1], [y, y], color='black', lw=lw(1), transform=trans, clip_on=False)
        ax.text((x0 + x1) / 2, y - .02, label, transform=trans, ha='right' if rotation else 'center',
                va='top', fontsize=fs(fontsize), rotation=rotation,
                rotation_mode='anchor' if rotation else 'default')
    if title:
        ax.text(-.02, y - .02, title, transform=ax.transAxes, ha='right', va='top', fontsize=fs(fontsize))


# ── Plot recipes ─────────────────────────────────────────────────────────────
def bar_strip(data, x, y, order=None, palette=None, hue=None, hue_order=None,
              ax=None, log=False, point_size=7, errorbar=('ci', 95),
              point_alpha=.55, point_color=None):
    """Bar = group mean, every replicate overlaid as a semi-transparent dot.

    Whenever replicates exist, plot both: the bar carries the mean, the points
    carry n and spread. A bar alone hides both. For vertical grouped bars the
    DEFAULT is now bar_dots() (Nature look: flat bars, SD caps, bar positions for
    brackets); use bar_strip for horizontal bars, log axes, or the older B-ALL look.

    Points are semi-transparent (alpha .55) so overlapping replicates stay
    readable — lower to .3-.4 when n is large. They sit above the bar and keep a
    thin black edge so a single point still reads as a mark.

    Horizontal when `y` is the categorical column, vertical otherwise. Without
    `hue`, bars are coloured per category (defaults to SAMPLE_COLORS when every
    category is a known sample type). With `hue`, bars/dots are dodged and you
    add the legend yourself with legend_patches(). point_color='black' matches
    the published Nature panels; None keeps the category colour."""
    horizontal = not pd.api.types.is_numeric_dtype(data[y])
    cat = y if horizontal else x
    if order is None:
        present = list(pd.unique(data[cat]))
        order = [s for s in SAMPLE_ORDER if s in present] if set(present) <= set(SAMPLE_COLORS) else present
    if hue is None:
        hue, hue_order, dodge = cat, order, False
        if palette is None and set(order) <= set(SAMPLE_COLORS):
            palette = SAMPLE_COLORS
    else:
        dodge = True
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(*((5, 4) if horizontal else (max(4, 1.1 * len(order)), 4))))
    common = dict(data=data, x=x, y=y, hue=hue, order=order, hue_order=hue_order,
                  palette=palette, ax=ax, legend=False)
    sns.barplot(**common, edgecolor='black', linewidth=lw(1), errorbar=errorbar,
                zorder=1)
    strip = dict(common)
    if point_color is not None:          # one colour for every point
        strip.pop('palette'); strip['color'] = point_color
    sns.stripplot(**strip, dodge=dodge, edgecolor='black', linewidth=lw(.5),
                  size=ms(point_size), alpha=point_alpha, zorder=3)
    if log:
        log_axis(ax, 'x' if horizontal else 'y')
    (ax.set_ylabel if horizontal else ax.set_xlabel)('')
    style_ax(ax)
    return ax


def bar_dots(data, x, y, hue=None, order=None, hue_order=None, palette=None, ax=None,
             error='sd', edge=False, dot_color='black', dot_size=5, width=.8, rotation=None,
             seed=0):
    """Nature-style bars (Nature 2024 La Fig b/d; Nat Biotech 2024 PEGG Fig b):
    mean bar, black error bar with caps (error='sd' | 'sem' | None), every
    replicate as a small dot on its bar.

    edge=False — flat fill, black dots (La paper). edge=True with dot_color=None
    — black outlines and dots in the bar colour with a black rim (PEGG paper).
    `palette` is {hue level: colour} (or {category: colour} without hue), a list,
    or one colour. Tick labels rotate 45° automatically when they are long.

    Drawn with plain matplotlib so every bar position is known. Returns
    {(category, hue level): (x, top)} — `top` is the highest point of the bar,
    error bar and dots — ready for group_bracket() / compare_bracket()."""
    order = list(pd.unique(data[x])) if order is None else list(order)
    levels = [None] if hue is None else (list(pd.unique(data[hue])) if hue_order is None else list(hue_order))
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(max(4, .8 * len(order) * len(levels) ** .5), 4))
    keys = order if hue is None else levels
    if palette is None:
        palette = PALETTE_NATURE[0] if hue is None else categorical(len(levels), 'nature')
    if isinstance(palette, str) or (isinstance(palette, tuple) and mcolors.is_color_like(palette)):
        colors = {k: palette for k in keys}
    elif isinstance(palette, dict):
        colors = palette
    else:
        colors = dict(zip(keys, palette))

    rng = np.random.default_rng(seed)
    bw = width / len(levels)
    pos = {}
    for i, cat in enumerate(order):
        for j, lev in enumerate(levels):
            sel = data[x] == cat if lev is None else (data[x] == cat) & (data[hue] == lev)
            vals = data.loc[sel, y].dropna().to_numpy(float)
            if len(vals) == 0:
                continue
            xc = i - width / 2 + bw * (j + .5)
            c = colors[cat if lev is None else lev]
            mean = vals.mean()
            err = 0.0
            if error == 'sd' and len(vals) > 1:
                err = vals.std(ddof=1)
            elif error == 'sem' and len(vals) > 1:
                err = vals.std(ddof=1) / np.sqrt(len(vals))
            ax.bar(xc, mean, width=bw * .9, color=c, edgecolor='black' if edge else 'none',
                   linewidth=lw(1) if edge else 0, zorder=1)
            if err:
                ax.errorbar(xc, mean, yerr=err, fmt='none', ecolor='black', elinewidth=lw(1),
                            capsize=ms(5), capthick=lw(1), zorder=2)
            dots = xc + rng.uniform(-bw * .18, bw * .18, len(vals))
            if dot_color is None:
                ax.scatter(dots, vals, s=area(dot_size ** 2), color=c, edgecolors='black',
                           linewidths=lw(.6), zorder=3)
            else:
                ax.scatter(dots, vals, s=area(dot_size ** 2), color=dot_color, linewidths=0, zorder=3)
            pos[(cat, lev)] = (xc, max(mean + err, vals.max()))

    if rotation is None:
        rotation = 45 if max(len(str(c)) for c in order) > 6 else 0
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=rotation, ha='right' if rotation else 'center',
                       rotation_mode='anchor' if rotation else 'default')
    ax.set_xlim(-.6, len(order) - .4)
    ax.set_xlabel('')
    style_ax(ax)
    return pos


def dot_mean(data, cat, value, hue=None, order=None, hue_order=None, palette=None, ax=None,
             dot_size=7, spread=.25, legend=True):
    """Horizontal dot plot with a mean marker (Nature 2024 La Fig e): one row per
    category (first on top), one dot per measurement coloured by `hue` (e.g.
    target site), each hue level at a fixed small vertical offset, and a short
    black vertical line at the row mean. Legend of dots below the axis, 2 columns."""
    order = list(pd.unique(data[cat])) if order is None else list(order)
    levels = [None] if hue is None else (list(pd.unique(data[hue])) if hue_order is None else list(hue_order))
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(5, 1 + .9 * len(order)))
    if palette is None:
        palette = PALETTE_NATURE_SITES if len(levels) <= 4 else categorical(len(levels), 'nature')
    colors = palette if isinstance(palette, dict) else dict(zip(levels, palette))
    offsets = np.linspace(-spread / 2, spread / 2, len(levels)) if len(levels) > 1 else [0]
    ys = np.arange(len(order))[::-1]
    for yc, c in zip(ys, order):
        rows = data[data[cat] == c]
        for off, lev in zip(offsets, levels):
            vals = rows[value] if lev is None else rows.loc[rows[hue] == lev, value]
            ax.scatter(vals, np.full(len(vals), yc + off), s=area(dot_size ** 2),
                       color=colors[lev] if lev is not None else PALETTE_NATURE[0],
                       linewidths=0, zorder=3)
        m = rows[value].mean()
        ax.plot([m, m], [yc - spread, yc + spread], color='black', lw=lw(2), zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels(order)
    ax.set_ylim(-.6, len(order) - .4)
    style_ax(ax)
    if legend and hue is not None:
        legend_top(ax, colors, levels, marker='o', ncol=2, loc='upper center',
                   bbox_to_anchor=(.5, -.28))
    return ax


def box_strip(data, x, y, order=None, palette=None, pairs=None, test='t-test_ind',
              correction='bonferroni', ax=None, points='white', point_size=8, stars=False):
    """Box (no fliers) + every point, optional comparison brackets
    (h2m Ext Fig 4c/4d). points='white': coloured boxes, white dots;
    points='colored': white boxes, coloured dots. `pairs` = [(a, b), …].
    Brackets print the exact corrected P (Nature style, format_p) by default;
    stars=True gives the older */**/*** labels."""
    order = list(pd.unique(data[x])) if order is None else order
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(max(3, 1.3 * len(order)), 5))
    if points == 'white':
        sns.boxplot(data=data, x=x, y=y, hue=x, order=order, hue_order=order, palette=palette,
                    showfliers=False, legend=False, linewidth=lw(1.2), ax=ax)
        sns.stripplot(data=data, x=x, y=y, order=order, color='white', alpha=.5,
                      edgecolor='black', linewidth=lw(1), size=ms(point_size), ax=ax)
    else:
        sns.boxplot(data=data, x=x, y=y, order=order, color='white', showfliers=False,
                    linewidth=lw(1.2), ax=ax)
        sns.stripplot(data=data, x=x, y=y, hue=x, order=order, hue_order=order, palette=palette,
                      edgecolor='black', linewidth=lw(1), size=ms(point_size), alpha=.8, legend=False, ax=ax)
    ax.set_xlabel('')
    style_ax(ax)
    if pairs:
        from statannotations.Annotator import Annotator
        ann = Annotator(ax, pairs, data=data, x=x, y=y, order=order, verbose=False)
        ann.configure(test=test, text_format='star', loc='outside', fontsize=fs(16),
                      line_width=lw(1), comparisons_correction=correction, verbose=0)
        if stars:
            ann.apply_and_annotate()
        else:
            ann.apply_test()             # p-values come back already corrected
            ann.set_custom_annotations([format_p(a.data.pvalue) for a in ann.annotations])
            ann.annotate()
    return ax


def volcano(df, lfc_col, fdr_col, ax=None, label_col=None, highlight=None,
            fdr_cutoff=0.1, lfc_cutoff=1.0, n_label=10, fontsize=8):
    """Median LFC vs −log10(FDR) (B-ALL Fig 2f).

    Non-hits and un-highlighted hits are light grey; highlighted hits are red
    (enriched) / blue (depleted) with a black edge. `highlight` is a boolean mask
    aligned to df (e.g. gene in a leukemia list); default = every hit. With
    `label_col`, the top `n_label` enriched and depleted highlighted hits are labelled."""
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(5, 5))
    sub = df.dropna(subset=[lfc_col, fdr_col]).copy()
    sub['_y'] = -np.log10(sub[fdr_col].clip(lower=1e-300))
    hit = (sub[fdr_col] < fdr_cutoff) & (sub[lfc_col].abs() >= lfc_cutoff)
    if highlight is None:
        hl = hit
    else:
        if not isinstance(highlight, pd.Series):
            highlight = pd.Series(np.asarray(highlight), index=df.index)
        hl = hit & highlight.loc[sub.index].fillna(False).astype(bool)
    up, down = hl & (sub[lfc_col] > 0), hl & (sub[lfc_col] < 0)

    ax.scatter(sub.loc[~hit, lfc_col], sub.loc[~hit, '_y'], color='lightgrey', s=area(10), linewidths=0, zorder=1)
    ax.scatter(sub.loc[hit & ~hl, lfc_col], sub.loc[hit & ~hl, '_y'], color='lightgrey', s=area(15), linewidths=0, zorder=2)
    for m, c in [(up, DIRECTION_COLORS['enriched']), (down, DIRECTION_COLORS['depleted'])]:
        ax.scatter(sub.loc[m, lfc_col], sub.loc[m, '_y'], color=c,
                   edgecolors='black', linewidths=lw(0.5), s=area(50), zorder=4)

    ax.axhline(-np.log10(fdr_cutoff), color='grey', lw=lw(0.8), ls='--')
    ax.axvline(lfc_cutoff, color='grey', lw=lw(0.8), ls='--')
    ax.axvline(-lfc_cutoff, color='grey', lw=lw(0.8), ls='--')

    if label_col is not None:
        lab = pd.concat([sub[up].nlargest(n_label, lfc_col), sub[down].nsmallest(n_label, lfc_col)])
        colors = [DIRECTION_COLORS['enriched'] if v > 0 else DIRECTION_COLORS['depleted'] for v in lab[lfc_col]]
        label_points(ax, lab[lfc_col], lab['_y'], lab[label_col], colors=colors, fontsize=fontsize,
                     avoid=(sub.loc[hl, lfc_col], sub.loc[hl, '_y']))

    ax.set_xlabel('Median LFC (vs. Input)', fontsize=fs(14))
    ax.set_ylabel('-log10(FDR)', fontsize=fs(14))
    style_ax(ax, labelsize=12)
    return ax


def waterfall(df, lfc_col, fdr_col, ax=None, color_col=None, colors=None,
              fdr_cutoff=0.1, label_col=None, n_label=10, label_fontsize=9, legend=True):
    """Guides ranked by LFC (B-ALL Fig 2e). Significant guides are sized by FDR
    bin and coloured by `color_col` (default map: COSMIC role); non-significant
    guides are near-invisible. With `label_col`, the top/bottom `n_label`
    significant guides are listed as coloured text in the empty corners."""
    colors = COSMIC_COLORS if colors is None else colors
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(4, 8))
    d = df.dropna(subset=[lfc_col, fdr_col]).sort_values(lfc_col, ascending=False).copy()
    d['Rank'] = np.arange(1, len(d) + 1)
    d['_c'] = d[color_col].map(colors).fillna('tab:grey') if color_col else 'tab:grey'

    ns = d[d[fdr_col] >= fdr_cutoff]
    ax.scatter(ns['Rank'], ns[lfc_col], color=ns['_c'], s=area(1), alpha=.01)
    handles = []
    for i, cut in enumerate(FDR_BINS):
        upper = min(cut, fdr_cutoff)
        lower = FDR_BINS[i + 1] if i + 1 < len(FDR_BINS) else -np.inf
        m = (d[fdr_col] < upper) & (d[fdr_col] >= lower)
        ax.scatter(d.loc[m, 'Rank'], d.loc[m, lfc_col], color=d.loc[m, '_c'], s=area(FDR_SIZES[i]),
                   linewidth=lw(.2), edgecolor='black', alpha=.7)
        handles.append(Line2D([], [], ls='', marker='o', markersize=np.sqrt(area(FDR_SIZES[i])),
                              markerfacecolor='lightgrey', markeredgecolor='black', markeredgewidth=lw(.2),
                              label=f'FDR < $10^{{{-(i + 1)}}}$'))
    ax.axhline(0, linestyle='dashed', color='black', lw=lw(1))

    if label_col is not None:
        sig = d[d[fdr_col] < fdr_cutoff]
        top = sig[sig[lfc_col] > 0].head(n_label)
        bottom = sig[sig[lfc_col] < 0].tail(n_label)
        step = 0.032
        for j, (_, r) in enumerate(top.iterrows()):
            ax.text(0.30, 0.98 - j * step, f"{r['Rank']}. {r[label_col]}", color=r['_c'],
                    transform=ax.transAxes, fontsize=fs(label_fontsize), va='top')
        for j, (_, r) in enumerate(bottom.iterrows()):
            ax.text(0.03, 0.02 + (len(bottom) - 1 - j) * step, f"{r['Rank']}. {r[label_col]}", color=r['_c'],
                    transform=ax.transAxes, fontsize=fs(label_fontsize), va='bottom')
    if legend:
        ax.legend(handles=handles, loc='center left', frameon=False, fontsize=fs(10))
    ax.set_xlabel('gRNA Rank', fontsize=fs(16))
    ax.set_ylabel('Median LFC (vs. Input)', fontsize=fs(16))
    style_ax(ax)
    return ax


def corr_scatter(x, y, ax=None, color=None, ctrl=None, s=30, alpha=.5, edge=False, equal=True):
    """Paired measurements (screen A vs B, replicate vs replicate) with y = x and
    zero reference lines, titled with Pearson and Spearman (B-ALL Ext Fig 13c).
    ctrl: optional (x, y) of control guides, drawn in red with their own R line.
    edge=True gives the h2m look (black-edged, larger points)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(5, 5))
    color = GUIDE_COLORS['targeting'] if color is None else color
    ekw = dict(edgecolors='black', linewidths=lw(.5)) if edge else dict(linewidths=0)
    if ctrl is not None:
        ax.scatter(ctrl[0], ctrl[1], color=GUIDE_COLORS['control'], s=area(s), alpha=alpha, label='Controls', **ekw)
    ax.scatter(x, y, color=color, s=area(s), alpha=alpha, label='Targeting guides' if ctrl is not None else None, **ekw)

    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    if lo < 0 < hi:
        lo, hi = -max(-lo, hi) * 1.05, max(-lo, hi) * 1.05
        ax.axhline(0, color='grey', lw=lw(0.8), ls=':')
        ax.axvline(0, color='grey', lw=lw(0.8), ls=':')
    else:
        pad = (hi - lo) * .05
        lo, hi = lo - pad, hi + pad
    if equal:
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect('equal')
    ax.plot([lo, hi], [lo, hi], color='grey', ls='--', lw=lw(.8), zorder=0)

    title = corr_text(x, y, prefix='Targeting: ' if ctrl is not None else '')
    if ctrl is not None:
        title += '\n' + corr_text(ctrl[0], ctrl[1], prefix='Controls: ')
        ax.legend(fontsize=fs(10), frameon=False)
    ax.set_title(title, fontsize=fs(12))
    style_ax(ax, labelsize=13)
    return ax


def hexbin_density(x, y, ax=None, logx=False, logy=False, gridsize=40, cmap='inferno',
                   cbar_label='Points in bin', title=None, corr=True):
    """Density scatter for thousands of points (Nat Biotech 2024 PEGG Fig c):
    hexbin with log-scaled counts in inferno, a slim colour bar labelled on its
    right (rotated), log axes as needed. Title = `title` followed by
    'ρ = … and R = …' when corr=True."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    if logx:
        keep &= x > 0
    if logy:
        keep &= y > 0
    x, y = x[keep], y[keep]
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(5.5, 4.5))
    hb = ax.hexbin(x, y, gridsize=gridsize, cmap=cmap, bins='log', mincnt=1, linewidths=0,
                   xscale='log' if logx else 'linear', yscale='log' if logy else 'linear')
    cb = ax.figure.colorbar(hb, ax=ax, pad=.03, fraction=.05, aspect=25)
    cb.set_label(cbar_label, rotation=270, va='bottom', fontsize=fs(14))
    cb.ax.tick_params(labelsize=fs(13), width=lw(.8), length=ms(4))
    cb.outline.set_linewidth(lw(.8))
    if corr:
        text = corr_text(x, y, style='nature')
        title = f'{title} – {text}' if title else text
    if title:
        ax.set_title(title, fontsize=fs(14))
    style_ax(ax, labelsize=13)
    return ax


def z_heatmap(data, ax=None, vlim=5, fmt='.2f', annot_size=11, cmap='seismic', cbar=True):
    """Annotated diverging heatmap for LFC / z-scores (B-ALL Fig 3c-d, Fig 4
    barcoding): seismic, symmetric ±vlim, white cell borders."""
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(1.5 + 1.2 * data.shape[1], 0.8 + 0.5 * data.shape[0]))
    sns.heatmap(data, cmap=cmap, annot=True, fmt=fmt, annot_kws={'size': fs(annot_size)},
                linewidth=lw(1), linecolor='white', vmin=-vlim, vmax=vlim, cbar=cbar, ax=ax)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fs(11))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=fs(12))
    ax.set_xlabel('')
    ax.set_ylabel('')
    return ax


def exon_heatmap(mats, fig=None, cmap='Blues', vmin=0, vmax=None, cbar_label='',
                 exon_labels='bottom', strands=('+', '-'), cbar_ticks=None):
    """Per-exon small-multiple heatmaps (Nat Biotech 2024 PEGG Fig d — protospacer
    bias along a gene). `mats` = {exon: {'+': 2-D array, '-': 2-D array}}; rows are
    a within-site coordinate (e.g. position in protospacer), columns positions
    along the exon. Panel widths follow each exon's column count; every panel is a
    thin black box with faint white row lines; one colour bar on the right.
    Counts → 'Blues', efficiencies → 'Reds'. Stack two calls on
    fig.subfigures(2, 1) for the counts-over-efficiency layout of the paper.
    Returns the figure (or subfigure)."""
    if fig is None:
        fig = plt.figure(figsize=_figsize(14, 3))
    exons = list(mats)
    widths = [np.asarray(mats[e][strands[0]]).shape[1] for e in exons]
    vmax = vmax if vmax is not None else np.nanmax([np.nanmax(mats[e][s]) for e in exons for s in strands])
    gs = fig.add_gridspec(len(strands), len(exons) + 1, width_ratios=widths + [max(1, sum(widths) * .012)],
                          wspace=.06, hspace=.04)
    im = None
    for i, e in enumerate(exons):
        for r, s in enumerate(strands):
            m = np.asarray(mats[e][s], float)
            ax = fig.add_subplot(gs[r, i])
            im = ax.imshow(m, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto', interpolation='nearest')
            ax.set_yticks(np.arange(-.5, m.shape[0], 1), minor=True)
            ax.grid(which='minor', axis='y', color='white', lw=lw(.3))
            ax.tick_params(which='both', length=0, labelbottom=False, labelleft=False)
            for sp in ax.spines.values():
                sp.set_visible(True)
                sp.set_linewidth(lw(.8))
            if i == 0:
                ax.set_ylabel('−' if s == '-' else s, rotation=0, ha='right', va='center',
                              fontsize=fs(16), labelpad=ms(8))
            if exon_labels == 'bottom' and r == len(strands) - 1:
                ax.set_xlabel(str(e), fontsize=fs(15), labelpad=ms(4))
            if exon_labels == 'top' and r == 0:
                ax.set_title(str(e), fontsize=fs(15), pad=ms(4))
    cax = fig.add_subplot(gs[:, -1])
    cb = fig.colorbar(im, cax=cax, ticks=cbar_ticks if cbar_ticks is not None else [vmin, vmax])
    cb.set_label(cbar_label, rotation=270, va='bottom', fontsize=fs(14))
    cb.ax.tick_params(labelsize=fs(13), length=0)
    cb.outline.set_linewidth(lw(.8))
    return fig


def corr_matrix(data, method='spearman', figsize=None, annot_size=10, vmin=0, vmax=1):
    """Sample × sample correlation matrix (B-ALL Ext Fig 8, 14): Reds, annotated,
    original order (no clustering). Returns the seaborn ClusterGrid; save with g.savefig."""
    corr = data.corr(method=method)
    g = sns.clustermap(corr, cmap='Reds', vmin=vmin, vmax=vmax, annot=True, fmt='.2f',
                       linewidth=lw(.5), annot_kws={'size': fs(annot_size)},
                       row_cluster=False, col_cluster=False, figsize=figsize or _figsize(8, 8))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=fs(11), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=fs(11), rotation=0)
    return g


def pca_samples(counts, samples=None, ax=None, method='spearman', colors=None, legend=True):
    """PCA of the sample × sample correlation of guide counts (B-ALL Fig 2c),
    points coloured by tissue via sample_color()."""
    samples = list(counts.columns) if samples is None else list(samples)
    corr = counts[samples].corr(method=method).values
    centered = corr - corr.mean(axis=0)
    U, S, _ = np.linalg.svd(centered, full_matrices=False)
    pcs, var = U[:, :2] * S[:2], S ** 2 / (S ** 2).sum()
    colors = [sample_color(s) for s in samples] if colors is None else colors
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(6, 6))
    ax.scatter(pcs[:, 0], pcs[:, 1], c=colors, edgecolor='black', linewidths=lw(1), s=area(200), alpha=.8)
    if legend and colors is not None:
        present = [sample_label(s) for s in samples]
        legend_patches(ax, SAMPLE_COLORS, [l for l in SAMPLE_ORDER if l in present], fontsize=fs(11))
    ax.set_xlabel(f'PC 1\n{var[0] * 100:.1f}% explained variance', fontsize=fs(16))
    ax.set_ylabel(f'PC 2\n{var[1] * 100:.1f}% explained variance', fontsize=fs(16))
    style_ax(ax)
    return ax


def pale_barh(values, labels, ax=None, cmap='Blues', height=0.6, linewidth=1.5):
    """h2m horizontal bars (Fig 1h, 2g): pale fill and dark outline of one hue,
    light vertical grid, value labels on the top axis, first label on top."""
    if ax is None:
        _, ax = plt.subplots(figsize=_figsize(6, 0.5 + 0.5 * len(values)))
    pos = np.arange(len(values))[::-1]
    ax.barh(pos, values, height=height, color=shade(cmap, .2), edgecolor=shade(cmap, .9), linewidth=lw(linewidth))
    ax.set_yticks(pos, labels)
    economist_ax(ax, grid='x')
    ax.tick_params(axis='both', labelsize=fs(14))
    ax.tick_params(axis='y', length=0)
    return ax


def domain_diagram(ax, length, domains=(), y=0, label=None, mutations=(), numbers=(),
                   height=.5, fontsize=13, backbone=None):
    """Protein / construct schematic (Nature 2024 La Fig a): a thin grey backbone
    of `length` residues, domains as coloured blocks, names above, residue numbers
    below, mutations as red bars with red labels, construct name on the left.

    domains   [(start, end, name, colour), …] — name may be '' ; colours from
              SCHEMATIC_COLORS ('light', 'dark', 'tag') keep the paper's look
    mutations [(position, label), …]
    numbers   residue positions to print under the backbone
    Call once per construct with a decreasing `y` to stack them; turns the axis off."""
    backbone = SCHEMATIC_COLORS['backbone'] if backbone is None else backbone
    ax.add_patch(Rectangle((0, y - height * .22), length, height * .44, color=backbone, lw=0))
    for start, end, name, color in domains:
        ax.add_patch(Rectangle((start, y - height / 2), end - start, height, color=color, lw=0))
        if name:
            ax.text((start + end) / 2, y + height / 2 + height * .12, name, ha='center', va='bottom',
                    fontsize=fs(fontsize))
    for p, text in mutations:
        ax.plot([p, p], [y - height / 2, y + height / 2], color=SCHEMATIC_COLORS['mutation'],
                lw=lw(2.5), solid_capstyle='butt')
        ax.text(p, y + height / 2 + height * .12, text, color=SCHEMATIC_COLORS['mutation'],
                ha='center', va='bottom', fontsize=fs(fontsize))
    for p in numbers:
        ax.text(p, y - height / 2 - height * .12, str(p), ha='center', va='top', fontsize=fs(fontsize))
    if label:
        ax.text(-length * .03, y, label, ha='right', va='center', fontsize=fs(fontsize))
    ax.update_datalim([(0, y - height), (length, y + height)])
    ax.autoscale_view()
    ax.axis('off')
    return ax
