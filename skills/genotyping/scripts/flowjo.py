"""FlowJo gate-statistics .xls -> tidy per-sample table, percent-of-ancestor
chaining, absolute counts, fold expansion, and the three composition panels.

Lifted from
  /Users/kexindong/Documents/GitHub/hsc-ch-mutation/20260508-WT-HSC-LSK-culture/flow_plot.ipynb
  cell 0  -- xlrd parse, the nested `get()` gate lookup, the Q2 x Lin- chaining,
             the absolute-count math and the grouped bar panel
  cell 2  -- linear stacked bar
  cell 4  -- log-space stacked bar  (NOT cell 1: cell 1 is an earlier draft of
             the same figure that stacks Lin+ at the bottom and writes the same
             flow_stacked_log.pdf, so whichever ran last won. Cell 4 stacks LSK
             at the bottom, matching the linear panel, and is the intended one.)
Do not edit the arithmetic here without editing there -- these must not drift.

Two changes to the lifted code, both correctness, neither arithmetic:
  * `gate()` raises on a gate path that does not match, instead of silently
    returning (nan, 0) and producing a plot of nothing;
  * the hand-entered hemocytometer numbers (initial_LSK = 24_000,
    final_total_live = 0.452*3*12*1e6) are arguments, not literals.

Styling is NOT defined here. Figure style belongs to the sci-plot skill: this
module sets up via kdplot and only decides what goes on which axis.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

# Figure style is sci-plot's job, not this skill's.
_SCIPLOT = os.path.expanduser('~/.claude/skills/sci-plot/scripts')
if os.path.isdir(_SCIPLOT) and _SCIPLOT not in sys.path:
    sys.path.insert(0, _SCIPLOT)

# Default segment colors, carried over from the source notebook so existing
# panels still match. Override with the `colors` argument.
C_LSK = "#762a83"
C_LINNEG = "#4dac26"
C_LINPOS = "tab:grey"


def _kdplot():
    try:
        import kdplot
        return kdplot
    except ImportError:
        return None


# ---------------------------------------------------------------------------
# Parsing -- lifted from flow_plot.ipynb cell 0
# ---------------------------------------------------------------------------

def read_flowjo_xls(path):
    """Rows of a FlowJo 'export group statistics' workbook: [depth, name, stat, ncells].

    xlrd 2.x reads ONLY legacy BIFF .xls. FlowJo's export must be saved as .xls;
    an .xlsx export raises XLRDError here and needs openpyxl instead.
    """
    if not os.path.exists(path):
        raise SystemExit("No such flow statistics file: " + str(path))
    if path.lower().endswith('.xlsx'):
        raise SystemExit(
            "xlrd 2.x cannot read .xlsx, only legacy .xls. Either re-export from "
            "FlowJo as .xls, or read it with pandas.read_excel(engine='openpyxl') "
            "and pass the rows in directly."
        )
    import xlrd
    wb = xlrd.open_workbook(path)
    sh = wb.sheet_by_index(0)
    rows = [sh.row_values(r) for r in range(1, sh.nrows)]
    return rows


def list_gates(rows):
    """Every gate path in the workbook, in file order.

    Run this before naming any gate. FlowJo's quadrant names contain exact
    spacing that nobody reproduces from memory -- 'Q2: APC-A+ , BV421-A+' has a
    space before the comma.
    """
    return [str(name).strip() for _depth, name, _stat, _n in rows]


def list_specimens(rows):
    """Distinct Specimen_*.fcs names in the workbook.

    The export contains single-stain compensation controls as separate
    specimens (Specimen_001_APC.fcs, _BV421.fcs, _BV510.fcs, _PE.fcs,
    _PEcy7.fcs) alongside the real samples (Specimen_001_Sample 1/2/3.fcs).
    Never glob all specimens: the controls have near-zero LSK and would
    silently drag down every mean. This is why the notebook hardcoded
    samples = ['Sample 1', 'Sample 2', 'Sample 3'].
    """
    seen = []
    for path in list_gates(rows):
        head = path.split('/')[0]
        if head and head not in seen:
            seen.append(head)
    return seen


def gate(rows, full_path, required=True):
    """(statistic, n_cells) for one full slash-delimited gate path.

    Matching is exact string equality after .strip(), as in the notebook. The
    notebook returned (nan, 0) on a miss, which turned a typo into NaN columns
    and an empty plot; this raises instead, and shows the near-misses.
    """
    want = full_path.strip()
    for _depth, name, stat, ncells in rows:
        if str(name).strip() == want:
            return float(stat), int(ncells)
    if not required:
        return float('nan'), 0
    near = [g for g in list_gates(rows) if g.split('/')[-1] == want.split('/')[-1]]
    msg = "Gate path not found:\n  " + want
    if near:
        msg += "\nPaths with the same leaf name:\n  " + "\n  ".join(near[:8])
    msg += "\nRun list_gates() to see them all; FlowJo quadrant names have exact spacing."
    raise SystemExit(msg)


def pct_of_ancestor(*parent_percents):
    """Chain parent-relative percentages into a percent of a chosen ancestor.

    FlowJo statistics are percent-of-PARENT, not percent-of-total. The
    notebook's `lsk_of_sc = (q2_pct/100)*(lsk_pct/100)*100` exists for exactly
    this reason, and skipping the chaining is the commonest flow readout error.
    Anything two or more gates deep must be chained through every intermediate.

        pct_of_ancestor(9.83, 69.0)  ->  6.7827
    """
    out = 1.0
    for p in parent_percents:
        out = out * (float(p) / 100.0)
    return out * 100.0


def summarize_lsk(rows, samples, specimen_prefix='Specimen_001_'):
    """The LSK panel from cell 0, as a DataFrame; writes nothing.

    Gate tree assumed (as exported by the 08-May-2026 workspace):
        <specimen>/Lymphocytes/Single Cells
        .../Lin-
        .../Q2: APC-A+ , BV421-A+
        .../Q2: APC-A+ , BV421-A+/Lin-
        .../Q2: APC-A+ , BV421-A+/Lin-/CD150+
        .../Q2: APC-A+ , BV421-A+/Lin-/Q7: BV510-A+ , PE-Cy7-A-
    A different panel needs its own gate spec -- call `gate()` directly rather
    than bending this one.

    Every % is carried alongside its #Cells on purpose. Deep gates here have
    tiny denominators: Sample 1's LSK CD150+ CD34- gate is 4 cells out of 260
    LSK, and a percentage off single-digit counts is noise.
    """
    records = []
    for sname in samples:
        base = specimen_prefix + sname + '.fcs/Lymphocytes'
        sc = base + '/Single Cells'
        q2 = sc + '/Q2: APC-A+ , BV421-A+'

        sc_pct, sc_n = gate(rows, sc)
        lin_pct, lin_n = gate(rows, sc + '/Lin-')
        q2_pct, q2_n = gate(rows, q2)
        lsk_pct, lsk_n = gate(rows, q2 + '/Lin-')
        cd150_pct, cd150_n = gate(rows, q2 + '/Lin-/CD150+')
        cd34neg_pct, cd34neg_n = gate(rows, q2 + '/Lin-/Q7: BV510-A+ , PE-Cy7-A-')

        # True LSK fraction of Single Cells = Q2% x Lin-%
        lsk_of_sc = pct_of_ancestor(q2_pct, lsk_pct)

        records.append({
            'Sample': sname,
            'Single Cells %': sc_pct,
            'Single Cells #': sc_n,
            'Lin- %': lin_pct,
            'Lin- #': lin_n,
            'Q2 (APC+ BV421+) %': q2_pct,
            'Q2 (APC+ BV421+) #': q2_n,
            'LSK (Q2xLin-) %': lsk_of_sc,
            'LSK (Q2xLin-) #': lsk_n,
            'LSK CD150+ %': cd150_pct,
            'LSK CD150+ #': cd150_n,
            'LSK CD150+ CD34- %': cd34neg_pct,
            'LSK CD150+ CD34- #': cd34neg_n,
        })
    return pd.DataFrame(records)


def absolute_counts(df, initial_lsk, final_total_live):
    """Absolute counts and fold expansion -- cell 0 section 2 and cell 3.

    `final_total_live` is viability x volume(mL) x cells/mL; in the source that
    was 0.452 * 3 * 12e6, hand-entered from the counter. Both it and
    `initial_lsk` are arguments here, not literals.

    Means are taken across samples, and n / SD are returned alongside so the
    panel can show replicate spread -- the notebook only ever took .mean().
    """
    frac = df['LSK (Q2xLin-) %'] / 100.0
    cd150 = df['LSK CD150+ %'] / 100.0
    cd34neg = df['LSK CD150+ CD34- %'] / 100.0

    mean_lsk_frac = frac.mean()
    mean_cd150_frac = cd150.mean()
    mean_cd34neg_frac = cd34neg.mean()

    final_lsk = final_total_live * mean_lsk_frac
    final_lsk_cd150 = final_lsk * mean_cd150_frac
    final_lsk_cd34neg = final_lsk_cd150 * mean_cd34neg_frac

    lin_neg_mean = df['Lin- %'].mean() / 100.0
    return {
        'initial_LSK': initial_lsk,
        'final_total_live': final_total_live,
        'n_samples': int(len(df)),
        'mean_lsk_frac': mean_lsk_frac,
        'sd_lsk_frac': float(frac.std(ddof=1)) if len(df) > 1 else float('nan'),
        'mean_cd150_frac': mean_cd150_frac,
        'mean_cd34neg_frac': mean_cd34neg_frac,
        'final_LSK': final_lsk,
        'final_LSK_CD150': final_lsk_cd150,
        'final_LSK_CD34neg': final_lsk_cd34neg,
        'fold_expansion': final_lsk / initial_lsk if initial_lsk else float('nan'),
        'per_sample_LSK_counts': [final_total_live * f for f in frac],
        'lsk_frac': mean_lsk_frac,
        'lin_neg_nonlsk_frac': lin_neg_mean - mean_lsk_frac,
        'lin_pos_frac': 1.0 - lin_neg_mean,
    }


# ---------------------------------------------------------------------------
# Panels. Mechanics only -- fonts, sizes and palette are sci-plot's.
# ---------------------------------------------------------------------------

def _fmt_count(v):
    """K/M tick formatter, from cell 0."""
    if v >= 1e6:
        return "%.1fM" % (v / 1e6)
    if v >= 1e3:
        return "%.0fK" % (v / 1e3)
    return "%.0f" % v


def bar_counts(labels, heights, ax=None, colors=None, dots=None):
    """Grouped bar of absolute counts -- cell 0 section 3.

    `dots` optionally overlays per-sample values on one bar index as
    {index: [values]}, which the lab's sci-plot house style asks for and the
    notebook omitted.
    """
    import matplotlib.pyplot as plt
    kd = _kdplot()
    if ax is None:
        _fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(range(len(heights)), heights, edgecolor='black', linewidth=1,
           color=colors, alpha=0.8)
    if dots:
        for i, vals in dots.items():
            ax.plot([i] * len(vals), vals, 'o', color='black', ms=4, zorder=5)
    tick_vals = np.linspace(0, max(heights), 6)
    ax.set_yticks(tick_vals)
    ax.set_yticklabels([_fmt_count(v) for v in tick_vals])
    ax.set_ylabel('Cell Count')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    if kd:
        kd.style_ax(ax)
    else:
        ax.spines[['top', 'right']].set_visible(False)
    return ax


def stacked_linear(initial, total, segments, ax=None, labels=('Initial\nLSK', 'Final\nTotal Live')):
    """Linear stacked composition bar -- cell 2.

    `segments` is [(count, color, legend_label), ...] bottom-up; the source
    order is LSK, Lin-(non-LSK), Lin+, so LSK sits at the bottom.
    In-segment percentage labels are drawn for segments above 1%.
    """
    import matplotlib.pyplot as plt
    kd = _kdplot()
    if ax is None:
        _fig, ax = plt.subplots(figsize=(3, 4))

    ax.bar(0, initial, color=segments[0][1], edgecolor='black', linewidth=1, alpha=0.85)
    ax.text(0, initial * 1.02, "100%", ha='center', va='bottom')

    bot = 0
    for count, color, _lab in segments:
        ax.bar(1, count, bottom=bot, color=color, edgecolor='black',
               linewidth=0.5, alpha=0.85)
        pct = count / total * 100
        if pct > 1:
            ax.text(1, bot + count / 2, "%.1f%%" % pct, ha='center', va='center',
                    color='white', fontweight='bold')
        bot += count

    ax.set_xticks([0, 1])
    ax.set_xticklabels(list(labels))
    ax.set_ylabel('Cell Count')
    ax.set_ylim(0, total * 1.1)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: _fmt_count(v)))
    if kd:
        kd.style_ax(ax)
        # kdplot.legend_patches takes {label: colour}.
        kd.legend_patches(ax, dict((s[2], s[1]) for s in segments), loc='upper left')
    else:
        ax.spines[['top', 'right']].set_visible(False)
    return ax


def stacked_log(initial, total, segments, ax=None,
                labels=('Initial\nLSK', 'Final\nTotal Live')):
    """Log-space stacked composition bar -- cell 4 (not cell 1).

    IMPORTANT, and worth saying in the caption: this is a presentational
    construction, not a log axis. Each segment's height is
    log10(total) * (count/total), so segment heights are proportional to
    composition while the bar's total height reads as log10(count). Segment
    boundaries do NOT correspond to their own counts on the 10^n tick labels.
    Keep it paired with the linear panel or a reader will misread segment sizes.
    """
    import matplotlib.pyplot as plt
    kd = _kdplot()
    if ax is None:
        _fig, ax = plt.subplots(figsize=(5, 6))

    ax.bar(0, np.log10(initial), color=segments[0][1], edgecolor='black',
           linewidth=1, alpha=0.85)
    ax.text(0, np.log10(initial) + 0.05, "100%", ha='center', va='bottom')

    total_log = np.log10(total)
    bot = 0
    for count, color, _lab in segments:
        h = total_log * (count / total)
        ax.bar(1, h, bottom=bot, color=color, edgecolor='black',
               linewidth=0.5, alpha=0.85)
        pct = count / total * 100
        if pct > 1:
            ax.text(1, bot + h / 2, "%.1f%%" % pct, ha='center', va='center',
                    color='white', fontweight='bold')
        bot += h

    ax.set_xticks([0, 1])
    ax.set_xticklabels(list(labels))
    ax.set_ylabel('log10(Cell Count)')
    ymax = total_log * 1.15
    ax.set_ylim(0, ymax)
    yticks = [y for y in [3, 4, 5, 6, 7] if y <= ymax]
    ax.set_yticks(yticks)
    ax.set_yticklabels(["10$^" + str(y) + "$" for y in yticks])
    if kd:
        kd.style_ax(ax)
        kd.legend_patches(ax, dict((s[2], s[1]) for s in segments), loc='upper left')
    else:
        ax.spines[['top', 'right']].set_visible(False)
    return ax


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('xls', help='FlowJo "export group statistics" .xls')
    p.add_argument('--samples', nargs='+',
                   help='sample specimen names, e.g. "Sample 1" "Sample 2". '
                        'Omit with --list to see what is in the file.')
    p.add_argument('--specimen-prefix', default='Specimen_001_')
    p.add_argument('--list', action='store_true',
                   help='print specimens and gate paths, then exit')
    p.add_argument('--initial-lsk', type=float,
                   help='input cell number (hemocytometer), e.g. 24000')
    p.add_argument('--viability', type=float, help='live fraction, e.g. 0.452')
    p.add_argument('--volume-ml', type=float, help='harvest volume in mL, e.g. 3')
    p.add_argument('--conc-per-ml', type=float,
                   help='cells per mL, e.g. 12e6')
    p.add_argument('--no-plots', action='store_true')
    p.add_argument('-o', '--outdir', default='flow_out')
    args = p.parse_args()

    rows = read_flowjo_xls(args.xls)

    if args.list or not args.samples:
        print('Specimens in this workbook:')
        for s in list_specimens(rows):
            print('  ' + s)
        print('\nGate paths:')
        for g in list_gates(rows):
            print('  ' + g)
        print('\nPick only the real sample specimens with --samples; the '
              'single-stain compensation controls are listed above too.')
        return

    df = summarize_lsk(rows, args.samples, args.specimen_prefix)
    os.makedirs(args.outdir, exist_ok=True)
    summary_path = os.path.join(args.outdir, 'flow_summary.csv')
    df.to_csv(summary_path, index=False)
    print(df.to_string(index=False))
    print('\nwrote ' + summary_path)

    tiny = df[df['LSK CD150+ CD34- #'] < 20]
    if len(tiny):
        print('\nWARNING: these samples have <20 cells in the deepest gate, so '
              'that percentage is noise:')
        print(tiny[['Sample', 'LSK CD150+ CD34- #', 'LSK CD150+ CD34- %']].to_string(index=False))

    if args.initial_lsk is None:
        return
    if None in (args.viability, args.volume_ml, args.conc_per_ml):
        raise SystemExit(
            "Absolute counts need --viability, --volume-ml and --conc-per-ml "
            "together with --initial-lsk (final total live = viability x volume x conc)."
        )
    final_total_live = args.viability * args.volume_ml * args.conc_per_ml
    ac = absolute_counts(df, args.initial_lsk, final_total_live)

    print('')
    print('Final total live:   %.2fM' % (ac['final_total_live'] / 1e6))
    print('Final LSK:          {:,.0f}'.format(ac['final_LSK']))
    print('Final LSK CD150+:   {:,.0f}'.format(ac['final_LSK_CD150']))
    print('Final LSK CD34-:    {:,.0f}'.format(ac['final_LSK_CD34neg']))
    print('LSK fold expansion: %.1fx  (n=%d)' % (ac['fold_expansion'], ac['n_samples']))

    if args.no_plots:
        return

    kd = _kdplot()
    if kd:
        kd.setup()
    import matplotlib.pyplot as plt

    total = ac['final_total_live']
    segments = [
        (total * ac['lsk_frac'], C_LSK, 'LSK'),
        (total * ac['lin_neg_nonlsk_frac'], C_LINNEG, 'Lin- (non-LSK)'),
        (total * ac['lin_pos_frac'], C_LINPOS, 'Lin+'),
    ]

    heights = [ac['initial_LSK'], total, ac['final_LSK'],
               ac['final_LSK_CD150'], ac['final_LSK_CD34neg']]
    bar_labels = ['Initial\nLSK', 'Final\nTotal Live', 'Final\nLSK',
                  'Final\nLSK\nCD150+', 'Final\nLSK\nCD150+ CD34-']
    ax = bar_counts(bar_labels, heights,
                    colors=['tab:blue', 'tab:gray', 'tab:purple', 'tab:orange', 'tab:green'],
                    dots={2: ac['per_sample_LSK_counts']})
    _save(ax, os.path.join(args.outdir, 'flow_bar_plot.pdf'), kd)

    ax = stacked_linear(ac['initial_LSK'], total, segments)
    _save(ax, os.path.join(args.outdir, 'flow_stacked_linear.pdf'), kd)

    ax = stacked_log(ac['initial_LSK'], total, segments)
    _save(ax, os.path.join(args.outdir, 'flow_stacked_log.pdf'), kd)

    print('\nwrote flow_bar_plot.pdf, flow_stacked_linear.pdf, flow_stacked_log.pdf in '
          + args.outdir)
    print('These are QC panels. For a paper figure, restyle through the sci-plot skill.')
    plt.close('all')


def _save(ax, path, kd):
    import matplotlib.pyplot as plt
    fig = ax.figure
    fig.tight_layout()
    if kd:
        kd.save(fig, path)
    else:
        fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
