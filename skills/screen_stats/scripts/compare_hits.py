#!/usr/bin/env python3
"""Diff the hit calls of two <screen>_LFC_FDR.csv tables.

This is the replicate-exclusion sensitivity analysis: run the screen twice, once
with a suspect replicate in and once with it out, then ask which hits moved. It
is also how you check a Fisher-vs-Stouffer disagreement, and how you compare two
tables you are not sure were produced the same way.

No lifted logic -- the hit rule is the one printed by the engine's own summary
(FDR < 0.1 by the chosen combiner, |LFC| > 1 by the chosen summary).

    python compare_hits.py A_LFC_FDR.csv B_LFC_FDR.csv
    python compare_hits.py A.csv B.csv --conditions bm,spleen --fdr 0.05
    python compare_hits.py A.csv B.csv --method stouffers --lfc-summary avg
    python compare_hits.py A.csv B.csv --out hit_diff.csv

House default, matching the published figures: --lfc-summary median --method
fishers --fdr 0.1 --lfc 1. Control guides are excluded from the hit sets unless
--keep-controls, because BH correction in this pipeline runs over controls too.
"""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pandas as pd  # noqa: E402

import screenstats as S  # noqa: E402


def conditions_in(df, method):
    """Conditions with a collapsed FDR column.

    Must NOT match the per-tail columns: the table carries FDR_<c>_<method>
    alongside FDR_high_<c>_<method> and FDR_low_<c>_<method>, so a greedy
    '^FDR_(.+)_method$' would invent conditions named 'high_bm' and 'low_bm'.
    """
    pat = re.compile(r'^FDR_(?!high_|low_)(.+)_' + method + '$')
    out = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            out.append(m.group(1))
    return out


def hit_set(df, cond, method, summary, fdr, lfc, keep_controls):
    fdr_col = 'FDR_' + cond + '_' + method
    lfc_col = 'LFC_' + summary + '_' + cond
    if fdr_col not in df.columns:
        return None, None
    d = df
    if not keep_controls and 'classification' in d.columns:
        d = d[~d['classification'].isin(S.CONTROL_CLASSES)]
    sig = d[fdr_col] < fdr
    if lfc_col in d.columns:
        # NaN in LFC_concordant_mean_* (ambiguous-direction guides) compares
        # False here, so those guides silently never become hits.
        up = d[sig & (d[lfc_col] > lfc)]
        down = d[sig & (d[lfc_col] < -lfc)]
    else:
        up, down = d[sig], d.iloc[0:0]
    return set(up['gRNA_id']), set(down['gRNA_id'])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('table_a')
    ap.add_argument('table_b')
    ap.add_argument('--conditions', help='comma-separated; default = all shared')
    ap.add_argument('--method', default='fishers', choices=['fishers', 'stouffers'])
    ap.add_argument('--lfc-summary', default='median',
                    choices=['median', 'avg', 'concordant_mean'])
    ap.add_argument('--fdr', type=float, default=0.1)
    ap.add_argument('--lfc', type=float, default=1.0)
    ap.add_argument('--keep-controls', action='store_true')
    ap.add_argument('--out', help='write the per-guide diff to this CSV')
    args = ap.parse_args(argv)

    a = pd.read_csv(args.table_a)
    b = pd.read_csv(args.table_b)
    na, nb = os.path.basename(args.table_a), os.path.basename(args.table_b)

    print(na + ': ' + str(a.shape[0]) + ' guides x ' + str(a.shape[1]) + ' columns')
    print(nb + ': ' + str(b.shape[0]) + ' guides x ' + str(b.shape[1]) + ' columns')

    # replicate membership differs -> the two tables are not the same analysis
    reps_a = set(a.columns) - set(b.columns)
    reps_b = set(b.columns) - set(a.columns)
    if reps_a or reps_b:
        print('')
        print('column sets differ:')
        if reps_a:
            print('  only in ' + na + ': ' + repr(sorted(reps_a)[:20]))
        if reps_b:
            print('  only in ' + nb + ': ' + repr(sorted(reps_b)[:20]))
        print('  If these are replicate columns, the two tables include different')
        print('  samples -- that is the comparison, not a bug. If you did not mean')
        print('  to compare different replicate sets, check which table you loaded.')

    ca = conditions_in(a, args.method)
    cb = conditions_in(b, args.method)
    conds = args.conditions.split(',') if args.conditions else [c for c in ca if c in cb]
    only_a = [c for c in ca if c not in cb]
    only_b = [c for c in cb if c not in ca]
    if only_a or only_b:
        print('conditions tested in only one table: ' + na + repr(only_a)
              + ' ' + nb + repr(only_b))

    rows = []
    print('')
    print('hit rule: FDR_' + args.method + ' < ' + str(args.fdr) + ' and |LFC_'
          + args.lfc_summary + '| > ' + str(args.lfc)
          + ('' if args.keep_controls else ', controls excluded'))
    for cond in conds:
        ha = hit_set(a, cond, args.method, args.lfc_summary, args.fdr, args.lfc,
                     args.keep_controls)
        hb = hit_set(b, cond, args.method, args.lfc_summary, args.fdr, args.lfc,
                     args.keep_controls)
        if ha[0] is None or hb[0] is None:
            print('  ' + cond + ': not tested in both tables, skipped')
            continue
        print('')
        print('  ' + cond + ':')
        for direction, sa, sb in [('enriched', ha[0], hb[0]),
                                  ('depleted', ha[1], hb[1])]:
            retained = sa & sb
            lost = sa - sb
            gained = sb - sa
            print('    ' + direction + ': ' + str(len(sa)) + ' in ' + na + ', '
                  + str(len(sb)) + ' in ' + nb + ' -> ' + str(len(retained))
                  + ' retained, ' + str(len(lost)) + ' lost, '
                  + str(len(gained)) + ' gained')
            if lost:
                print('      lost:   ' + ', '.join(sorted(lost)[:8])
                      + ('...' if len(lost) > 8 else ''))
            if gained:
                print('      gained: ' + ', '.join(sorted(gained)[:8])
                      + ('...' if len(gained) > 8 else ''))
            for g in sorted(retained):
                rows.append((cond, direction, g, 'retained'))
            for g in sorted(lost):
                rows.append((cond, direction, g, 'lost'))
            for g in sorted(gained):
                rows.append((cond, direction, g, 'gained'))

        # Fisher vs Stouffer agreement within each table
        for lbl, df in [(na, a), (nb, b)]:
            f_col = 'FDR_' + cond + '_fishers'
            s_col = 'FDR_' + cond + '_stouffers'
            if f_col in df.columns and s_col in df.columns:
                d = df if args.keep_controls or 'classification' not in df.columns \
                    else df[~df['classification'].isin(S.CONTROL_CLASSES)]
                f = set(d[d[f_col] < args.fdr]['gRNA_id'])
                st = set(d[d[s_col] < args.fdr]['gRNA_id'])
                print('    [' + lbl + '] Fisher ' + str(len(f)) + ', Stouffer '
                      + str(len(st)) + ', overlap ' + str(len(f & st))
                      + ' (Fisher-only ' + str(len(f - st)) + ' = hits carried by '
                      'one extreme replicate)')

    if args.out and rows:
        out = pd.DataFrame(rows, columns=['condition', 'direction', 'gRNA_id', 'status'])
        out.to_csv(args.out, index=False)
        print('')
        print('wrote ' + args.out + ' (' + str(len(out)) + ' rows)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
