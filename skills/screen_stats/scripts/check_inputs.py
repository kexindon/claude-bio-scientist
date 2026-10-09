#!/usr/bin/env python3
"""Preflight a screen config before running the pipeline, and check the output
contract afterwards.

No lifted logic here -- this is a wrapper that asserts, before a run, the
conditions the notebook assumed silently and the ones the gotchas in SKILL.md
describe:

  * the count matrix's first two columns are exactly sgRNA, gene
    (rpm() treats df.columns[2:] as samples, so a stray index column from a
    to_csv/read_csv round trip silently corrupts every RPM value)
  * every replicate name in samps exists in the count matrix
  * every conditions_null entry is a key of samps
  * baseline mode is unambiguous (not both baseline_samps and counts_bc)
  * classification contains the exact control strings
  * how many control guides and pooled control LFCs each tested condition gets
  * how many guides the library and MLE inner joins would drop
  * single-replicate conditions in conditions_null

Usage:
    python check_inputs.py --config screens.yaml
    python check_inputs.py --config screens.yaml --screen ABE_OG
    python check_inputs.py --check-downstream LFC-FDR/ABE_OG_LFC_FDR.csv \
        --notebooks figure2_CLEANED.ipynb volcano_plots.ipynb

The --check-downstream mode greps figure notebooks for the column names they
reference and reports which ones the new table does not supply -- the actual
failure mode behind "renaming a condition breaks every figure notebook".

Exit status 1 if any ERROR was found; warnings alone exit 0.
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pandas as pd  # noqa: E402

import screenstats as S  # noqa: E402
from run_screen import load_config, read_counts  # noqa: E402

ERRORS = []
WARNINGS = []


def err(msg):
    ERRORS.append(msg)
    print('ERROR   ' + msg)


def warn(msg):
    WARNINGS.append(msg)
    print('WARN    ' + msg)


def ok(msg):
    print('ok      ' + msg)


# ---------------------------------------------------------------------------

def check_screen(name, spec):
    print('')
    print('=== ' + name + ' ===')

    counts_path = spec.get('counts')
    if not counts_path or not os.path.exists(counts_path):
        err(name + ': counts file not found: ' + repr(counts_path))
        return
    try:
        df = read_counts(counts_path)
    except SystemExit as e:
        err(name + ': ' + str(e))
        return
    ok('counts ' + counts_path + ': ' + str(df.shape[0]) + ' guides, '
       + str(df.shape[1] - 2) + ' sample columns, leading columns sgRNA/gene')

    samps = spec.get('samps') or {}
    if not samps:
        err(name + ': no samps block')
        return
    known = set(df.columns)
    for cond in samps:
        missing = [r for r in samps[cond] if r not in known]
        if missing:
            err(name + '/' + cond + ': replicate columns not in count matrix: '
                + repr(missing))
    unused = sorted(known - set(['sgRNA', 'gene'])
                    - set(r for g in samps.values() for r in g)
                    - set(spec.get('baseline_samps') or []))
    if unused:
        warn(name + ': sample columns present but never used: ' + repr(unused))

    conditions_null = spec.get('conditions_null') or []
    bad = [c for c in conditions_null if c not in samps]
    if bad:
        err(name + ': conditions_null entries missing from samps: ' + repr(bad))
    single = [c for c in conditions_null if c in samps and len(samps[c]) == 1]
    if single:
        warn(name + ': single-replicate conditions in conditions_null: ' + repr(single)
             + ' -- Fisher and Stouffer pass a lone p-value through unchanged, so '
             'the "combined" p carries no replicate-consistency evidence.')

    # baseline mode
    if spec.get('counts_bc'):
        if spec.get('baseline_samps'):
            err(name + ': both counts_bc and baseline_samps given -- ambiguous. '
                'The original notebook silently ignored baseline_samps here.')
        if not spec.get('samps_bc'):
            err(name + ': counts_bc given but no samps_bc')
        else:
            missing_c = [c for c in samps if c not in spec['samps_bc']]
            if missing_c:
                err(name + ': samps_bc is missing conditions present in samps: '
                    + repr(missing_c))
            if os.path.exists(spec['counts_bc']):
                dbc = read_counts(spec['counts_bc'])
                if len(dbc) != len(df):
                    err(name + ': barcoding count table has ' + str(len(dbc))
                        + ' rows vs ' + str(len(df)) + ' in the screen table; the '
                        'barcoding baseline is taken positionally (np.median over '
                        'rows), so the two tables must be the same guides in the '
                        'same order.')
                elif not df['sgRNA'].equals(dbc['sgRNA']):
                    err(name + ': barcoding table sgRNA order differs from the screen '
                        'table. The barcoding baseline is positional, not merged -- '
                        'mismatched order silently pairs the wrong guides.')
                else:
                    ok('barcoding baseline table aligns row-for-row')
        ok('baseline mode: barcoding')
    else:
        base = spec.get('baseline_samps')
        if not base:
            err(name + ': neither baseline_samps nor counts_bc')
        else:
            missing = [b for b in base if b not in known]
            if missing:
                err(name + ': baseline_samps not in count matrix: ' + repr(missing))
            else:
                ok('baseline mode: t0 from ' + repr(base))

    # library
    lib_path = spec.get('library')
    if not lib_path:
        err(name + ': no library')
        return
    try:
        lib = S.load_library(lib_path)
    except SystemExit as e:
        err(name + ': ' + str(e))
        return
    classes = sorted(lib['classification'].dropna().unique())
    present = [c for c in S.CONTROL_CLASSES if c in classes]
    if not present:
        err(name + ': library has none of the control classes ' + repr(S.CONTROL_CLASSES)
            + '. Present classifications: ' + repr(classes) + '. Matching is on exact '
            'strings, so a difference in case or hyphenation yields an empty control '
            'set and bootstrap_null fails.')
    elif len(present) < len(S.CONTROL_CLASSES):
        warn(name + ': only ' + repr(present) + ' of ' + repr(S.CONTROL_CLASSES)
             + ' present in the library')

    merged = pd.merge(df[['sgRNA']].rename(columns={'sgRNA': 'gRNA_id'}),
                      lib, on='gRNA_id')
    if len(merged) != len(df):
        warn(name + ': library inner join would drop ' + str(len(df) - len(merged))
             + ' of ' + str(len(df)) + ' guides')
    else:
        ok('library join keeps all ' + str(len(df)) + ' guides')

    n_ctrl = int(merged['classification'].isin(S.CONTROL_CLASSES).sum())
    if n_ctrl == 0:
        err(name + ': 0 control guides survive the library join -- the null cannot '
            'be built')
    else:
        ok(str(n_ctrl) + ' control guides after the library join')
        for cond in conditions_null:
            if cond in samps:
                n = n_ctrl * len(samps[cond])
                line = cond + ': ' + str(n) + ' pooled control LFCs (' + str(n_ctrl) \
                    + ' guides x ' + str(len(samps[cond])) + ' reps)'
                if n < 100:
                    warn(name + '/' + line + ' -- thin. The p-value floor is still '
                         '1/(n_bootstrap+1), so the bootstrap will look precise '
                         'whether or not the controls support it.')
                else:
                    ok(line)

    # MLE files
    mle_dir = spec.get('mle_dir')
    if not mle_dir:
        ok('no mle_dir -- sensor editing merge will be skipped')
    elif not os.path.isdir(mle_dir):
        err(name + ': mle_dir not found: ' + mle_dir)
    else:
        mle_names = spec.get('mle_names') or list(samps.keys())
        affix = spec.get('mle_affix', '')
        print('        MLE folder for this screen: ' + mle_dir)
        ids = set(merged['gRNA_id'])
        for cond in mle_names:
            fname = cond + '_' + affix + '.csv' if affix else cond + '.csv'
            fpath = os.path.join(mle_dir, fname)
            if not os.path.exists(fpath):
                err(name + ': missing MLE file ' + fpath)
                continue
            e = pd.read_csv(fpath)
            missing_cols = [c for c in S.MLE_COLS if c not in e.columns]
            if missing_cols:
                err(name + '/' + fname + ': missing columns ' + repr(missing_cols))
                continue
            kept = ids & set(e['Guide_ID'])
            lost = len(ids) - len(kept)
            if lost:
                warn(name + '/' + fname + ': inner join would drop ' + str(lost)
                     + ' guides. The editing merges cascade, so this shrinks the '
                     'whole output table, not just these columns.')
            else:
                ok(fname + ': covers all ' + str(len(ids)) + ' guides')
            ids = kept


# ---------------------------------------------------------------------------

# Only fully-spelled, condition-suffixed column names. A bare 'LFC' is almost
# always a locally renamed column, and a trailing-underscore stem like
# 'sensor_reads_' is the prefix half of an f-string -- neither is a real miss.
COL_PAT = re.compile(
    r"""['"]((?:LFC_(?:avg|median|concordant_mean|directional)"""
    r"""|z_score_(?:avg|median|concordant_mean)"""
    r"""|FDR(?:_high|_low)?|p_(?:high|low)_unadjusted"""
    r"""|sensor_reads|corr_perc|target_base_edit_perc"""
    r"""|byproduct_INDEL_perc|byproduct_sub_perc)_[A-Za-z0-9_]+)['"]""")


def check_downstream(table_path, notebooks):
    """Report column names the notebooks reference that the table does not have.

    Advisory, not authoritative: reported as warnings because most figure code
    builds column names with f-strings this grep cannot see, in both directions
    (a miss it cannot find, and a locally renamed column it wrongly flags).
    """
    cols = set(pd.read_csv(table_path, nrows=1).columns)
    print('')
    print('=== downstream contract: ' + table_path + ' ('
          + str(len(cols)) + ' columns) ===')
    for nb_path in notebooks:
        if not os.path.exists(nb_path):
            warn('notebook not found: ' + nb_path)
            continue
        if nb_path.endswith('.ipynb'):
            nb = json.load(open(nb_path))
            text = '\n'.join(''.join(c.get('source', []))
                             for c in nb.get('cells', [])
                             if c.get('cell_type') == 'code')
        else:
            text = open(nb_path).read()
        wanted = set(COL_PAT.findall(text))
        missing = sorted(w for w in wanted if w not in cols)
        if missing:
            warn(os.path.basename(nb_path) + ' references ' + str(len(missing))
                 + ' column(s) this table does not supply: ' + repr(missing[:12])
                 + ('...' if len(missing) > 12 else ''))
        else:
            ok(os.path.basename(nb_path) + ': all ' + str(len(wanted))
               + ' referenced columns present')
    print('')
    print('Note: f-string-built column names (f"LFC_median_{cond}") cannot be')
    print('matched by this grep, so a clean pass is necessary, not sufficient.')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config')
    ap.add_argument('--screen', action='append')
    ap.add_argument('--check-downstream', metavar='LFC_FDR.csv')
    ap.add_argument('--notebooks', nargs='*', default=[])
    args = ap.parse_args(argv)

    if args.config:
        screens = load_config(args.config)
        names = args.screen or list(screens)
        for n in names:
            if n not in screens:
                err('not in config: ' + n)
            else:
                check_screen(n, screens[n])

    if args.check_downstream:
        check_downstream(args.check_downstream, args.notebooks)

    if not args.config and not args.check_downstream:
        ap.error('give --config and/or --check-downstream')

    print('')
    print(str(len(ERRORS)) + ' error(s), ' + str(len(WARNINGS)) + ' warning(s)')
    return 1 if ERRORS else 0


if __name__ == '__main__':
    sys.exit(main())
