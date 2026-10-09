#!/usr/bin/env python3
"""Run the LFC / empirical-FDR pipeline for one or many screens, from a config.

The statistics live in `screenstats.py`, lifted verbatim from
  fsrlab-sensor-analysis/NGS-scripts/STEP9_LFC_FDR_calculation/LFC_FDR_calculation.ipynb
This file is only the config loader, the CLI and the provenance sidecar -- it
replaces the notebook's cells 17/18/20/21/22 (load counts, define replicate
groups, define `screens`, loop, save), so that nobody has to edit a notebook to
add a screen or to drop a replicate.

    python run_screen.py --config screens.yaml
    python run_screen.py --config screens.yaml --screen ABE_OG --out LFC-FDR
    python run_screen.py --counts counts/og/OG_ABE_COUNTS.txt \
        --library MBESv2_CORRECTED.csv \
        --samps 'bm=bonemarrow1,bonemarrow2,bonemarrow3' \
        --samps 'spleen=spleen1,spleen2,spleen3,spleen4' \
        --baseline-samps plasmidlib --conditions-null bm,spleen \
        --name ABE_OG --out LFC-FDR

Writes, per screen, under --out:
    <screen>_LFC_FDR.csv              the wide per-guide table
    <screen>_null_dist.pdf            null QC panel (plain; paper figures -> sci-plot)
    <screen>_LFC_FDR.provenance.json  seed, n_bootstrap, replicate lists, MLE folder,
                                      control counts, input file hashes
"""

import argparse
import datetime
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pandas as pd  # noqa: E402

import screenstats as S  # noqa: E402


# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------

CONFIG_DOC = """
Config format (YAML if PyYAML is installed, otherwise JSON -- same keys).

defaults:                       # optional; merged into every screen
  library: MBESv2_focused.csv   # or  lab:MBESv2_focused  -> $LAB_DATA_ROOT/libraries/
  mle_affix: ''
  n_bootstrap: 10000
  seed: 0

screens:
  ABE_FC_EPO_T0:
    counts: counts/focused/FC_ABE_EPO_COUNTS.txt
    library: MBESv2_focused.csv
    samps:                      # condition -> replicate column names.
      bm: [bm1, bm3, bm4, bm5]  #   bm2 excluded (abnormal skewing, global_analysis QC)
      d15: [d15-rep1, d15-rep2, d15-rep3, d15-rep4, d15-rep5]
      d5: [d5]
      input: [input]
      lib: [lib]
      spleen: [spleen1, spleen3, spleen4, spleen5]   # spleen2 excluded
    baseline_samps: [input]     # T0 mode
    conditions_null: [d5, d15, spleen, bm]
    mle_dir: MLE/focused/ABE_EPO     # omit entirely to skip the editing merge
    mle_names: [input, d5, d15, bm, spleen]
    mle_affix: ''

  ABE_FC_EPO_TO_BC:             # barcoding baseline instead of T0
    counts: counts/focused/FC_ABE_EPO_COUNTS.txt
    counts_bc: counts/focused/FC_ABE_BC_COUNTS.txt
    samps: {...}
    samps_bc: {...}             # same condition keys as samps
    conditions_null: [d5, d15, spleen, bm]
    # NOTE: do NOT also give baseline_samps -- that combination is rejected.
"""


def load_config(path):
    with open(path) as fh:
        text = fh.read()
    if path.endswith(('.yaml', '.yml')):
        try:
            import yaml
        except ImportError:
            raise SystemExit('PyYAML is not installed; use a .json config '
                             '(same keys) or `pip install pyyaml`')
        cfg = yaml.safe_load(text)
    else:
        cfg = json.loads(text)
    if not isinstance(cfg, dict) or 'screens' not in cfg:
        raise SystemExit("config must be a mapping with a 'screens' key.\n" + CONFIG_DOC)
    defaults = cfg.get('defaults') or {}
    screens = {}
    for name, spec in cfg['screens'].items():
        merged = dict(defaults)
        merged.update(spec or {})
        screens[name] = merged
    return screens


def parse_samps(pairs):
    """['bm=bm1,bm3', 'spleen=spleen1'] -> {'bm': ['bm1','bm3'], ...}"""
    out = {}
    for p in pairs:
        if '=' not in p:
            raise SystemExit("--samps expects 'condition=rep1,rep2', got " + repr(p))
        cond, reps = p.split('=', 1)
        out[cond.strip()] = [r.strip() for r in reps.split(',') if r.strip()]
    return out


def read_counts(path):
    sep = ',' if path.endswith('.csv') else '\t'
    df = pd.read_csv(path, sep=sep)
    cols = list(df.columns)
    if cols[:2] != ['sgRNA', 'gene']:
        raise SystemExit(
            path + ": first two columns must be exactly ['sgRNA', 'gene'], got "
            + repr(cols[:2]) + ". rpm() treats every column from the third onward "
            "as a sample, so a stray index column from a to_csv/read_csv round "
            "trip corrupts every RPM value. Re-export without the index.")
    return df


def file_digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()[:16]


# ---------------------------------------------------------------------------
# one screen
# ---------------------------------------------------------------------------

def run_one(name, spec, out_dir, show=False, legacy_loop=False):
    counts_path = spec['counts']
    df_counts = read_counts(counts_path)
    library = S.load_library(spec['library'])

    samps = spec['samps']
    new_names = list(samps.keys())
    samps_to_merge = [samps[k] for k in new_names]

    known = set(df_counts.columns)
    missing = [r for g in samps_to_merge for r in g if r not in known]
    if missing:
        raise SystemExit('replicate columns absent from ' + counts_path + ': '
                         + repr(missing))

    kwargs = dict(
        screen_name=name,
        df_counts=df_counts,
        library=library,
        samps_to_merge=samps_to_merge,
        new_names=new_names,
        conditions_null=spec['conditions_null'],
        n_bootstrap=spec.get('n_bootstrap', 10000),
        seed=spec.get('seed', 0),
        ambiguity_threshold=spec.get('ambiguity_threshold', 2/3),
        null_color=spec.get('null_color', 'tab:purple'),
        save_path=os.path.join(out_dir, name + '_null_dist.pdf'),
        show=show,
        legacy_loop=legacy_loop,
    )

    counts_bc_path = spec.get('counts_bc')
    if counts_bc_path:
        df_bc = read_counts(counts_bc_path)
        samps_bc = spec['samps_bc']
        kwargs['df_counts_bc'] = df_bc
        kwargs['new_names_bc'] = list(samps_bc.keys())
        kwargs['samps_to_merge_bc'] = [samps_bc[k] for k in samps_bc]
        if spec.get('baseline_samps'):
            raise SystemExit(name + ': both counts_bc and baseline_samps given. '
                             'The notebook silently ignored the T0 baseline in '
                             'that case; pick one.')
    else:
        if not spec.get('baseline_samps'):
            raise SystemExit(name + ": needs either baseline_samps (T0 mode) or "
                             "counts_bc + samps_bc (barcoding mode)")
        kwargs['to_comparison_samps'] = spec['baseline_samps']

    mle_dir = spec.get('mle_dir')
    if mle_dir:
        kwargs['fp_editing'] = mle_dir
        kwargs['mle_names'] = spec.get('mle_names') or new_names
        kwargs['mle_affix'] = spec.get('mle_affix', '')

    table, null_dict, controls = S.run_screen_pipeline(**kwargs)

    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, name + '_LFC_FDR.csv')
    table.to_csv(out_csv, index=False)
    print('  wrote ' + out_csv + '  (' + str(table.shape[0]) + ' guides x '
          + str(table.shape[1]) + ' columns)')

    prov = {
        'screen': name,
        'written': datetime.datetime.now().isoformat(timespec='seconds'),
        'engine': 'screenstats.py (lifted from STEP9 LFC_FDR_calculation.ipynb)',
        'baseline': 'barcoding' if counts_bc_path else 't0',
        'baseline_samps': spec.get('baseline_samps'),
        'counts': counts_path,
        'counts_sha256_16': file_digest(counts_path),
        'counts_bc': counts_bc_path,
        'library': spec['library'],
        'replicates_per_condition': samps,
        'conditions_null': spec['conditions_null'],
        'mle_dir': mle_dir,
        'mle_names': spec.get('mle_names'),
        'n_bootstrap': spec.get('n_bootstrap', 10000),
        'seed': spec.get('seed', 0),
        'ambiguity_threshold': spec.get('ambiguity_threshold', 2/3),
        'control_classes': S.CONTROL_CLASSES,
        'n_guides_out': int(table.shape[0]),
        'n_columns_out': int(table.shape[1]),
        'n_control_guides': int(len(controls)),
        'n_control_lfcs_per_condition': dict(
            (c, int(len(controls) * len(samps[c]))) for c in spec['conditions_null']),
        'legacy_loop': bool(legacy_loop),
    }
    prov_path = os.path.join(out_dir, name + '_LFC_FDR.provenance.json')
    with open(prov_path, 'w') as fh:
        json.dump(prov, fh, indent=2)
    print('  wrote ' + prov_path)
    return table


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=CONFIG_DOC)
    ap.add_argument('--config', help='YAML/JSON config with a screens: block')
    ap.add_argument('--screen', action='append',
                    help='only run this screen from the config (repeatable)')
    ap.add_argument('--out', default='LFC-FDR', help='output directory')
    ap.add_argument('--show', action='store_true',
                    help='display the null QC panel instead of only writing the PDF')
    ap.add_argument('--legacy-loop', action='store_true',
                    help='use the original pure-Python empirical-p double loop '
                         '(minutes per condition) instead of the verified '
                         'searchsorted form -- for equivalence checking')
    ap.add_argument('--print-config-doc', action='store_true',
                    help='print the config format and exit')

    # single-screen flag form
    ap.add_argument('--counts')
    ap.add_argument('--counts-bc')
    ap.add_argument('--library')
    ap.add_argument('--samps', action='append', default=[],
                    metavar='COND=rep1,rep2', help='repeatable')
    ap.add_argument('--samps-bc', action='append', default=[], metavar='COND=rep1,rep2')
    ap.add_argument('--baseline-samps', help='comma-separated T0 baseline columns')
    ap.add_argument('--conditions-null', help='comma-separated conditions to test')
    ap.add_argument('--mle-dir')
    ap.add_argument('--mle-names', help='comma-separated')
    ap.add_argument('--mle-affix', default='')
    ap.add_argument('--name', default='screen')
    ap.add_argument('--n-bootstrap', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args(argv)

    if args.print_config_doc:
        print(CONFIG_DOC)
        return 0

    if args.config:
        screens = load_config(args.config)
        if args.screen:
            unknown = [s for s in args.screen if s not in screens]
            if unknown:
                raise SystemExit('not in config: ' + repr(unknown)
                                 + '; available: ' + repr(sorted(screens)))
            screens = dict((s, screens[s]) for s in args.screen)
    else:
        if not (args.counts and args.library and args.samps and args.conditions_null):
            ap.error('either --config, or all of --counts --library --samps '
                     '--conditions-null')
        spec = {
            'counts': args.counts,
            'library': args.library,
            'samps': parse_samps(args.samps),
            'conditions_null': [c.strip() for c in args.conditions_null.split(',')],
            'n_bootstrap': args.n_bootstrap,
            'seed': args.seed,
        }
        if args.baseline_samps:
            spec['baseline_samps'] = [c.strip() for c in args.baseline_samps.split(',')]
        if args.counts_bc:
            spec['counts_bc'] = args.counts_bc
            spec['samps_bc'] = parse_samps(args.samps_bc)
        if args.mle_dir:
            spec['mle_dir'] = args.mle_dir
            spec['mle_affix'] = args.mle_affix
            if args.mle_names:
                spec['mle_names'] = [c.strip() for c in args.mle_names.split(',')]
        screens = {args.name: spec}

    for name in screens:
        run_one(name, screens[name], args.out, show=args.show,
                legacy_loop=args.legacy_loop)
    return 0


if __name__ == '__main__':
    sys.exit(main())
