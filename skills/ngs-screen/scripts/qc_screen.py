#!/usr/bin/env python3
"""Diagnostics for a sensor-screen run: which tasks finished, and why reads were lost.

Three subcommands, each answering a question that comes up on every real run and
that the lab's scripts do not answer:

    qc_screen.py report  --run-dir DIR
        Per sample: % good quality, % no_match_bc, recombination rate, zero-count
        guides, count skew. Reads classification/ and counts/, which the docs
        describe column by column but nothing ever summarizes.

    qc_screen.py failed  --run-dir DIR --config CONFIG
        Which ArrayTaskIDs have no counts/<folder>_count_df.csv, printed as a
        ready-to-paste --array=3,17,58 list. A 188-task array on a shared
        filesystem routinely loses a few tasks, and "which samples failed,
        rerun those" is the most common question about this pipeline.

    qc_screen.py strand  --library LIB --r1 R1.fastq [--layout BE|PE]
        Takes the first N reads of one R1 and reports how many match the
        library's Hamming_BC as given, versus reverse-complemented. The
        extraction script reverse-complements Hamming_BC before matching, so a
        library that already stores the R1-strand barcode sends every read to
        no_match_bc and the counts come back all zero. This says which strand
        the library is on instead of just reporting a bad number.

No code is lifted here -- the lab's scripts have no diagnostics. The column
names and outcome rows come from
  fsrlab-sensor-analysis/docs/source/extraction.md ("Actions before running STEP7")
and from the dataframes built in STEP4's extraction_filtration.

Plots: this script prints tables only. If you want these as a figure panel, hand
the numbers to the sci-plot skill -- no skill in this repo defines its own
figure style.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

OUTCOMES = ['good quality', 'low_qual_r1', 'low_qual_r2', 'low_qual_r12',
            'no_match_bc', 'bc_identified', 'proto_identified_perfect',
            'proto_identified_imperfect', 'proto_identified_recombined',
            'no_match_proto']

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(seq):
    return seq.translate(_COMP)[::-1]


def _pct(num, den):
    if not den:
        return float('nan')
    return 100.0 * float(num) / float(den)


def cmd_report(args):
    cls_dir = os.path.join(args.run_dir, 'classification')
    cnt_dir = os.path.join(args.run_dir, 'counts')
    if not os.path.isdir(cls_dir):
        raise SystemExit('no classification/ under ' + args.run_dir)

    rows = []
    for fname in sorted(os.listdir(cls_dir)):
        if not fname.endswith('_classification_df.csv'):
            continue
        sample = fname.replace('_classification_df.csv', '')
        cls = pd.read_csv(os.path.join(cls_dir, fname)).set_index('classification')['count']

        # 'good quality' + the three low_qual rows are the read-level partition.
        # no_match_bc and the proto_* rows are counted INSIDE good quality, so
        # they are percentages of good/bc_identified, not of the read total.
        total = float(cls.reindex(OUTCOMES[:4]).fillna(0).sum())
        good = float(cls.get('good quality', 0))
        bc_id = float(cls.get('bc_identified', 0))
        row = {
            'sample': sample,
            'reads': int(total),
            'good_%': _pct(good, total),
            'no_match_bc_%': _pct(cls.get('no_match_bc', 0), good),
            'bc_id_%': _pct(bc_id, good),
            'perfect_%': _pct(cls.get('proto_identified_perfect', 0), bc_id),
            'recomb_%': _pct(cls.get('proto_identified_recombined', 0), bc_id),
        }

        cpath = os.path.join(cnt_dir, sample + '_count_df.csv')
        if os.path.exists(cpath):
            cnt = pd.read_csv(cpath)
            col = args.count_col if args.count_col in cnt.columns else 'bc_count'
            v = np.asarray(cnt[col], dtype=float)
            row['guides'] = len(v)
            row['zero_count'] = int((v == 0).sum())
            tot = v.sum()
            if tot > 0 and len(v) >= 10:
                top = np.sort(v)[::-1][:max(1, len(v) // 10)].sum()
                row['top10%_share'] = 100.0 * top / tot
        rows.append(row)

    if not rows:
        raise SystemExit('no classification CSVs found -- has STEP4 run?')
    rep = pd.DataFrame(rows)
    with pd.option_context('display.width', 200, 'display.max_columns', 50):
        print(rep.to_string(index=False, float_format=lambda x: '%.1f' % x))

    print('')
    bad_bc = rep[rep['no_match_bc_%'] > 50]
    if len(bad_bc):
        print('>50% of good-quality reads matched no barcode in '
              + str(len(bad_bc)) + ' sample(s). Run:  qc_screen.py strand  -- the usual '
              'cause is that Hamming_BC is stored on the R1 strand already, while the '
              'extraction script reverse-complements it before matching.')
    hi_recomb = rep[rep['recomb_%'] > 20]
    if len(hi_recomb):
        print(str(len(hi_recomb)) + ' sample(s) over 20% protospacer-barcode recombination. '
              'bc_count, matched_guide_count and total_guide_count will diverge badly; pick '
              'matched_guide_count for the stringent read, and look at confusion_mats/ '
              'off-diagonal mass before trusting PE counts.')
    print('bc_count / matched_guide_count / total_guide_count are not interchangeable. '
          'Pick one and keep it fixed across a paper.')
    return 0


def cmd_failed(args):
    config = pd.read_csv(args.config, sep=' ')
    cnt_dir = os.path.join(args.run_dir, 'counts')
    agg_dir = os.path.join(args.run_dir, 'crispresso')

    missing_counts = []
    missing_agg = []
    for _, r in config.iterrows():
        tag = os.path.basename(str(r['folder_name']).rstrip('/'))
        if not os.path.exists(os.path.join(cnt_dir, tag + '_count_df.csv')):
            missing_counts.append(int(r['ArrayTaskID']))
        elif not os.path.exists(os.path.join(agg_dir, tag + '_crispresso_aggregated.csv')):
            missing_agg.append(int(r['ArrayTaskID']))

    n = len(config)
    print(str(n - len(missing_counts)) + '/' + str(n) + ' samples have counts/ (STEP4 done)')
    if missing_counts:
        print('STEP4 incomplete, resubmit with:')
        print('  #SBATCH --array=' + ','.join(str(i) for i in missing_counts))
        print('  Per-guide FASTQs are appended to, so clear each failed sample\'s folder '
              'first, or pass --force to sensor_extraction.py. Rerunning into a populated '
              'folder double-counts reads.')
    if missing_agg:
        print('STEP4 done but no aggregated CSV (STEP5/6 incomplete):')
        print('  #SBATCH --array=' + ','.join(str(i) for i in missing_agg))
        print('  STEP5/6 are safe to rerun; they overwrite.')
    if not missing_counts and not missing_agg:
        print('Nothing missing. Copy counts/ and crispresso/*_crispresso_aggregated.csv '
              'locally, then run compile_editing.py.')
    return 0


def _r1_reads(path, n, gz=None):
    if gz is None:
        gz = path.endswith('.gz')
    _open = gzip.open if gz else open
    mode = 'rt' if gz else 'r'
    out = []
    with _open(path, mode) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.strip())
                if len(out) >= n:
                    break
    return out


def cmd_strand(args):
    lib = pd.read_csv(args.library)
    if 'Hamming_BC' not in lib.columns:
        raise SystemExit('library has no Hamming_BC column')
    bcs = [str(x).upper() for x in lib['Hamming_BC']]
    bc_lens = sorted(set(len(b) for b in bcs))
    if len(bc_lens) != 1:
        raise SystemExit('Hamming_BC has mixed lengths ' + repr(bc_lens))
    bc_len = bc_lens[0]

    as_given = set(bcs)
    as_rc = set(rc(b) for b in bcs)

    reads = _r1_reads(args.r1, args.n_reads)
    if not reads:
        raise SystemExit('no reads read from ' + args.r1)

    hit_given = 0
    hit_rc = 0
    for r in reads:
        head = r[:bc_len].upper()
        if head in as_given:
            hit_given += 1
        if head in as_rc:
            hit_rc += 1

    print('R1: ' + args.r1)
    print(str(len(reads)) + ' reads, barcode taken as R1[:' + str(bc_len) + ']')
    print('  matches Hamming_BC as given:            %6.2f%%' % _pct(hit_given, len(reads)))
    print('  matches revcomp(Hamming_BC):            %6.2f%%' % _pct(hit_rc, len(reads)))
    print('')
    print('The extraction script matches against revcomp(Hamming_BC), so the second number '
          'is the one that becomes bc_identified.')
    if hit_rc > hit_given * 2 and _pct(hit_rc, len(reads)) > 5:
        print('=> Library is on the sense strand, as the pipeline expects. Good.')
    elif hit_given > hit_rc * 2 and _pct(hit_given, len(reads)) > 5:
        print('=> Library already stores the R1-strand barcode. The script will '
              'reverse-complement it again, so every read lands in no_match_bc and counts '
              'come back all zero. Reverse-complement the Hamming_BC column in the library, '
              'or drop the revcomp in extraction_filtration -- not both.')
    else:
        print('=> Neither strand matches well. Check instead: (a) is there an adapter or '
              'index before the barcode in R1, (b) is bc_len right, (c) is this the right '
              'library file for this run. For PE, Hamming_BC must already be the 16 nt '
              'matching barcode (8 nt real BC + first 8 nt of sensor).')

    if 'Protospacer' in lib.columns and args.r2:
        protos = set(str(x)[1:].upper() for x in lib['Protospacer'])
        r2 = _r1_reads(args.r2, args.n_reads)
        hit = sum(1 for r in r2 if r[:19].upper() in protos)
        print('')
        print('R2: ' + args.r2)
        print('  R2[:19] matches Protospacer[1:]:        %6.2f%%' % _pct(hit, len(r2)))
        p_obs = sorted(set(len(str(x)) for x in lib['Protospacer']))
        print('  Protospacer lengths in the library: ' + repr(p_obs) + ' (must be 20)')
        if p_obs != [20]:
            print('  => Protospacer[1:] is hamming-compared against R2[:19]. '
                  'Levenshtein.hamming counts a length difference as a mismatch, so with '
                  'anything but 20 nt the best achievable distance is nonzero and '
                  'proto_identified_perfect stays at zero for the whole run.')
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)

    r = sub.add_parser('report', help='per-sample read accounting and count skew')
    r.add_argument('--run-dir', required=True,
                   help='directory holding classification/ and counts/')
    r.add_argument('--count-col', default='bc_count')
    r.set_defaults(func=cmd_report)

    f = sub.add_parser('failed', help='which array tasks did not finish')
    f.add_argument('--run-dir', required=True)
    f.add_argument('--config', required=True)
    f.set_defaults(func=cmd_failed)

    s = sub.add_parser('strand', help='is Hamming_BC on the strand the pipeline expects')
    s.add_argument('--library', required=True)
    s.add_argument('--r1', required=True)
    s.add_argument('--r2', default=None, help='also check Protospacer against R2[:19]')
    s.add_argument('--n-reads', type=int, default=10000)
    s.set_defaults(func=cmd_strand)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
