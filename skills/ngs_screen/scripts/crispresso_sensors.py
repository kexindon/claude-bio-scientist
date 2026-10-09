#!/usr/bin/env python3
"""STEP5 + STEP6: run CRISPResso2 per guide, then aggregate one sample's output.

Two subcommands, because they are two array jobs in two different environments
but share the library-row -> guide-id logic:

    crispresso_sensors.py run       LIBRARY SAMPLE_FOLDER --layout BE|PE
    crispresso_sensors.py aggregate LIBRARY SAMPLE_FOLDER --blank crispresso_quant_blank.csv

Lifted verbatim from:
  fsrlab-sensor-analysis/NGS-scripts/STEP5_sensor_analysis/crispresso_analysis_w_qwc_42nt.py
    -- the BE flag string, character for character:
       --quantification_window_coordinates 6-35,6-35 --exclude_bp_from_left 0
       --exclude_bp_from_right 0 --base_editor_output --plot_window_size 5
       --suppress_report --suppress_plots
  fsrlab-sensor-analysis/NGS-scripts-prime-editing/STEP5_SENSOR_ANALYSIS/crispresso_analysis_55nt.py
    -- the PE flag string: per-row {wt_qwc},{edit_qwc}, --prime_editor_output,
       --plot_window_size 15
  fsrlab-sensor-analysis/NGS-scripts/STEP6_sensor_aggregation/crispresso_analysis_aggregation.py
    -- the aggregation loop, the Reference/HDR blank-row substitution, and the
       <sample>_crispresso_aggregated.csv output
Do not change a CRISPResso flag here without changing it there -- editing
percentages are flag-dependent and the two must not drift.

Changed from the originals, deliberately:

  * os.system(command) -> subprocess.run(list). The originals built one shell
    string and noted "apparently this is bad practice (but subprocess module
    didn't work...)". A list argv needs no quoting and a non-zero return code
    is visible instead of discarded.
  * Guides whose per-guide FASTQ is missing or empty are skipped with a count,
    rather than spending a CRISPResso startup on a file with no reads. STEP6's
    blank-row fallback already reports these as 0 reads.
  * STEP6 originally did os.listdir(.../CRISPResso_on_<guide_id>) with no
    existence check, so ONE guide without an output directory raised
    FileNotFoundError and the whole sample's aggregation was lost -- after all
    the CRISPResso compute had already been paid for. A missing directory now
    takes the same all-zero blank rows as a missing quantification file, and the
    count of substituted guides is printed so it cannot pass unnoticed.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

import pandas as pd

LAYOUTS = {
    'BE': {'id_col': 'gRNA_id', 'editor_flag': '--base_editor_output', 'plot_window': '5'},
    'PE': {'id_col': 'pegRNA_id', 'editor_flag': '--prime_editor_output', 'plot_window': '15'},
}

# The BE quantification window, hardcoded in the lab's STEP5 script: 5 nt of
# flank either side of a 20 nt protospacer inside the 42 nt sensor. Change the
# sensor design and this window quantifies the wrong bases, with no error.
BE_QWC = '6-35,6-35'


def _crispresso_cmd(fastq, wt, edited, proto, qwc, editor_flag, plot_window,
                    out_dir, name):
    """One CRISPResso invocation, as the lab's STEP5 scripts build it."""
    return [
        'CRISPResso',
        '--fastq_r1', fastq,
        '--amplicon_seq', wt,
        '--expected_hdr_amplicon_seq', edited,
        '--guide_seq', proto,
        '--quantification_window_coordinates', qwc,
        '--exclude_bp_from_left', '0',
        '--exclude_bp_from_right', '0',
        editor_flag,
        '--plot_window_size', plot_window,
        '--suppress_report',
        '--suppress_plots',
        '-o', out_dir,
        '-n', name,
    ]


def cmd_run(args):
    spec = LAYOUTS[args.layout]
    id_col = args.id_col or spec['id_col']
    df = pd.read_csv(args.library)
    if id_col not in df.columns:
        raise SystemExit('library has no column ' + id_col)
    if args.layout == 'PE':
        for col in ('wt_qwc', 'edit_qwc'):
            if col not in df.columns:
                raise SystemExit(
                    'PE layout needs a ' + col + ' column: the quantification window is '
                    'per-pegRNA for PE, unlike BE where it is the fixed ' + BE_QWC + '.')

    sample_name = args.sample_folder.rstrip('/')
    out_dir = os.path.join(args.crispresso_dir, os.path.basename(sample_name))

    n_run = 0
    n_empty = 0
    n_failed = 0
    for _, val in df.iterrows():
        guide_id = val[id_col]
        wt = val['sensor_wt']
        edited = val['sensor_alt']
        proto = val['Protospacer'][1:]  # 19 bp (excluding G start)

        fastq = os.path.join(sample_name, str(guide_id) + '.fastq')
        if (not os.path.exists(fastq)) or os.path.getsize(fastq) == 0:
            n_empty += 1
            continue

        if args.layout == 'PE':
            qwc = str(val['wt_qwc']).strip() + ',' + str(val['edit_qwc']).strip()
        else:
            qwc = BE_QWC

        cmd = _crispresso_cmd(fastq, wt, edited, proto, qwc,
                              spec['editor_flag'], spec['plot_window'],
                              out_dir, str(guide_id))
        if args.dry_run:
            print(' '.join(cmd))
            n_run += 1
            continue
        proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            n_failed += 1
            if n_failed <= 5:
                sys.stderr.write('CRISPResso failed on ' + str(guide_id) + ': '
                                 + proc.stderr.decode('utf-8', 'replace').strip()[-300:] + '\n')
        else:
            n_run += 1

    sys.stderr.write(os.path.basename(sample_name) + ': ' + str(n_run) + ' guides run, '
                     + str(n_empty) + ' skipped (no reads), ' + str(n_failed) + ' failed\n')
    if n_failed:
        sys.stderr.write(
            'If these are permission errors on the shared filesystem: '
            'chmod -R g+w $LABROOT/<run>/crispresso\n')
    return 0


def cmd_aggregate(args):
    spec = LAYOUTS[args.layout]
    id_col = args.id_col or spec['id_col']
    df = pd.read_csv(args.library)
    if id_col not in df.columns:
        raise SystemExit('library has no column ' + id_col)

    if not os.path.exists(args.blank):
        raise SystemExit(
            args.blank + ' not found. crispresso_quant_blank.csv -- two all-zero rows, '
            'Reference and HDR -- is what makes a guide with no CRISPResso output appear '
            'as 0 reads instead of vanishing from the table. Ship it next to this script; '
            'a copy is in the skill\'s assets/.')
    quant_zero = pd.read_csv(args.blank)

    sample_name = os.path.basename(args.sample_folder.rstrip('/'))
    fp = args.crispresso_dir  # path to output folders

    # --------run it for all of the guides
    rows = []
    n_blank_missing_dir = 0
    n_blank_missing_file = 0
    for _, val in df.iterrows():
        guide_id = val[id_col]
        guide_dir = os.path.join(fp, sample_name, 'CRISPResso_on_' + str(guide_id))

        if not os.path.isdir(guide_dir):
            # The original os.listdir() here raised FileNotFoundError and killed
            # the whole sample. Same blank rows as a missing quantification file.
            quant = quant_zero.copy()
            quant['Guide_ID'] = guide_id
            n_blank_missing_dir += 1
            rows.append(quant)
            continue

        output_folder_x = os.listdir(guide_dir)

        if 'CRISPResso_quantification_of_editing_frequency.txt' in output_folder_x:
            quant = pd.read_csv(
                os.path.join(guide_dir, 'CRISPResso_quantification_of_editing_frequency.txt'),
                sep='\t')
            quant['Guide_ID'] = guide_id

        else:
            quant = quant_zero.copy()
            quant['Guide_ID'] = guide_id
            n_blank_missing_file += 1

        rows.append(quant)

    output = pd.concat(rows)
    out_path = os.path.join(fp, sample_name + '_crispresso_aggregated.csv')
    output.to_csv(out_path, index=False)

    n_blank = n_blank_missing_dir + n_blank_missing_file
    sys.stderr.write(sample_name + ': ' + str(len(df)) + ' guides -> ' + out_path + '\n')
    if n_blank:
        sys.stderr.write(
            '  ' + str(n_blank) + ' guides got all-zero blank rows ('
            + str(n_blank_missing_dir) + ' with no CRISPResso output directory, '
            + str(n_blank_missing_file) + ' with a directory but no quantification file). '
            'A handful is normal for guides with no reads; a large fraction means STEP5 '
            'did not finish for this sample.\n')
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--layout', choices=['BE', 'PE'], default='BE')
    p.add_argument('--id-col', default=None,
                   help='guide id column; default gRNA_id (BE) / pegRNA_id (PE)')
    p.add_argument('--crispresso-dir', default='./crispresso',
                   help='directory holding one subdirectory per sample (default ./crispresso)')
    sub = p.add_subparsers(dest='cmd', required=True)

    r = sub.add_parser('run', help='STEP5: one CRISPResso run per guide for one sample')
    r.add_argument('library')
    r.add_argument('sample_folder', help='the per-guide FASTQ folder STEP4 wrote')
    r.add_argument('--dry-run', action='store_true', help='print the commands, run nothing')
    r.set_defaults(func=cmd_run)

    a = sub.add_parser('aggregate', help='STEP6: one aggregated CSV per sample')
    a.add_argument('library')
    a.add_argument('sample_folder')
    a.add_argument('--blank', default='crispresso_quant_blank.csv')
    a.set_defaults(func=cmd_aggregate)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
