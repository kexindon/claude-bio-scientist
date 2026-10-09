#!/usr/bin/env python3
"""STEP3: stage a sequencing run -- config file, library check, folders, sbatch scripts.

Everything that must be right BEFORE any cluster time is spent. Run it, read the
report, then submit.

Lifted from:
  fsrlab-sensor-analysis/NGS-scripts/STEP3_generate_config_file/03_generate_config_file.ipynb
    cell 4 -- the config layout (ArrayTaskID R1_FILE R2_FILE folder_name,
    space-separated, 1-indexed, paths as ./<sample_dir>/<fastq>) and the
    `_guide_split_BARCODES` folder suffix; cells 7/11/12 -- the required-column
    check and the bc_len / sensor_len length reads.
  NGS-scripts/STEP4_sensor_extraction/sensor_extraction_42nt_trimming_quality.sh,
  NGS-scripts/STEP5_sensor_analysis/crispresso_analysis_42nt.sh,
  NGS-scripts/STEP6_sensor_aggregation/crispresso_analysis_aggregation.sh,
  NGS-scripts-prime-editing/STEP4_SENSOR_extraction/sensor_extraction_peg_counts_ext_bc_16nt.sh
    -- the sbatch templates below, header for header.
Do not edit the sbatch headers or the parameter blocks here without editing
those -- these must not drift.

Two things the notebook did that are deliberately NOT reproduced:

  * `k = [i for i in k if i[:3]=='D25']` -- a hardcoded sample-name prefix that
    silently drops every sample in a run whose ids start with anything else.
    Here the prefix is --sample-prefix and defaults to "any directory that
    contains a pair of FASTQs".
  * `if len(j)==2: j[0],j[1] else: j[2],j[3]` -- a positional hack that depends
    on how many files happen to sit in the sample directory. R1/R2 are matched
    on the _1_sequence / _2_sequence token instead.

The shipped .sh files also contain `srun --pty bash` in the body (a copy-paste
artifact from the interactive setup notes -- in a batch script it opens an
interactive shell the job cannot use) and a hardcoded run path, library name and
--mail-user. None of that is emitted here.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import pandas as pd

R1_TOKENS = ('_1_sequence.fastq', '_R1_001.fastq', '_1.fastq')
R2_TOKENS = ('_2_sequence.fastq', '_R2_001.fastq', '_2.fastq')

OUT_FOLDERS = ('counts', 'classification', 'confusion_mats', 'crispresso')


def _lab_root():
    """$LABROOT = the lab volume on Luria. Never a literal /.../kexindong/."""
    root = os.environ.get('LABROOT')
    if not root:
        user = os.environ.get('USER', 'YOURNAME')
        root = '/net/bmc-lab2/data/lab/sanchezrivera/' + user
    return root


# --------------------------------------------------------------------------
# config file
# --------------------------------------------------------------------------

def find_pairs(run_dir, sample_prefix=None):
    """One (sample_dir, r1, r2) per sample directory under run_dir."""
    pairs = []
    problems = []
    for name in sorted(os.listdir(run_dir)):
        sub = os.path.join(run_dir, name)
        if not os.path.isdir(sub):
            continue
        if sample_prefix and not name.startswith(sample_prefix):
            continue
        files = sorted(f for f in os.listdir(sub) if f.endswith('.fastq') or f.endswith('.fastq.gz'))
        r1 = [f for f in files if any(t in f for t in R1_TOKENS)]
        r2 = [f for f in files if any(t in f for t in R2_TOKENS)]
        if len(r1) == 1 and len(r2) == 1:
            pairs.append((name, r1[0], r2[0]))
        elif files:
            problems.append((name, files))
    return pairs, problems


def write_config(run_dir, out_path, sample_prefix=None, suffix='_guide_split_BARCODES',
                 folder_prefix=''):
    pairs, problems = find_pairs(run_dir, sample_prefix)
    if not pairs:
        raise SystemExit('no paired FASTQs found under ' + run_dir
                         + ' -- expected one subdirectory per sample containing '
                           '*_1_sequence.fastq and *_2_sequence.fastq')
    for name, files in problems:
        sys.stderr.write('WARNING: skipped ' + name + ' -- could not identify one R1 and '
                         'one R2 among ' + ', '.join(files) + '\n')

    gz = [n for n, r1, r2 in pairs if r1.endswith('.gz') or r2.endswith('.gz')]
    if gz:
        sys.stderr.write('WARNING: ' + str(len(gz)) + ' sample(s) have gzipped FASTQs. '
                         'Gunzip them, or pass --gz to sensor_extraction.py. The lab\'s '
                         'original scripts hardcode GZ=False and would read a .gz as text, '
                         'giving zero matches with no error.\n')

    rows = []
    for i, (name, r1, r2) in enumerate(pairs, 1):
        rows.append({
            'ArrayTaskID': i,
            'R1_FILE': './' + name + '/' + r1,
            'R2_FILE': './' + name + '/' + r2,
            'folder_name': folder_prefix + name + suffix,
        })
    config_df = pd.DataFrame(rows, columns=['ArrayTaskID', 'R1_FILE', 'R2_FILE', 'folder_name'])
    config_df.to_csv(out_path, index=False, sep=' ')
    return config_df


# --------------------------------------------------------------------------
# library validation
# --------------------------------------------------------------------------

def check_library(lib_path, layout='BE'):
    """Return (lib, bc_len, sensor_len, problems). problems is a list of strings."""
    lib = pd.read_csv(lib_path)
    problems = []
    notes = []

    id_col = 'pegRNA_id' if layout == 'PE' else 'gRNA_id'
    required = [id_col, 'Protospacer', 'Hamming_BC', 'sensor_wt', 'sensor_alt']
    if layout == 'PE':
        required += ['wt_qwc', 'edit_qwc']
    for col in required:
        if col not in lib.columns:
            problems.append('missing required column: ' + col)
    if problems:
        return lib, None, None, problems

    bc_lens = sorted(set(len(str(x)) for x in lib['Hamming_BC']))
    wt_lens = sorted(set(len(str(x)) for x in lib['sensor_wt']))
    alt_lens = sorted(set(len(str(x)) for x in lib['sensor_alt']))

    if len(bc_lens) != 1:
        problems.append('Hamming_BC has mixed lengths: ' + repr(bc_lens))
    if len(wt_lens) != 1:
        problems.append('sensor_wt has mixed lengths: ' + repr(wt_lens)
                        + ' -- sensor_len is one number for the whole run, so every '
                          'sensor_wt must be the same length')

    bc_len = bc_lens[0] if len(bc_lens) == 1 else None
    sensor_len = wt_lens[0] if len(wt_lens) == 1 else None

    # sensor_alt is allowed to differ in length from sensor_wt: for PE, an
    # insertion or deletion edit changes the amplicon length, and the real
    # EA_LUAD library has sensor_alt from 32 to 64 nt against a flat 55 nt
    # sensor_wt. Only sensor_wt sets sensor_len (it is what R1 is trimmed to).
    if alt_lens != wt_lens:
        notes.append('sensor_alt lengths ' + repr(alt_lens) + ' vary against sensor_wt '
                     + repr(wt_lens) + '. Expected for indel edits -- CRISPResso aligns to '
                     'both amplicons independently. Only sensor_wt sets sensor_len.')

    # The scripts compare `Protospacer[1:]` against `r2[:19]`, and lv.hamming on
    # unequal-length strings counts the length difference as mismatches -- so a
    # 21 nt protospacer makes the best possible score 1 and NOTHING ever matches
    # perfectly. Protospacer must be 20 nt (leading G + 19), as the lab's real
    # libraries are.
    proto_lens = sorted(set(len(str(x)) for x in lib['Protospacer']))
    if proto_lens != [20]:
        problems.append(
            'Protospacer lengths are ' + repr(proto_lens) + ', expected 20 (leading G + 19 nt). '
            'The scripts use Protospacer[1:] both as the string compared against R2[:19] and '
            'as CRISPResso\'s --guide_seq. Levenshtein.hamming counts a length difference as '
            'mismatches, so with 21 nt the best achievable distance is 1 and '
            'proto_identified_perfect stays at zero for the whole run.')
    bad_g = [str(x) for x in lib['Protospacer'] if not str(x).upper().startswith('G')]
    if bad_g:
        problems.append(str(len(bad_g)) + ' protospacer(s) do not start with G, e.g. '
                        + bad_g[0] + ' -- the first base is stripped by position, not by '
                        'base, so a protospacer without it loses a real base')

    dup_proto = lib['Protospacer'].duplicated(keep=False)
    if dup_proto.any():
        notes.append(
            str(int(dup_proto.sum())) + ' guide rows share a protospacer with another row '
            '(duplicate_sgRNA=True in counts/). total_guide_count increments EVERY row with '
            'that protospacer, so it is inflated for these; bc_count and '
            'matched_guide_count are per-barcode and unaffected.')

    dup_bc = lib['Hamming_BC'].duplicated(keep=False)
    if dup_bc.any():
        problems.append(str(int(dup_bc.sum())) + ' rows share a Hamming_BC. Barcode matching '
                        'takes the FIRST row with that barcode, so the later ones never '
                        'receive a read.')

    if layout == 'PE':
        qwc_bad = [str(x) for x in lib['wt_qwc'] if not re.match(r'^\d+-\d+$', str(x).strip())]
        if qwc_bad:
            problems.append('wt_qwc values must look like "6-35"; found e.g. ' + qwc_bad[0])

    for col in ('gene_name_m', 'classification', 'Editor'):
        if col not in lib.columns:
            notes.append('no ' + col + ' column -- count_matrix.py needs it '
                         '(gene names, control collapse, subpool filter)')
    if 'classification' in lib.columns:
        ctrl = set(str(x) for x in lib['classification'].unique())
        recognized = {'non-targeting control', 'safe-targeting control'}
        looks_like = [c for c in ctrl if 'control' in c.lower() or 'ntc' in c.lower()]
        unrecognized = [c for c in looks_like if c not in recognized]
        if unrecognized:
            problems.append(
                'control classifications not spelled as expected: ' + repr(unrecognized)
                + '. Only the exact strings "non-targeting control" and "safe-targeting '
                  'control" get collapsed to the gene name controlsgRNA; anything else '
                  'leaves each control guide as its own single-guide gene and breaks the '
                  'control-guide null in screen_stats.')

    return lib, bc_len, sensor_len, problems + ['NOTE: ' + n for n in notes]


# --------------------------------------------------------------------------
# sbatch templates -- faithful to the shipped .sh files
# --------------------------------------------------------------------------

_HEADER = """#!/bin/bash
#SBATCH -N 1                      # Number of nodes. Always 1 for these jobs.
#SBATCH -n 1                      # Number of tasks (CPU cores/task).
#SBATCH --array=1-{n}             # must equal the number of config rows
#SBATCH --mail-type=END
#SBATCH --mail-user={mail}
#SBATCH --exclude=c[5-22]
#SBATCH --nice=100000
#############################################
# Generated by the ngs_screen skill from the lab's STEP{step} template.
# Templates: fsrlab-sensor-analysis/NGS-scripts{pe}/STEP{step}_*/*.sh
#############################################

eval "$(micromamba shell hook --shell bash)"
micromamba activate {env}

cd {run_dir}

config=./{config}

# Extract file names for this array task
R1_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $2}}' $config)
R2_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $3}}' $config)
folder_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $4}}' $config)
"""

_STEP4_BODY = """
# Parameters for the script (NO SPACES ALLOWED)
splitby='{splitby}'                      # 'protospacer' or 'barcode'; use 'protospacer'
proto_mismatches_allowed={mm}
bc_len={bc_len}                          # == len(Hamming_BC), not a free parameter
{real_bc_line}sensor_len={sensor_len}                    # == len(sensor_wt)
quality_check_mode='{qc}'                # full_average | region_average | region_average_and_threshold

python3 sensor_extraction.py {library} ${{R1_FILE}} ${{R2_FILE}} \\
    -o ${{folder_name}} \\
    --layout {layout} \\
    --bc-len ${{bc_len}} {real_bc_arg}--sensor-len ${{sensor_len}} \\
    --splitby ${{splitby}} \\
    --proto-mismatches-allowed ${{proto_mismatches_allowed}} \\
    --quality-check-mode ${{quality_check_mode}}
"""

_STEP5_BODY = """
python3 crispresso_sensors.py run {library} ${{folder_name}} --layout {layout}
"""

_STEP6_BODY = """
python3 crispresso_sensors.py aggregate {library} ${{folder_name}} \\
    --blank crispresso_quant_blank.csv --layout {layout}
"""


def write_sbatch(out_dir, n, run_dir, config, library, mail, layout,
                 bc_len, real_bc_len, sensor_len, splitby, mm, qc):
    pe = '-prime-editing' if layout == 'PE' else ''
    written = []

    real_bc_line = ''
    real_bc_arg = ''
    if layout == 'PE':
        real_bc_line = 'real_bc_len=' + str(real_bc_len) + \
                       '                     # where the sensor really starts in R1\n'
        real_bc_arg = '--real-bc-len ${real_bc_len} '

    specs = [
        ('step4_sensor_extraction.sh', '4', 'sensor_env',
         _STEP4_BODY.format(splitby=splitby, mm=mm, bc_len=bc_len, sensor_len=sensor_len,
                            qc=qc, library=library, layout=layout,
                            real_bc_line=real_bc_line, real_bc_arg=real_bc_arg)),
        ('step5_crispresso.sh', '5', 'crispresso_env',
         _STEP5_BODY.format(library=library, layout=layout)),
        ('step6_aggregation.sh', '6', 'crispresso_env',
         _STEP6_BODY.format(library=library, layout=layout)),
    ]
    for fname, step, env, body in specs:
        text = _HEADER.format(n=n, mail=mail, step=step, pe=pe, env=env,
                              run_dir=run_dir, config=config) + body
        path = os.path.join(out_dir, fname)
        with open(path, 'w') as fh:
            fh.write(text)
        written.append(path)
    return written


# --------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description='Stage a sensor-screen sequencing run: config, library check, '
                    'output folders, sbatch scripts.')
    p.add_argument('--run-dir', required=True,
                   help='the sequencing run folder (one subdirectory per NGS sample). '
                        'Mounted locally this is /Volumes/sanchezrivera/<user>/<run>.')
    p.add_argument('--library', required=True, help='library CSV')
    p.add_argument('--layout', choices=['BE', 'PE'], default='BE')
    p.add_argument('-o', '--out-dir', default='.',
                   help='where to write the config and sbatch scripts (default: cwd)')
    p.add_argument('--config-name', default=None,
                   help='config file name (default CONFIG_<run folder name>.txt)')
    p.add_argument('--sample-prefix', default=None,
                   help='only use sample directories starting with this (e.g. D25)')
    p.add_argument('--folder-suffix', default='_guide_split_BARCODES')
    p.add_argument('--folder-prefix', default='',
                   help='prefix for folder_name, e.g. ./extended_bc_16nt_analysis/ as the '
                        'PE run did')
    p.add_argument('--mail-user', default=None,
                   help='--mail-user for the sbatch headers (default $USER@mit.edu)')
    p.add_argument('--cluster-run-dir', default=None,
                   help='the run directory as seen ON the cluster; default '
                        '$LABROOT/<run folder name>')
    p.add_argument('--splitby', choices=['protospacer', 'barcode'], default='protospacer')
    p.add_argument('--proto-mismatches-allowed', type=int, default=2)
    p.add_argument('--quality-check-mode',
                   choices=['full_average', 'region_average', 'region_average_and_threshold'],
                   default='region_average')
    p.add_argument('--make-folders', action='store_true',
                   help='also create counts/ classification/ confusion_mats/ crispresso/ '
                        'inside --run-dir')
    args = p.parse_args(argv)

    if not os.path.isdir(args.run_dir):
        raise SystemExit('not a directory: ' + args.run_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    run_name = os.path.basename(args.run_dir.rstrip('/'))
    config_name = args.config_name or ('CONFIG_' + run_name + '.txt')
    mail = args.mail_user or (os.environ.get('USER', 'YOURNAME') + '@mit.edu')
    cluster_run_dir = args.cluster_run_dir or (_lab_root() + '/' + run_name)

    print('== library check: ' + args.library)
    lib, bc_len, sensor_len, problems = check_library(args.library, args.layout)
    for msg in problems:
        print('  ' + msg)
    hard = [m for m in problems if not m.startswith('NOTE:')]
    if bc_len is None or sensor_len is None:
        raise SystemExit('cannot determine bc_len / sensor_len from the library; fix the '
                         'problems above first.')

    real_bc_len = bc_len
    if args.layout == 'PE':
        # PE: the matched barcode runs 8 nt into the sensor, so the sensor starts
        # at bc_len - 8. The shipped run used bc_len=16, real_bc_len=8.
        real_bc_len = bc_len - 8
        print('  PE layout: bc_len=' + str(bc_len) + ' (matching), real_bc_len='
              + str(real_bc_len) + ' (sensor start), sensor_len=' + str(sensor_len))
    else:
        print('  BE layout: bc_len=' + str(bc_len) + ', sensor_len=' + str(sensor_len))
    print('  ' + str(len(lib)) + ' guides')

    if hard:
        print('\n' + str(len(hard)) + ' problem(s) above would corrupt the run. Fix the '
              'library, then rerun. Nothing was written.')
        return 1

    print('\n== config file')
    config_df = write_config(args.run_dir, os.path.join(args.out_dir, config_name),
                             args.sample_prefix, args.folder_suffix, args.folder_prefix)
    n = len(config_df)
    print('  ' + os.path.join(args.out_dir, config_name) + ': ' + str(n) + ' samples')
    print('  --array=1-' + str(n))

    if args.make_folders:
        for sub in OUT_FOLDERS:
            os.makedirs(os.path.join(args.run_dir, sub), exist_ok=True)
        print('  created ' + ', '.join(OUT_FOLDERS) + ' in ' + args.run_dir)

    print('\n== sbatch scripts')
    for path in write_sbatch(args.out_dir, n, cluster_run_dir, config_name,
                             os.path.basename(args.library), mail, args.layout,
                             bc_len, real_bc_len, sensor_len, args.splitby,
                             args.proto_mismatches_allowed, args.quality_check_mode):
        print('  ' + path)

    print('\nCopy to ' + cluster_run_dir + ': the config, the library CSV, the three .sh '
          'files, sensor_extraction.py, crispresso_sensors.py, and '
          'assets/crispresso_quant_blank.csv (STEP6 fails without it).')
    print('Then: sbatch step4_sensor_extraction.sh')
    return 0


if __name__ == '__main__':
    sys.exit(main())
