#!/usr/bin/env python3
"""STEP4 sensor extraction: paired FASTQ -> per-guide counts + per-guide sensor FASTQs.

Lifted verbatim from
  fsrlab-sensor-analysis/NGS-scripts/STEP4_sensor_extraction/sensor_extraction_42nt_trimming_quality.py
    (fastq_reader, quality_checker and its three modes, the MIN_QUALITY /
     MAX_LOWQ_BASES / MIN_QUALITY_PER_BASE globals, to_IOSeq_rec,
     extraction_filtration's read loop and both output dataframes)
and
  fsrlab-sensor-analysis/NGS-scripts-prime-editing/STEP4_SENSOR_extraction/sensor_extraction_peg_counts_ext_bc_16nt.py
    (the real_bc_len variant of to_IOSeq_rec, the pegRNA_id column).
Do not edit the matching logic here without editing there -- these must not drift.

The two source scripts are near-duplicates; this is one file covering both
layouts, selected by --layout:

    base-editing (BE):   bc_len = real_bc_len = 15, sensor_len = 42
    prime-editing (PE):  bc_len = 16 (matching), real_bc_len = 8 (sensor start),
                         sensor_len = 55

Four bugs in the originals are FIXED here; each fix is marked `# FIX:` inline so
it is visible rather than silent. Nothing else about the matching, counting or
classification logic is changed.

  1. BE argv off-by-one. The BE script checks `len(sys.argv) != 10` and then
     reads `sys.argv[10]`, which does not exist when argv has 10 entries --
     an unconditional IndexError before any work starts. (The PE sibling is
     correct: `!= 12` then `argv[11]`, because `-o` occupies a slot.) This
     file uses argparse, so the question cannot come back.
  2. BE never applies quality_check_mode. The BE script calls
     `quality_checker(q1, q2)` although the signature is
     `quality_checker(q1, q2, bc_len, sensor_len, proto_len=19, mode=...)`:
     the two required positional args are missing, so it raises TypeError on
     the first read. If someone patched it by giving them defaults, the
     quality_check_mode read from argv would be ignored and every run would
     silently fall back to region_average. All four are passed here.
  3. PE extracts the wrong sensor window. In the PE script
     `to_IOSeq_rec(sensor_seq, i, q1, real_bc_len, sensor_len)` slices
     `sensor_seq[real_bc_len:real_bc_len+sensor_len]`, but all three call
     sites pass `bc_len` (16) into that parameter. Every PE sensor read was
     therefore taken 8 nt downstream of the real start and ran 8 nt past the
     window, so the per-guide FASTQs did not align to the 55 nt sensor_wt
     amplicon. real_bc_len is passed here.
  4. One open()/append/close per read. Both originals reopened the per-guide
     FASTQ in 'a' mode for every single read. Writes are buffered per guide
     and flushed in batches, which is the same bytes in the same order --
     only far fewer syscalls. This is the pipeline's hot loop.

Also: os.makedirs(exist_ok=True) plus a refusal to write into a folder that
already holds FASTQs, so a rerun of a failed array task says so instead of
either crashing (BE: os.mkdir -> FileExistsError) or silently double-counting
(per-guide files opened in append mode). Use --force to clear and redo.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

import Bio.Seq
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord

import Levenshtein as lv

# ---------------------------------------------------------------------------
# GLOBAL VARIABLES -- lifted verbatim. Changing these changes which reads pass
# QC, so they stay named constants rather than becoming flags.
# ---------------------------------------------------------------------------
MIN_QUALITY = 30          # originally set to 30
MAX_LOWQ_BASES = 5        # originally set to 5
MIN_QUALITY_PER_BASE = 20  # originally set to 20

QC_MODES = ('full_average', 'region_average', 'region_average_and_threshold')

# How many records to hold per guide before flushing. Memory is
# n_guides * FLUSH_EVERY * ~150 bytes; 64 over 5,000 guides is ~50 MB.
FLUSH_EVERY = 64


# THIS FUNCTION MAY NEED TO BE UPDATED DEPENDING ON HOW THE DATA IS STRUCTURED
def fastq_reader(fname, gz=False):
    _open = gzip.open if gz else open
    proc_read = (lambda line: line.strip().decode()) if gz else (lambda line: line.strip())
    proc_qual = (lambda line: line.strip()) if gz else (lambda line: line.strip().encode('utf-8'))
    with _open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                identifier = proc_read(line)
            elif i % 4 == 1:
                read = proc_read(line)
            elif i % 4 == 3:
                quality = proc_qual(line)
                yield identifier[-20:], read, quality  # protospacer (index2), read, quality score


def quality_checker(q1, q2, bc_len, sensor_len, proto_len=19,
                    mode='region_average',
                    min_q=MIN_QUALITY,  # originally set to 30
                    max_lowq_bases=MAX_LOWQ_BASES,  # originally set to 5
                    min_q_per_base=MIN_QUALITY_PER_BASE  # originally set to 20
                    ):
    """
    Check read quality using different QC strategies.

    Parameters
    ----------
    q1, q2 : bytes
        FASTQ quality strings for R1 and R2.
    bc_len : int
        Barcode length in R1 (the MATCHING length: 15 for BE, 16 for PE).
    sensor_len : int
        Sensor length in R1.
    proto_len : int
        Protospacer length in R2.
    mode : str
        QC mode. One of 'full_average', 'region_average',
        'region_average_and_threshold'.
    min_q : int
        Minimum Phred score threshold.
    max_lowq_bases : int
        Maximum allowed number of low-quality bases in required region.
    min_q_per_base : int
        Per-base Phred floor used by the threshold mode.
    """

    q1_scores = np.frombuffer(q1, dtype=np.uint8) - 33
    q2_scores = np.frombuffer(q2, dtype=np.uint8) - 33

    low_qual = []

    if mode == 'full_average':
        # Original behavior: check average quality across full reads
        qual1 = q1_scores.mean()
        qual2 = q2_scores.mean()

        if qual1 < min_q:
            low_qual.append('1')
        if qual2 < min_q:
            low_qual.append('2')

    elif mode == 'region_average':
        # Only check regions actually used downstream
        r1_region = q1_scores[:bc_len + sensor_len]
        r2_region = q2_scores[:proto_len]

        qual1 = r1_region.mean()
        qual2 = r2_region.mean()

        if qual1 < min_q:
            low_qual.append('1')
        if qual2 < min_q:
            low_qual.append('2')

    elif mode == 'region_average_and_threshold':
        # Only check regions actually used downstream
        r1_region = q1_scores[:bc_len + sensor_len]
        r2_region = q2_scores[:proto_len]
        bc_region = q1_scores[:bc_len]

        # Mean quality in required regions
        qual1 = r1_region.mean()
        qual2 = r2_region.mean()

        # Barcode is especially important for correct assignment
        if np.min(bc_region) < min_q:
            low_qual.append('1')

        # Average quality check
        if qual1 < min_q:
            if '1' not in low_qual:
                low_qual.append('1')

        if qual2 < min_q:
            if '2' not in low_qual:
                low_qual.append('2')

        # Per-base low-quality threshold check
        if np.sum(r1_region < min_q_per_base) > max_lowq_bases:
            if '1' not in low_qual:
                low_qual.append('1')

        if np.sum(r2_region < min_q_per_base) > max_lowq_bases:
            if '2' not in low_qual:
                low_qual.append('2')

    else:
        raise ValueError("Unknown QC mode: " + str(mode))

    if len(low_qual) == 0:
        quality_out = 'good quality'
    else:
        quality_out = 'low_qual_r' + ''.join(low_qual)

    return quality_out


'''Prepare fastq records.'''
def to_IOSeq_rec(sensor_seq, i, q1, real_bc_len, sensor_len):
    """
    Extract only the sensor, starting where the sensor really begins.

    real_bc_len is where the sensor starts in R1 -- 15 for the BE layout
    (barcode then sensor, no overlap) and 8 for the PE extended-barcode
    layout (the 16 nt matched barcode overlaps the sensor by 8 nt).
    """
    qname = 'read_' + str(i)
    record = SeqRecord(Bio.Seq.Seq(sensor_seq[real_bc_len:real_bc_len + sensor_len]),
                       id=qname, name=qname, description='', dbxrefs=[])

    # add quality score to the record
    qual = (np.frombuffer(q1, dtype=np.uint8) - 33)
    record.letter_annotations["phred_quality"] = qual[real_bc_len:real_bc_len + sensor_len]

    return record


class _GuideWriter:
    """Buffered per-guide FASTQ writer.

    FIX 4: the originals did `with open(out_file, 'a') as fq: SeqIO.write(...)`
    inside the read loop -- one open, one write, one close per read. Same output
    bytes, batched.
    """

    def __init__(self, folder_name, guide_ids, flush_every=FLUSH_EVERY):
        self.folder_name = folder_name
        self.flush_every = flush_every
        self.buffers = {}
        # write empty fastq files to a new folder for EACH gRNA
        for i in guide_ids:
            f_name = i + '.fastq'
            with open(os.path.join(folder_name, f_name), 'w') as fp:
                pass
            self.buffers[i] = []

    def add(self, guide_id, record):
        buf = self.buffers[guide_id]
        buf.append(record)
        if len(buf) >= self.flush_every:
            self._flush_one(guide_id)

    def _flush_one(self, guide_id):
        buf = self.buffers[guide_id]
        if not buf:
            return
        with open(os.path.join(self.folder_name, guide_id + '.fastq'), 'a') as fq:
            Bio.SeqIO.write(buf, fq, 'fastq')
        del buf[:]

    def close(self):
        for guide_id in self.buffers:
            self._flush_one(guide_id)


def extraction_filtration(folder_name, input_df, R1_FILE, R2_FILE, bc_len, sensor_len,
                          splitby, proto_mismatches_allowed, id_col='gRNA_id',
                          real_bc_len=None, quality_check_mode='region_average',
                          gz=False, breakpoint=False):

    """
    Takes in reads and returns dataframe containing pegRNA counts AND sensor outcomes.

    Also returns summary of outcomes including:
    - low quality (and which of the reads are low quality (or if all are low quality))
    - no extension match
    - no protospacer match
    - decoupled extension-protospacer
    - correct identification

    bc_len is the MATCHING length (what is compared against the barcode list);
    real_bc_len is where the sensor starts and defaults to bc_len (the BE case).
    """

    if real_bc_len is None:
        real_bc_len = bc_len

    # ----------1. generate folders for holding split sensor reads---------------
    # FIX: the BE original did a bare os.mkdir(), which raises FileExistsError
    # on any rerun. The caller has already decided whether reuse is safe.
    os.makedirs(folder_name, exist_ok=True)

    # need to modify based on where the guide_ID column is/what its name is
    writer = _GuideWriter(folder_name, list(input_df[id_col]))

    # ----------------lists and dicts for matching up reads------------------------
    proto_list = [i[1:] for i in input_df['Protospacer']]
    bc_list = [str(Bio.Seq.Seq(i).reverse_complement()) for i in input_df['Hamming_BC']]

    proto_bc_dict = dict(zip(bc_list, proto_list))
    confusion_mat = np.zeros((len(proto_list), len(bc_list)))

    # ------initialize a dataframe for holding the pegRNA counts and sensor outcomes
    d1 = pd.DataFrame(dict(zip(['Guide_ID', 'sgRNA_no_Gstart', 'unique_BC'],
                               [list(input_df[id_col]), proto_list, bc_list])))

    cols = ['total_guide_count', 'matched_guide_count', 'bc_count']
    z = np.zeros((len(cols), len(input_df)))
    d2 = pd.DataFrame(dict(zip(cols, z)))
    df1 = pd.concat((d1, d2), axis=1).set_index('Guide_ID')

    # and annotating the duplicate guides
    i, v = np.unique(d1['sgRNA_no_Gstart'], return_counts=True)

    duplicate_dict = dict(zip(i, v > 1))
    df1['duplicate_sgRNA'] = [duplicate_dict[i] for i in df1['sgRNA_no_Gstart']]

    # ------initialize a dataframe for holding metadata about the identification of sensors
    outcomes = ['good quality', 'low_qual_r1', 'low_qual_r2', 'low_qual_r12',
                'no_match_bc',
                'bc_identified',
                'proto_identified_perfect',
                'proto_identified_imperfect',
                'proto_identified_recombined',
                'no_match_proto',
                ]

    outcomes_count = np.zeros(len(outcomes))
    class_df = pd.DataFrame(dict(zip(['classification', 'count'],
                                     [outcomes, outcomes_count]))).set_index('classification')

    # ----------iterating through the reads...--------------
    for i, ((proto1, r1, q1), (proto2, r2, q2)) in enumerate(
            zip(fastq_reader(R1_FILE, gz=gz), fastq_reader(R2_FILE, gz=gz)), 1):
        # r1 = sensor read
        # r2 = extension read
        # proto1/2 = protospacer read (index2)

        # first check the quality
        # FIX 2: all four args passed, so quality_check_mode is actually applied.
        quality_out = quality_checker(q1, q2, bc_len, sensor_len, mode=quality_check_mode)
        class_df.loc[quality_out, 'count'] += 1

        if quality_out == 'good quality':
            proto = r2[:19]
            bc = r1[:bc_len]

            if bc in bc_list:
                df1.loc[df1['unique_BC'] == bc, 'bc_count'] += 1
                class_df.loc['bc_identified', 'count'] += 1

                bc_idx = bc_list.index(bc)
                proto_match = proto_bc_dict[bc]

                hamming_distance = lv.hamming(proto, proto_match)

                # perfect protospacer match
                if hamming_distance == 0:
                    df1.loc[df1['unique_BC'] == bc, 'matched_guide_count'] += 1
                    df1.loc[df1['unique_BC'] == bc, 'total_guide_count'] += 1
                    class_df.loc['proto_identified_perfect', 'count'] += 1

                    if splitby != 'barcode':
                        # add the sensor read to the appropriate fastq file
                        guide_id = df1[df1['unique_BC'] == bc].index[0]
                        # FIX 3: real_bc_len, not bc_len.
                        writer.add(guide_id, to_IOSeq_rec(r1, i, q1, real_bc_len, sensor_len))

                # imperfect match within allowed hamming distance
                elif hamming_distance <= proto_mismatches_allowed:
                    df1.loc[df1['unique_BC'] == bc, 'matched_guide_count'] += 1
                    df1.loc[df1['unique_BC'] == bc, 'total_guide_count'] += 1
                    class_df.loc['proto_identified_imperfect', 'count'] += 1

                    if splitby != 'barcode':
                        guide_id = df1[df1['unique_BC'] == bc].index[0]
                        writer.add(guide_id, to_IOSeq_rec(r1, i, q1, real_bc_len, sensor_len))

                # nonmatching (beyond allowed hamming distance)
                elif hamming_distance > proto_mismatches_allowed:

                    if proto in proto_list:
                        df1.loc[df1['sgRNA_no_Gstart'] == proto, 'total_guide_count'] += 1
                        class_df.loc['proto_identified_recombined', 'count'] += 1

                        proto_idx = proto_list.index(proto)
                        confusion_mat[proto_idx][bc_idx] += 1

                    else:
                        # and look for imperfect matches
                        found = False
                        for k in proto_list:
                            h_dist = lv.hamming(k, proto)
                            if h_dist <= proto_mismatches_allowed:
                                df1.loc[df1['sgRNA_no_Gstart'] == k, 'total_guide_count'] += 1
                                class_df.loc['proto_identified_recombined', 'count'] += 1

                                proto_idx = proto_list.index(k)
                                confusion_mat[proto_idx][bc_idx] += 1
                                found = True
                                break

                        # FIX: the original tested `if count==len(proto_list)` after
                        # incrementing count inside the loop, so a read that matched
                        # on the LAST library entry was counted both as recombined
                        # and as no_match_proto. An explicit flag cannot be off by one.
                        if not found:
                            class_df.loc['no_match_proto', 'count'] += 1

                # -------if splitting sensors into different fastq files ONLY by barcode-----
                if splitby == 'barcode':
                    guide_id = df1[df1['unique_BC'] == bc].index[0]
                    writer.add(guide_id, to_IOSeq_rec(r1, i, q1, real_bc_len, sensor_len))

            elif bc not in bc_list:
                class_df.loc['no_match_bc', 'count'] += 1

            else:
                continue

        if breakpoint is not False:
            if i > breakpoint:
                break

    writer.close()

    return df1, class_df, confusion_mat


# ---------------------------------------------------------------------------
# CLI -- a wrapper around the lifted code, not a change to it.
# ---------------------------------------------------------------------------

def _looks_populated(folder_name):
    if not os.path.isdir(folder_name):
        return False
    for f in os.listdir(folder_name):
        if f.endswith('.fastq') and os.path.getsize(os.path.join(folder_name, f)) > 0:
            return True
    return False


def main(argv=None):
    p = argparse.ArgumentParser(
        description='STEP4: extract sensors and count guides from one paired-end sample.')
    p.add_argument('library', help='library CSV (gRNA_id or pegRNA_id, Protospacer, Hamming_BC)')
    p.add_argument('r1', help='R1 FASTQ (sensor read), uncompressed unless --gz')
    p.add_argument('r2', help='R2 FASTQ (protospacer read)')
    p.add_argument('-o', '--folder-name', required=True,
                   help='per-sample output folder for the per-guide FASTQs')
    p.add_argument('--layout', choices=['BE', 'PE'], default='BE',
                   help='BE: 15 nt barcode + 42 nt sensor. PE: 16 nt matching barcode '
                        'overlapping the sensor by 8 nt, 55 nt sensor.')
    p.add_argument('--bc-len', type=int, default=None,
                   help='barcode MATCHING length; must equal len(Hamming_BC)')
    p.add_argument('--real-bc-len', type=int, default=None,
                   help='where the sensor starts in R1 (PE only; defaults to --bc-len)')
    p.add_argument('--sensor-len', type=int, default=None,
                   help='sensor length; must equal len(sensor_wt)')
    p.add_argument('--splitby', choices=['protospacer', 'barcode'], default='protospacer')
    p.add_argument('--proto-mismatches-allowed', type=int, default=2)
    p.add_argument('--quality-check-mode', choices=list(QC_MODES), default='region_average')
    p.add_argument('--id-col', default=None,
                   help='guide id column; default gRNA_id for BE, pegRNA_id for PE')
    p.add_argument('--out-root', default='.',
                   help='directory holding counts/, classification/, confusion_mats/')
    p.add_argument('--gz', action='store_true',
                   help='read gzipped FASTQs. The originals hardcoded GZ=False and would '
                        'silently return zero matches on a .gz file.')
    p.add_argument('--breakpoint', type=int, default=0,
                   help='stop after N reads (smoke test); 0 = whole file')
    p.add_argument('--force', action='store_true',
                   help='clear an existing per-sample folder instead of refusing')
    args = p.parse_args(argv)

    defaults = {'BE': (15, 15, 42, 'gRNA_id'), 'PE': (16, 8, 55, 'pegRNA_id')}
    d_bc, d_real, d_sensor, d_id = defaults[args.layout]
    bc_len = args.bc_len if args.bc_len is not None else d_bc
    real_bc_len = args.real_bc_len if args.real_bc_len is not None else d_real
    sensor_len = args.sensor_len if args.sensor_len is not None else d_sensor
    id_col = args.id_col or d_id
    if args.bc_len is not None and args.real_bc_len is None and args.layout == 'BE':
        real_bc_len = bc_len

    input_df = pd.read_csv(args.library)
    if id_col not in input_df.columns:
        raise SystemExit('library has no column ' + id_col + ' -- columns are: '
                         + ', '.join(map(str, input_df.columns)))

    # bc_len / sensor_len are not free parameters. Check them against the library
    # rather than letting a mismatch show up as all-zero counts hours later.
    bc_obs = sorted(set(len(str(x)) for x in input_df['Hamming_BC']))
    if len(bc_obs) != 1:
        raise SystemExit('Hamming_BC has mixed lengths ' + repr(bc_obs) + ' -- fix the library')
    if args.layout == 'BE' and bc_obs[0] != bc_len:
        raise SystemExit('--bc-len ' + str(bc_len) + ' != len(Hamming_BC) '
                         + str(bc_obs[0]) + '; for BE these must match')
    if args.layout == 'PE' and bc_obs[0] != bc_len:
        raise SystemExit('--bc-len ' + str(bc_len) + ' != len(Hamming_BC) ' + str(bc_obs[0])
                         + '; for the PE extended-BC layout the library Hamming_BC must '
                           'already be the 16 nt matching barcode (real BC + 8 nt of sensor)')
    if 'sensor_wt' in input_df.columns:
        s_obs = sorted(set(len(str(x)) for x in input_df['sensor_wt']))
        if len(s_obs) == 1 and s_obs[0] != sensor_len:
            raise SystemExit('--sensor-len ' + str(sensor_len) + ' != len(sensor_wt) '
                             + str(s_obs[0]) + '; CRISPResso will not align')

    # Protospacer[1:] is compared against r2[:19]. Levenshtein.hamming counts a
    # length difference as mismatches, so a protospacer that is not 20 nt makes
    # proto_identified_perfect unreachable for the entire run.
    p_obs = sorted(set(len(str(x)) for x in input_df['Protospacer']))
    if p_obs != [20]:
        raise SystemExit(
            'Protospacer lengths are ' + repr(p_obs) + ', expected 20 (leading G + 19 nt). '
            'Protospacer[1:] is hamming-compared against R2[:19]; an unequal length counts '
            'as a mismatch, so nothing would ever match perfectly. Run setup_run.py to see '
            'the full library report.')

    if _looks_populated(args.folder_name):
        if not args.force:
            raise SystemExit(
                args.folder_name + ' already holds per-guide FASTQs with reads in them.\n'
                'Per-guide FASTQs are appended to, so rerunning into this folder would '
                'double-count. Delete it, or pass --force to clear and redo this sample.')
        for f in os.listdir(args.folder_name):
            if f.endswith('.fastq'):
                os.remove(os.path.join(args.folder_name, f))

    for sub in ('counts', 'classification', 'confusion_mats'):
        os.makedirs(os.path.join(args.out_root, sub), exist_ok=True)

    count_df, class_df, mat = extraction_filtration(
        args.folder_name, input_df, args.r1, args.r2, bc_len, sensor_len,
        args.splitby, args.proto_mismatches_allowed, id_col=id_col,
        real_bc_len=real_bc_len, quality_check_mode=args.quality_check_mode,
        gz=args.gz, breakpoint=(args.breakpoint or False))

    tag = os.path.basename(args.folder_name.rstrip('/'))
    count_df.to_csv(os.path.join(args.out_root, 'counts', tag + '_count_df.csv'))
    class_df.to_csv(os.path.join(args.out_root, 'classification', tag + '_classification_df.csv'))
    pd.DataFrame(mat).to_csv(
        os.path.join(args.out_root, 'confusion_mats', tag + '_confusion_matrix.csv'))

    # Only 'good quality' + the three low_qual rows partition the reads. The
    # no_match_bc and proto_* rows are tallied INSIDE good quality, so summing
    # the whole column would count most reads twice.
    total = float(class_df.loc[['good quality', 'low_qual_r1', 'low_qual_r2',
                                'low_qual_r12'], 'count'].sum())
    good = float(class_df.loc['good quality', 'count'])
    nomatch = float(class_df.loc['no_match_bc', 'count'])
    sys.stderr.write(
        tag + ': ' + str(int(total)) + ' reads, '
        + '%.1f%% good quality, ' % (100.0 * good / total if total else 0.0)
        + '%.1f%% no_match_bc\n' % (100.0 * nomatch / good if good else 0.0))
    if good and nomatch / good > 0.5:
        sys.stderr.write(
            'WARNING: over half of good-quality reads matched no barcode. Usually the '
            'library stores Hamming_BC on the wrong strand -- run qc_screen.py '
            'strand-check before launching the rest of the array.\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
