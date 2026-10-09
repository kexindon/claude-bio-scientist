"""ABI .ab1 trace arithmetic, trace-quality guards, and the codon -> amino-acid
collapse, as library functions.

The arithmetic is re-expressed from
  /Users/kexindong/Documents/GitHub/hsc-ch-mutation/20260512-HSC-PE7-KP-PK-editing/endogenous_editing.ipynb
  cell 2, Step 1 (base proportions) and Step 4 (codon collapse).
`sanger_logo.sanger_weblogo` in this same directory is the verbatim lift and is
the reference; this module is the independent re-expression so the two can be
compared. Do not change the arithmetic here without changing it there -- the
check is `tests/test_trace_agrees.py`, which asserts both produce identical
frames on one .ab1.

What is NOT from the notebook, and is new here: `window_quality`, `indel_scan`
and `assert_interpretable`. The notebook reported peak-height proportions
unconditionally; these refuse to.
"""

import os
from itertools import product

import numpy as np
import pandas as pd
from Bio import Seq, SeqIO

BASES = ('A', 'C', 'G', 'T')
_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


# ---------------------------------------------------------------------------
# Reading the trace
# ---------------------------------------------------------------------------

def read_ab1(ab1_path):
    """Open an .ab1 and return (record, abif_raw). Clear error if missing."""
    if not os.path.exists(ab1_path):
        raise SystemExit("No such .ab1 file: " + str(ab1_path))
    record = SeqIO.read(ab1_path, "abi")
    return record, record.annotations['abif_raw']


def called_sequence(ab1_path):
    """The AB1-called base sequence (PBAS1) and its per-base phred scores.

    This is the coordinate system `start`/`end` are expressed in. It is what
    SnapGene or Benchling shows for the read -- NOT a genomic or transcript
    coordinate, and NOT the parallel .seq FASTA trimmed by the facility.
    """
    record, _ = read_ab1(ab1_path)
    return str(record.seq), list(record.letter_annotations['phred_quality'])


def base_proportions(ab1_path, start=None, end=None, reverse_complement=False):
    """Per-position A/C/G/T proportions from the analyzed trace channels.

    `start`/`end` are 1-based inclusive coordinates in the called sequence
    (PBAS1). Omit both for the whole read.

    Channel mapping is read off `FWO_1`, never assumed: on these instruments
    FWO_1 is b'GATC', so DATA9=G, DATA10=A, DATA11=T, DATA12=C. All four
    channels are sampled, keyed by base letter -- an alphabetical DATA9..12 ->
    A,C,G,T assumption would silently transpose G/A and T/C. DATA9-12 are the
    analyzed (baselined, mobility-corrected) traces; DATA1-4 are raw and are
    not what this reads.

    Each channel is sampled at `PLOC1[idx]`, the trace index of the called
    peak for base `idx` -- the peak position, not the base index.
    """
    record, abif = read_ab1(ab1_path)

    channel_order = abif['FWO_1'].decode()
    positions = abif['PLOC1']
    channel_data = {b: abif['DATA' + str(i + 9)] for i, b in enumerate(channel_order)}

    n = len(positions)
    if start is None:
        start = 1
    if end is None:
        end = n
    if start < 1:
        raise SystemExit("start must be >= 1 (it is 1-based); got " + str(start))
    if end < start:
        raise SystemExit("end must be >= start; got start=" + str(start) + " end=" + str(end))
    if end > n:
        raise SystemExit(
            "Window end " + str(end) + " is past the end of the read (" + str(n) +
            " called bases in " + os.path.basename(ab1_path) + ").\n"
            "The notebook had no bounds check here and raised an opaque IndexError "
            "from positions[idx] instead."
        )

    rows = []
    for idx in range(start - 1, end):
        pos = positions[idx]
        intensities = {b: channel_data[b][pos] for b in BASES}
        total = sum(intensities.values())
        rows.append({b: (intensities[b] / total if total > 0 else 0) for b in BASES})

    nt_df = pd.DataFrame(rows, columns=list(BASES))

    if reverse_complement:
        nt_df = nt_df.iloc[::-1].reset_index(drop=True)
        nt_df = nt_df.rename(columns=_COMPLEMENT)[list(BASES)]

    return nt_df


# ---------------------------------------------------------------------------
# Trace quality and indel guards -- NOT in the notebook.
# ---------------------------------------------------------------------------

def window_quality(ab1_path, start, end):
    """Quality summary for one window, as a dict.

    Peak-height proportions are only as good as the trace under them. A window
    inside the leading N-run of a read gives a perfectly plausible-looking logo
    built on noise.

    Reference numbers measured on KP_Sort-p53_F.ab1: in the clean interior
    (bases 100-200) mean phred is ~58 and the mean top-base proportion ~0.91;
    in the leading N-run the top-base proportion sits at 0.5-0.7 with phred
    below 15.
    """
    seq, phred = called_sequence(ab1_path)
    nt_df = base_proportions(ab1_path, start, end)

    sub_phred = phred[start - 1:end]
    sub_seq = seq[start - 1:end]
    top = nt_df.max(axis=1)

    return {
        'ab1': os.path.basename(ab1_path),
        'start': start,
        'end': end,
        'called_bases': sub_seq,
        'n_ambiguous_calls': sum(1 for c in sub_seq if c not in BASES),
        'mean_phred': float(np.mean(sub_phred)) if len(sub_phred) else float('nan'),
        'min_phred': int(np.min(sub_phred)) if len(sub_phred) else -1,
        'frac_phred_below_20': float(np.mean([p < 20 for p in sub_phred])) if len(sub_phred) else float('nan'),
        'mean_top_base_prop': float(top.mean()),
        'min_top_base_prop': float(top.min()),
    }


def indel_scan(ab1_path, window_start, window_end, flank=60,
               mixed_threshold=0.75, run_length=8):
    """Look for an indel upstream of the window that would corrupt it.

    An indel in a non-clonal pool puts the two alleles out of register, so from
    the indel onward every position is a superposition of two different bases.
    The trace does not get noisy-looking -- it gets *consistently mixed*, and
    every proportion downstream becomes meaningless while still looking like a
    clean heterozygote.

    Detection: scan the `flank` bases before the window plus the window itself,
    and flag a run of >= `run_length` consecutive positions whose top-base
    proportion is below `mixed_threshold`. Returns a dict with `suspect_indel`
    and the position where the run starts.

    This is a heuristic, deliberately on the paranoid side. A positive means
    "do not read a percentage off this trace, run amplicon NGS (ngs_screen)" --
    it does not localize or size the indel, and it cannot distinguish an indel
    from a genuinely very mixed pool of substitutions.
    """
    scan_start = max(1, window_start - flank)
    nt_df = base_proportions(ab1_path, scan_start, window_end)
    top = nt_df.max(axis=1).to_numpy()

    best_run = 0
    best_at = None
    run = 0
    for i, v in enumerate(top):
        if v < mixed_threshold:
            run += 1
            if run > best_run:
                best_run = run
                best_at = scan_start + i - run + 1
        else:
            run = 0

    return {
        'suspect_indel': best_run >= run_length,
        'longest_mixed_run': int(best_run),
        'mixed_run_starts_at': best_at,
        'scanned_from': scan_start,
        'scanned_to': window_end,
        'mixed_threshold': mixed_threshold,
        'run_length_threshold': run_length,
    }


def assert_interpretable(ab1_path, start, end, min_mean_phred=30,
                         max_frac_low_phred=0.2, allow_indel=False):
    """Raise unless the window is clean enough to quote a proportion from.

    Call this before reporting any editing percentage. It refuses rather than
    printing a confident wrong number. Returns the merged quality + indel dict
    when it passes.

    `allow_indel=True` downgrades the indel check to a returned flag, for when
    you only want the logo as a picture and are not going to quote a number.
    """
    q = window_quality(ab1_path, start, end)
    ind = indel_scan(ab1_path, start, end)
    q.update(ind)

    problems = []
    if q['mean_phred'] < min_mean_phred:
        problems.append(
            "mean phred over the window is " + ("%.1f" % q['mean_phred']) +
            " (< " + str(min_mean_phred) + ")"
        )
    if q['frac_phred_below_20'] > max_frac_low_phred:
        problems.append(
            ("%.0f%%" % (100 * q['frac_phred_below_20'])) + " of window bases have phred < 20"
        )
    if q['n_ambiguous_calls'] > 0:
        problems.append(
            str(q['n_ambiguous_calls']) + " base(s) in the window are called N"
        )
    if q['suspect_indel'] and not allow_indel:
        problems.append(
            "a run of " + str(q['longest_mixed_run']) + " consecutively mixed positions "
            "starting at base " + str(q['mixed_run_starts_at']) +
            " suggests an indel upstream; everything downstream of an indel is a "
            "superposition of two out-of-register alleles and its proportions are "
            "not an allele fraction"
        )

    if problems:
        raise SystemExit(
            "Refusing to report proportions for " + q['ab1'] +
            " window " + str(start) + ".." + str(end) + ":\n  - " +
            "\n  - ".join(problems) +
            "\nCalled bases in window: " + q['called_bases'] +
            "\nQuantify this one by amplicon NGS instead (the ngs_screen skill), "
            "or pass --allow-indel / --min-mean-phred to override if you only want "
            "the logo as a picture."
        )
    return q


# ---------------------------------------------------------------------------
# Codon -> amino acid collapse
# ---------------------------------------------------------------------------

def codon_aa_proportions(nt_df):
    """Collapse per-base proportions into per-residue amino-acid proportions.

    Re-expressed from endogenous_editing.ipynb cell 2, Step 4: for each codon,
    take the product of the three per-base marginals over all 64 codons,
    translate, and accumulate per amino acid.

    This assumes the three positions of a codon are independent, which is FALSE
    for a real edited pool -- a cell carrying the edit usually carries all of
    its bystanders together. With more than one edited base in a codon the
    product-of-marginals can show amino acids corresponding to haplotypes that
    do not exist. Use `haplotype_warning` to detect that case.
    """
    n = len(nt_df)
    if n % 3 != 0:
        raise SystemExit(
            "Window length " + str(n) + " is not a multiple of 3, so it does not "
            "cover whole codons. Adjust start/end."
        )
    aa_rows = []
    for i in range(n // 3):
        codon_df = nt_df.iloc[i * 3:(i + 1) * 3]
        aa_probs = {}
        for c1, c2, c3 in product('ACGT', repeat=3):
            p = codon_df.iloc[0][c1] * codon_df.iloc[1][c2] * codon_df.iloc[2][c3]
            if p == 0:
                continue
            aa = str(Seq.Seq(c1 + c2 + c3).translate())
            aa_probs[aa] = aa_probs.get(aa, 0) + p
        aa_rows.append(aa_probs)
    return pd.DataFrame(aa_rows).fillna(0)


def haplotype_warning(nt_df, mixed_threshold=0.85):
    """Codon indices (0-based within the window) with >1 mixed base.

    For these codons the independence assumption behind `codon_aa_proportions`
    is doing real work and the amino-acid proportions are an approximation that
    can invent haplotypes. Report them alongside any number taken from such a
    codon.
    """
    top = nt_df.max(axis=1).to_numpy()
    out = []
    for i in range(len(nt_df) // 3):
        n_mixed = int(np.sum(top[i * 3:(i + 1) * 3] < mixed_threshold))
        if n_mixed > 1:
            out.append({'codon_index': i, 'n_mixed_bases': n_mixed})
    return out


def label_aa_frame(aa_df, aa_first):
    """Index an AA proportion frame by real residue number."""
    out = aa_df.copy()
    out.index = list(range(aa_first, aa_first + len(out)))
    out.index.name = 'residue'
    return out
