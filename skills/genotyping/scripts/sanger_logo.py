"""Sanger .ab1 trace -> nucleotide / amino-acid sequence logo.

Lifted verbatim from
  /Users/kexindong/Documents/GitHub/hsc-ch-mutation/20260512-HSC-PE7-KP-PK-editing/endogenous_editing.ipynb
  cell 2 (the whole `sanger_weblogo` engine, docstring included).
Do not edit the body of `sanger_weblogo` here without editing there -- these two
must not drift. It is the reference implementation: `trace.py` re-expresses the
same arithmetic and `tests` assert the two agree.

Two deliberate differences from the notebook cell, both about usability and
neither touching the arithmetic:
  * the 'nc' branch returns `nt_df` instead of `None` (in the notebook,
    `aa_df = sanger_weblogo(..., 'nc', ...)` silently binds None);
  * the stray `print(nt_df)` in the 'nc' branch is gone.

logomaker is imported lazily so that `trace.py` / `call_edit.py` consumers can
import this module for its arithmetic on a machine where logomaker is missing.
"""

import os
from itertools import product
from typing import Optional, Tuple

import pandas as pd
from Bio import Seq, SeqIO


def _logomaker():
    """Import logomaker with an error that says what to install."""
    try:
        import logomaker as lm
    except ImportError:
        raise SystemExit(
            "logomaker is not importable in this interpreter.\n"
            "It lives in the system Python 3.9 (/usr/bin/python3, logomaker 0.8.7) "
            "and is absent from every conda env on this machine (base, myenv, pegg, "
            "pegg3, pegg_env, stat_env).\n"
            "Fix with:  pip install 'logomaker>=0.8,<1.0'\n"
            "or run under /usr/bin/python3, or see envs/genotyping_env.yml."
        )
    return lm


# ---------------------------------------------------------------------------
# VERBATIM from endogenous_editing.ipynb cell 2 -- see module docstring.
# ---------------------------------------------------------------------------

def sanger_weblogo(file_path: str,
                   start: int,
                   end: int,
                   logo_type: str = 'nc',
                   reverse_complement: bool = False,
                   aa_anno_range: Optional[Tuple[int, int]] = None,
                   save_path: Optional[str] = None):
    """
    Generate a nucleotide or amino acid WebLogo based on Sanger (.ab1) trace intensities,
    with optional reverse-complement and amino-acid numbering on the X-axis.

    Parameters
    ----------
    file_path : str
        Path to the .ab1 Sanger sequencing file.

    start : int
        1-based start coordinate of the region (inclusive) in AB1-called sequence (PBAS1).

    end : int
        1-based end coordinate of the region (inclusive) in AB1-called sequence (PBAS1).

    logo_type : {'nc','aa'}, default 'nc'
        'nc' : Nucleotide WebLogo. If `aa_anno_range=(a,b)` is given, the X-axis shows AA indices
               (a..b) placed at the center of each codon (every 3 nt).
        'aa' : Amino-acid WebLogo. If `aa_anno_range=(a,b)` is given, the X-axis shows AA indices
               (a..b), one per codon column.

    reverse_complement : bool, default False
        If True, reverse-complement the selected region before plotting.
        - In 'nc', bases are RC and X-axis (if not using aa_anno_range) goes from end->start (nt).
        - In 'aa', codons are RC, translated on reverse strand; ticks reflect RC order.

    aa_anno_range : Optional[Tuple[int,int]], default None
        Amino-acid numbering (inclusive) for X-axis labeling.
        Requirements when provided:
        - Region length (end - start + 1) must be a multiple of 3.
        - Number of codons must equal (b - a + 1).
        Behavior:
        - 'nc': place ticks at the center of each codon and label a..b.
        - 'aa': label each AA column as a..b.
        If None, fallback to nucleotide coordinates:
            'nc' -> real nt coordinates (start..end or reversed if RC)
            'aa' -> codon start nt coordinates (or reversed if RC)

    Notes
    -----
    * Base-wise probabilities are computed from peak trace intensities at PLOC1 using channels A/C/G/T.
    * In 'aa', codon probabilities (all 64) are collapsed to amino-acid probabilities per codon.
    * AA logo uses Logomaker's "chemistry" color scheme (20 AA colors).
    """
    lm = _logomaker()
    import matplotlib.pyplot as plt

    # ---------- Step 1: read .ab1 and get trace-derived nucleotide probabilities ----------
    record = SeqIO.read(file_path, "abi")
    abif = record.annotations['abif_raw']

    channel_order = abif['FWO_1'].decode()  # e.g., "GATC"
    positions = abif['PLOC1']
    channel_data = {b: abif['DATA' + str(i + 9)] for i, b in enumerate(channel_order)}

    # Build per-position base proportions for the requested region
    region_intensities = []
    for idx in range(start - 1, end):  # convert to 0-based
        pos = positions[idx]  # peak index into traces
        intensities = {b: channel_data[b][pos] for b in 'ACGT'}
        total = sum(intensities.values())
        freqs = {b: intensities[b] / total if total > 0 else 0 for b in 'ACGT'}
        region_intensities.append(freqs)

    nt_df = pd.DataFrame(region_intensities)  # rows correspond to start..end (forward order)
    # ---------- Step 2: handle reverse complement if requested ----------
    if reverse_complement:
        # reverse row order and complement base columns
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        nt_df = nt_df.iloc[::-1].reset_index(drop=True)
        nt_df = nt_df.rename(columns=complement)[['A', 'C', 'G', 'T']]  # enforce column order

    region_len = end - start + 1
    if region_len <= 0:
        raise ValueError("`end` must be >= `start`.")

    # ---------- Step 3: nucleotide logo ----------
    if logo_type == 'nc':
        logo = lm.Logo(nt_df, color_scheme="classic")
        # Default coordinate labeling (nt coordinates, forward or reversed)
        if aa_anno_range is None:
            if reverse_complement:
                coords = list(range(end, start - 1, -1))
            else:
                coords = list(range(start, end + 1))
            logo.ax.set_xticks(range(len(coords)))
            logo.ax.set_xticklabels(coords)
        else:
            # Validate that length is multiple of 3 and matches (b-a+1)
            if region_len % 3 != 0:
                raise ValueError(
                    "Region length (" + str(region_len) +
                    ") must be a multiple of 3 when using aa_anno_range."
                )
            a, b = aa_anno_range
            n_codons = region_len // 3
            if (b - a + 1) != n_codons:
                raise ValueError(
                    "aa_anno_range (" + str(a) + "," + str(b) + ") implies " +
                    str(b - a + 1) + " codons, but region has " + str(n_codons) + " codons."
                )
            # Place ticks at the center of each codon: positions 1,4,7,... in 0-based -> 3*i+1
            tick_positions = [3 * i + 1 for i in range(n_codons)]
            tick_labels = list(range(a, b + 1))
            logo.ax.set_xticks(tick_positions)
            logo.ax.set_xticklabels(tick_labels)
        logo.ax.set_ylabel("Proportion")
        plt.tight_layout()
        if save_path:
            _ensure_parent(save_path)
            plt.savefig(save_path, bbox_inches='tight')
        plt.close(logo.ax.figure)
        return nt_df

    # ---------- Step 4: amino acid logo ----------
    if logo_type == 'aa':
        # Must be multiple of 3
        if region_len % 3 != 0:
            raise ValueError(
                "Selected region length (" + str(region_len) + ") is not a multiple of 3. "
                "Please adjust start/end to cover complete codons."
            )

        # Build AA probability rows codon-by-codon (left->right in the plotted orientation)
        aa_rows = []
        n_codons = region_len // 3
        for i in range(n_codons):
            codon_df = nt_df.iloc[i * 3:(i + 1) * 3]  # 3 rows per codon
            aa_probs = {}
            # Accumulate probabilities for all 64 codons
            for c1, c2, c3 in product('ACGT', repeat=3):
                p = codon_df.iloc[0][c1] * codon_df.iloc[1][c2] * codon_df.iloc[2][c3]
                if p == 0:
                    continue
                aa = str(Seq.Seq(c1 + c2 + c3).translate())
                aa_probs[aa] = aa_probs.get(aa, 0) + p
            aa_rows.append(aa_probs)

        aa_df = pd.DataFrame(aa_rows).fillna(0)
        # X-axis labeling
        if aa_anno_range is None:
            # fallback: show nucleotide coordinates of codon starts (respect RC orientation)
            if reverse_complement:
                codon_starts = list(range(end - 2, start - 1, -3))
            else:
                codon_starts = list(range(start, end + 1, 3))
            xticks = range(len(codon_starts))
            xlabels = codon_starts
        else:
            a, b = aa_anno_range
            if (b - a + 1) != n_codons:
                raise ValueError(
                    "aa_anno_range (" + str(a) + "," + str(b) + ") implies " +
                    str(b - a + 1) + " codons, but region has " + str(n_codons) + " codons."
                )
            xticks = range(n_codons)
            xlabels = list(range(a, b + 1))

        # Plot AA logo
        logo = lm.Logo(aa_df, color_scheme="chemistry")
        logo.ax.set_xticks(xticks)
        logo.ax.set_xticklabels(xlabels)
        logo.ax.set_ylabel("Proportion")
        plt.tight_layout()
        if save_path:
            _ensure_parent(save_path)
            plt.savefig(save_path, bbox_inches='tight')
        plt.close(logo.ax.figure)
        return aa_df

    raise ValueError("logo_type must be 'nc' or 'aa'")


# ---------------------------------------------------------------------------
# Wrapper, not part of the lift.
# ---------------------------------------------------------------------------

def _ensure_parent(save_path):
    """The notebook had no makedirs, so a missing `endogenous_plot/` killed a
    long run at the very last savefig. Create it instead."""
    parent = os.path.dirname(os.path.abspath(save_path))
    if parent:
        os.makedirs(parent, exist_ok=True)
