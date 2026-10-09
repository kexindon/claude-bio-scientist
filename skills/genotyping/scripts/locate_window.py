"""Find the window in an .ab1 read: start, end, frame and strand, by aligning
the read against an expected coding sequence.

Nothing here is lifted: this is new code, written to replace a manual step in
  /Users/kexindong/Documents/GitHub/hsc-ch-mutation/20260512-HSC-PE7-KP-PK-editing/endogenous_editing.ipynb
The notebook required the user to supply
`start`, `end` and `reverse_complement` by hand for all 15 variants, which is
both the dominant time cost and the source of its three worst failure modes:
a window shifted by a few bases produces a plausible-looking logo of the wrong
codons; a frame offset of 1 nt produces a complete amino-acid logo that is
entirely wrong; and the wrong `reverse_complement` is not detectable from the
primer name.

The alignment makes all three derivable instead of guessed.

Usage
-----
    python locate_window.py --ab1 KP_Sort-p53_F.ab1 \
        --cds-file trp53_exon7.fa --residue 245 --flank-codons 3

    python locate_window.py --ab1 read.ab1 --cds ATGGAGGAG... --residue 245

`--cds` is the in-frame coding sequence of the amplified region, i.e. its
first base is codon position 1. `--cds-first-residue` says which residue that
first codon is (default 1), so residue numbering comes out in protein
coordinates.

Reference CDS can come from $LAB_DATA_ROOT; see docs/data-layout.md.
"""

import argparse
import os

import trace as tr

_COMP = str.maketrans('ACGTUNacgtun', 'TGCAANtgcaan')


def rc(seq):
    return seq.translate(_COMP)[::-1]


def lab_data_root():
    """Reference data lives under $LAB_DATA_ROOT, never a hardcoded path.

    The source notebook hardcoded /Users/kexindong/Desktop/PhD/Experiments/...
    in 15 cells, which is exactly what must not be reproduced here.
    """
    root = os.environ.get("LAB_DATA_ROOT")
    if not root:
        raise SystemExit("LAB_DATA_ROOT is not set -- see docs/data-layout.md")
    return root


def read_fasta_seq(path):
    """First record of a plain FASTA, uppercased, whitespace stripped."""
    if not os.path.exists(path):
        raise SystemExit("No such sequence file: " + str(path))
    lines = []
    started = False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if started:
                    break
                started = True
                continue
            lines.append(line)
    seq = ''.join(lines).upper().replace(' ', '')
    if not seq:
        raise SystemExit("No sequence found in " + str(path))
    return seq


def _best_offset(read, ref, min_anchor=18):
    """Crude but adequate: find where `ref` sits in `read` using exact anchors.

    Sanger reads carry miscalls and a leading N-run, so an exact full-length
    match will not happen. Instead slide a `min_anchor`-mer taken from several
    places in the reference and keep the offset that the most anchors agree on.
    Returns (read_offset_of_ref_base_0, n_anchors_agreeing) or (None, 0).

    An indel between read and reference breaks this by design -- the anchors
    disagree, you get a low agreement count, and you should be using amplicon
    NGS anyway.
    """
    if len(ref) < min_anchor:
        raise SystemExit("Reference is shorter than the anchor length; give more CDS.")
    votes = {}
    step = max(1, (len(ref) - min_anchor) // 12)
    for ref_i in range(0, len(ref) - min_anchor + 1, step):
        anchor = ref[ref_i:ref_i + min_anchor]
        pos = read.find(anchor)
        if pos < 0:
            continue
        offset = pos - ref_i
        votes[offset] = votes.get(offset, 0) + 1
    if not votes:
        return None, 0
    best = max(votes.items(), key=lambda kv: kv[1])
    return best[0], best[1]


def locate(ab1_path, cds, cds_first_residue=1, min_anchor=18):
    """Work out strand and where the CDS starts in the read.

    Returns a dict with `reverse_complement` (what to pass downstream),
    `read_offset`, `n_anchors` and the orientation evidence for both strands.
    """
    read, _phred = tr.called_sequence(ab1_path)
    read = read.upper()
    cds = cds.upper().replace(' ', '').replace('\n', '')

    fwd_off, fwd_n = _best_offset(read, cds, min_anchor)
    rev_off, rev_n = _best_offset(read, rc(cds), min_anchor)

    if fwd_n == 0 and rev_n == 0:
        raise SystemExit(
            "Could not locate the reference in " + os.path.basename(ab1_path) + ".\n"
            "Neither the CDS nor its reverse complement shares an exact " +
            str(min_anchor) + "-mer with the read. Either the reference is the wrong "
            "region, the read is too poor, or there is an indel. Check the called "
            "sequence:\n  " + read[:80] + "..."
        )

    if rev_n > fwd_n:
        # The read is antisense to the CDS: the window must be
        # reverse-complemented before translating.
        return {
            'reverse_complement': True,
            'read_offset': rev_off,
            'n_anchors': rev_n,
            'fwd_anchors': fwd_n,
            'rev_anchors': rev_n,
            'cds_len': len(cds),
            'read_len': len(read),
            'cds_first_residue': cds_first_residue,
        }
    return {
        'reverse_complement': False,
        'read_offset': fwd_off,
        'n_anchors': fwd_n,
        'fwd_anchors': fwd_n,
        'rev_anchors': rev_n,
        'cds_len': len(cds),
        'read_len': len(read),
        'cds_first_residue': cds_first_residue,
    }


def window_for_residue(loc, residue, flank_codons=2):
    """1-based inclusive start/end in the read, for a residue +/- flank codons.

    Frame is derived, not assumed: the CDS's own codon boundaries are carried
    through the offset, so the returned window is guaranteed in frame and a
    multiple of 3.
    """
    first_res = loc['cds_first_residue']
    r0 = residue - flank_codons
    r1 = residue + flank_codons
    if r0 < first_res:
        raise SystemExit(
            "Residue " + str(residue) + " minus " + str(flank_codons) +
            " codons falls before the first residue of the reference (" +
            str(first_res) + "). Give more CDS or fewer flanking codons."
        )

    cds_start0 = (r0 - first_res) * 3
    cds_end0 = (r1 - first_res) * 3 + 2
    if cds_end0 >= loc['cds_len']:
        raise SystemExit(
            "Residue " + str(residue) + " plus flank runs past the end of the "
            "reference CDS. Give more CDS."
        )

    if loc['reverse_complement']:
        # Offsets were computed against rc(cds); convert the rc-coordinate span
        # back to the read's own forward coordinates.
        rc_start0 = loc['cds_len'] - 1 - cds_end0
        rc_end0 = loc['cds_len'] - 1 - cds_start0
        read_start0 = loc['read_offset'] + rc_start0
        read_end0 = loc['read_offset'] + rc_end0
    else:
        read_start0 = loc['read_offset'] + cds_start0
        read_end0 = loc['read_offset'] + cds_end0

    if read_start0 < 0 or read_end0 >= loc['read_len']:
        raise SystemExit(
            "The window for residue " + str(residue) + " falls outside the read "
            "(computed " + str(read_start0 + 1) + ".." + str(read_end0 + 1) +
            " in a read of " + str(loc['read_len']) + " bases). This read does not "
            "cover that residue."
        )

    return {
        'start': read_start0 + 1,
        'end': read_end0 + 1,
        'reverse_complement': loc['reverse_complement'],
        'aa_first': r0,
        'aa_last': r1,
        'n_codons': r1 - r0 + 1,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ab1', required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--cds', help='in-frame coding sequence of the amplified region')
    g.add_argument('--cds-file', help='FASTA holding that sequence')
    p.add_argument('--cds-first-residue', type=int, default=1,
                   help='residue number of the first codon of --cds (default 1)')
    p.add_argument('--residue', type=int, required=True)
    p.add_argument('--flank-codons', type=int, default=2)
    p.add_argument('--min-anchor', type=int, default=18)
    args = p.parse_args()

    cds = args.cds if args.cds else read_fasta_seq(args.cds_file)
    loc = locate(args.ab1, cds, args.cds_first_residue, args.min_anchor)
    win = window_for_residue(loc, args.residue, args.flank_codons)

    read, _ = tr.called_sequence(args.ab1)
    observed = read[win['start'] - 1:win['end']].upper()

    print('ab1                 : ' + os.path.basename(args.ab1))
    print('strand evidence     : ' + str(loc['fwd_anchors']) + ' fwd anchors vs ' +
          str(loc['rev_anchors']) + ' rev anchors')
    print('reverse_complement  : ' + str(win['reverse_complement']))
    print('start..end (1-based): ' + str(win['start']) + '..' + str(win['end']) +
          '  (' + str(win['n_codons']) + ' codons, in frame)')
    print('aa_anno_range       : (' + str(win['aa_first']) + ', ' + str(win['aa_last']) + ')')
    print('called bases        : ' + observed)
    if win['reverse_complement']:
        print('as coding strand    : ' + rc(observed))
    print('')
    print('Confirm the coding-strand bases above translate to the expected WT')
    print('residues before trusting any proportion. Then:')
    print('  python call_edit.py --ab1 ' + args.ab1 +
          ' --start ' + str(win['start']) + ' --end ' + str(win['end']) +
          ' --aa-first ' + str(win['aa_first']) + ' --aa-last ' + str(win['aa_last']) +
          ('  --rc' if win['reverse_complement'] else '') +
          ' --residue ' + str(args.residue) + ' --wt <WT> --alt <ALT> --control <WT.ab1>')


if __name__ == '__main__':
    main()
