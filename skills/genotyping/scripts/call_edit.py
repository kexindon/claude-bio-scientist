"""Turn trace proportions into a number: editing fraction, bystanders, and a
homozygous / heterozygous / WT / ambiguous clone call. Batches a plate from a
variant table.

The notebook (endogenous_editing.ipynb) rendered `df_aa` and read it by eye;
there was no calling code. The arithmetic under this is `trace.py`
(re-expressed from cell 2) and `sanger_logo.sanger_weblogo` (verbatim lift) --
this module only thresholds and tabulates.

Two things here exist because peak-height proportions are not calibrated:
  * --control: a matched unedited trace, whose proportion at the same residue
    is subtracted as background. Without one, an "editing %" from a trace is
    not interpretable.
  * every reported number carries the quality/indel verdict from
    trace.assert_interpretable, which refuses rather than guessing.

Usage
-----
One sample:
    python call_edit.py --ab1 KP_Sort-p53_F.ab1 --start 118 --end 138 \
        --aa-first 249 --aa-last 255 --rc \
        --residue 252 --wt L --alt P \
        --control KP_Pre-p53_F.ab1 -o OUT

A plate, from a variant table (the 15 copy-pasted notebook blocks become rows):
    python call_edit.py --table plate.csv -o OUT

Variant table columns:
    ab1, variant, start, end, reverse_complement, aa_first, aa_last
  optional:
    residue, wt_aa, alt_aa, control_ab1
"""

import argparse
import os

import pandas as pd

import trace as tr

try:
    from sanger_logo import sanger_weblogo
except ImportError:  # running from another cwd
    sanger_weblogo = None


# ---------------------------------------------------------------------------
# Calling
# ---------------------------------------------------------------------------

def call_residue(aa_df_labeled, residue, wt_aa=None, alt_aa=None,
                 bystander_threshold=0.05, control_aa_df=None):
    """Read one residue out of a labeled AA proportion frame.

    `aa_df_labeled` is indexed by real residue number (trace.label_aa_frame).
    `control_aa_df` is the same frame from a matched unedited trace; when given,
    each amino-acid proportion is reported both raw and background-subtracted.
    """
    if residue not in aa_df_labeled.index:
        raise SystemExit(
            "Residue " + str(residue) + " is not in the window (" +
            str(aa_df_labeled.index.min()) + ".." + str(aa_df_labeled.index.max()) + ")."
        )
    row = aa_df_labeled.loc[residue]
    ctrl = None
    if control_aa_df is not None:
        if residue not in control_aa_df.index:
            raise SystemExit(
                "Residue " + str(residue) + " is not in the control window; the "
                "control read must be given the same start/end/frame as the sample."
            )
        ctrl = control_aa_df.loc[residue]

    def prop(aa):
        if aa is None:
            return float('nan')
        return float(row.get(aa, 0.0))

    def prop_corrected(aa):
        if aa is None or ctrl is None:
            return float('nan')
        return float(row.get(aa, 0.0)) - float(ctrl.get(aa, 0.0))

    top_aa = row.idxmax()
    bystanders = {}
    for aa in row.index:
        if aa in (wt_aa, alt_aa):
            continue
        v = float(row[aa])
        if ctrl is not None:
            v = v - float(ctrl.get(aa, 0.0))
        if v >= bystander_threshold:
            bystanders[aa] = round(v, 4)

    return {
        'residue': residue,
        'wt_aa': wt_aa,
        'alt_aa': alt_aa,
        'top_aa': top_aa,
        'top_aa_prop': round(float(row[top_aa]), 4),
        'wt_prop': round(prop(wt_aa), 4) if wt_aa else float('nan'),
        'alt_prop': round(prop(alt_aa), 4) if alt_aa else float('nan'),
        'wt_prop_ctrl_corrected': round(prop_corrected(wt_aa), 4) if wt_aa else float('nan'),
        'alt_prop_ctrl_corrected': round(prop_corrected(alt_aa), 4) if alt_aa else float('nan'),
        'has_control': ctrl is not None,
        'bystanders': bystanders,
    }


def genotype_call(call, hom_min=0.85, het_range=(0.25, 0.75), wt_max=0.10):
    """Classify one clone from a residue call.

    Thresholds are on the control-corrected alt proportion when a control was
    given, otherwise on the raw proportion -- and in that case the call is
    marked low-confidence, because an uncalibrated peak-height ratio is not an
    allele fraction.

    Classes: 'homozygous_edited', 'heterozygous', 'wild_type', 'ambiguous'.
    Anything between the bands is 'ambiguous' on purpose: a clone at 0.80 is
    not a clean homozygote and should not be filed as one.
    """
    alt = call['alt_prop_ctrl_corrected'] if call['has_control'] else call['alt_prop']
    if alt != alt:  # NaN
        return {'genotype': 'ambiguous', 'basis_alt_prop': alt,
                'confidence': 'none', 'reason': 'no alt amino acid specified'}

    conf = 'calibrated-vs-control' if call['has_control'] else 'uncalibrated'
    if alt >= hom_min:
        g = 'homozygous_edited'
    elif het_range[0] <= alt <= het_range[1]:
        g = 'heterozygous'
    elif alt <= wt_max:
        g = 'wild_type'
    else:
        g = 'ambiguous'

    reason = ''
    if g == 'ambiguous':
        reason = ('alt proportion ' + ("%.3f" % alt) + ' falls between the het band ' +
                  str(het_range) + ' and the hom cutoff ' + str(hom_min))
    if not call['has_control']:
        reason = (reason + '; ' if reason else '') + \
            'no matched WT control trace, so this proportion is a peak-height ratio, not an allele fraction'
    return {'genotype': g, 'basis_alt_prop': round(float(alt), 4),
            'confidence': conf, 'reason': reason}


# ---------------------------------------------------------------------------
# One sample, end to end
# ---------------------------------------------------------------------------

def analyze_sample(ab1, start, end, aa_first=None, reverse_complement=False,
                   residue=None, wt_aa=None, alt_aa=None, control_ab1=None,
                   outdir=None, variant=None, make_logos=True,
                   bystander_threshold=0.05, min_mean_phred=30, allow_indel=False):
    """Quality-gate, quantify, optionally plot. Returns (summary_dict, aa_df)."""
    variant = variant or os.path.splitext(os.path.basename(ab1))[0]

    qc = tr.assert_interpretable(ab1, start, end,
                                 min_mean_phred=min_mean_phred,
                                 allow_indel=allow_indel)

    nt_df = tr.base_proportions(ab1, start, end, reverse_complement)
    aa_df = tr.codon_aa_proportions(nt_df)
    mixed = tr.haplotype_warning(nt_df)

    if aa_first is None:
        aa_first = 1
    aa_lab = tr.label_aa_frame(aa_df, aa_first)

    control_lab = None
    if control_ab1:
        c_nt = tr.base_proportions(control_ab1, start, end, reverse_complement)
        control_lab = tr.label_aa_frame(tr.codon_aa_proportions(c_nt), aa_first)

    summary = {
        'variant': variant,
        'ab1': os.path.basename(ab1),
        'start': start,
        'end': end,
        'reverse_complement': reverse_complement,
        'called_bases_in_window': qc['called_bases'],
        'mean_phred': round(qc['mean_phred'], 1),
        'mean_top_base_prop': round(qc['mean_top_base_prop'], 3),
        'suspect_indel': qc['suspect_indel'],
        'multi_mixed_codons': ';'.join(str(m['codon_index']) for m in mixed),
        'control_ab1': os.path.basename(control_ab1) if control_ab1 else '',
    }

    if residue is not None:
        call = call_residue(aa_lab, residue, wt_aa, alt_aa,
                            bystander_threshold, control_lab)
        gt = genotype_call(call)
        summary.update({
            'residue': call['residue'],
            'wt_aa': call['wt_aa'],
            'alt_aa': call['alt_aa'],
            'top_aa': call['top_aa'],
            'wt_prop': call['wt_prop'],
            'alt_prop': call['alt_prop'],
            'alt_prop_ctrl_corrected': call['alt_prop_ctrl_corrected'],
            'bystanders': ';'.join(k + '=' + str(v) for k, v in call['bystanders'].items()),
            'genotype': gt['genotype'],
            'confidence': gt['confidence'],
            'call_note': gt['reason'],
        })

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        aa_lab.round(4).to_csv(os.path.join(outdir, variant + '_aa_proportions.csv'))
        nt_df.round(4).to_csv(os.path.join(outdir, variant + '_nt_proportions.csv'),
                              index=False)
        if make_logos and sanger_weblogo is not None:
            aa_last = aa_first + len(aa_df) - 1
            rng = (aa_first, aa_last)
            sanger_weblogo(ab1, start, end, 'nc', reverse_complement, rng,
                           os.path.join(outdir, variant + '_nt_logo.pdf'))
            sanger_weblogo(ab1, start, end, 'aa', reverse_complement, rng,
                           os.path.join(outdir, variant + '_aa_logo.pdf'))

    return summary, aa_lab


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _truthy(v):
    return str(v).strip().lower() in ('1', 'true', 't', 'yes', 'y')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--table', help='variant table CSV (batches a whole plate)')
    p.add_argument('--ab1')
    p.add_argument('--start', type=int)
    p.add_argument('--end', type=int)
    p.add_argument('--aa-first', type=int)
    p.add_argument('--aa-last', type=int,
                   help='checked against the codon count, as the notebook did')
    p.add_argument('--rc', action='store_true',
                   help='reverse-complement the window. For a minus-strand gene '
                        'such as TP53 this goes with the FW primer read -- verify '
                        'against sequence, not the filename')
    p.add_argument('--residue', type=int)
    p.add_argument('--wt')
    p.add_argument('--alt')
    p.add_argument('--control', help='matched unedited .ab1, same window/frame')
    p.add_argument('--bystander-threshold', type=float, default=0.05)
    p.add_argument('--min-mean-phred', type=float, default=30)
    p.add_argument('--allow-indel', action='store_true',
                   help='downgrade the indel refusal to a flag (logo only, do not quote a %%)')
    p.add_argument('--no-logos', action='store_true')
    p.add_argument('-o', '--outdir', default='genotyping_out')
    args = p.parse_args()

    rows = []
    if args.table:
        tab = pd.read_csv(args.table)
        required = ['ab1', 'start', 'end']
        missing = [c for c in required if c not in tab.columns]
        if missing:
            raise SystemExit("Variant table is missing column(s): " + ', '.join(missing))
        for _, r in tab.iterrows():
            aa_first = r.get('aa_first')
            aa_last = r.get('aa_last')
            start = int(r['start'])
            end = int(r['end'])
            if aa_first == aa_first and aa_last == aa_last:
                n_codons = (end - start + 1) // 3
                if int(aa_last) - int(aa_first) + 1 != n_codons:
                    raise SystemExit(
                        "Row " + str(r.get('variant', r['ab1'])) + ": aa_first/aa_last imply " +
                        str(int(aa_last) - int(aa_first) + 1) + " residues but the window has " +
                        str(n_codons) + " codons."
                    )
            summary, _ = analyze_sample(
                ab1=r['ab1'], start=start, end=end,
                aa_first=int(aa_first) if aa_first == aa_first else None,
                reverse_complement=_truthy(r.get('reverse_complement', False)),
                residue=int(r['residue']) if r.get('residue') == r.get('residue') and
                        pd.notna(r.get('residue')) else None,
                wt_aa=r.get('wt_aa') if pd.notna(r.get('wt_aa')) else None,
                alt_aa=r.get('alt_aa') if pd.notna(r.get('alt_aa')) else None,
                control_ab1=r.get('control_ab1') if pd.notna(r.get('control_ab1')) else None,
                outdir=args.outdir, variant=r.get('variant'),
                make_logos=not args.no_logos,
                bystander_threshold=args.bystander_threshold,
                min_mean_phred=args.min_mean_phred, allow_indel=args.allow_indel)
            rows.append(summary)
    else:
        if not (args.ab1 and args.start and args.end):
            raise SystemExit("Give either --table, or --ab1 with --start and --end.")
        if args.aa_first and args.aa_last:
            n_codons = (args.end - args.start + 1) // 3
            if args.aa_last - args.aa_first + 1 != n_codons:
                raise SystemExit(
                    "--aa-first/--aa-last imply " + str(args.aa_last - args.aa_first + 1) +
                    " residues but the window has " + str(n_codons) + " codons."
                )
        summary, aa_lab = analyze_sample(
            ab1=args.ab1, start=args.start, end=args.end, aa_first=args.aa_first,
            reverse_complement=args.rc, residue=args.residue,
            wt_aa=args.wt, alt_aa=args.alt, control_ab1=args.control,
            outdir=args.outdir, make_logos=not args.no_logos,
            bystander_threshold=args.bystander_threshold,
            min_mean_phred=args.min_mean_phred, allow_indel=args.allow_indel)
        rows.append(summary)
        print(aa_lab.round(3).to_string())

    out = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    calls_path = os.path.join(args.outdir, 'calls.csv')
    out.to_csv(calls_path, index=False)
    print(out.to_string(index=False))
    print('\nwrote ' + calls_path)
    if 'confidence' in out.columns and (out['confidence'] == 'uncalibrated').any():
        print('\nNOTE: rows marked "uncalibrated" had no matched WT control trace. '
              'Those proportions are peak-height ratios, not allele fractions. '
              'Pass --control (or a control_ab1 column) before quoting a percentage.')


if __name__ == '__main__':
    main()
