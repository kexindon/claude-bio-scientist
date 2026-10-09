"""Build a matched synonymous-codon control variant table.

Lifted verbatim from early-leukemia-evolution/ch_aml_library_design.ipynb
cell 29 (silent_controls and its two assertions). Do not edit here without
editing there -- these two must not drift.

Added as wrappers, not edits: argparse, the LAB_DATA_ROOT genome resolution,
and writing the skipped list out as a table instead of only printing it.

WHY A SECOND LIBRARY and not pegg's library.neutral_substitutions(): these
controls are matched ONE PER VARIANT, at the same codon, so each real design
has a synonymous twin at the same locus with the same geometry. Design them
through a SECOND prime.run() call with IDENTICAL parameters (including
silent_bystander), so the controls span the same efficiency range as the real
designs -- see design.py and SKILL.md.
"""

import argparse
import os

import Bio.Seq
import numpy as np
import pandas as pd

import prep_input


def silent_controls(design_input, chrom_dict, seed=0, single_aa_only=True):
    """Synonymous-codon controls, one per variant.

    single_aa_only restricts them to variants that change exactly ONE amino acid
    -- missense and nonsense SNPs/ONPs within a single codon. A frameshift or a
    multi-residue indel has no single 'corresponding' codon to swap: its control
    would be a synonymous change at one arbitrary residue while the real variant
    rewrites the protein from that point on, so the two are not comparable as an
    activity control. Those variants are skipped and reported.
    """
    from pegg import bystander

    rng = np.random.default_rng(seed)
    rows, skipped = [], []
    for _, val in design_input.iterrows():
        if single_aa_only:
            ref = '' if pd.isna(val['Reference_Allele']) else str(val['Reference_Allele']).replace('-', '')
            alt = '' if pd.isna(val['Tumor_Seq_Allele2']) else str(val['Tumor_Seq_Allele2']).replace('-', '')
            # an indel (unequal lengths) shifts or rewrites the reading frame;
            # a substitution longer than a codon spans more than one residue
            if len(ref) != len(alt):
                skipped.append((val['Hugo_Symbol'], 'indel -- not a single-aa change')); continue
            if len(ref) > 3:
                skipped.append((val['Hugo_Symbol'], 'spans >1 codon')); continue
        blocks = val.get('start_end_cds')
        strand = val.get('transcript_strand')
        if blocks is None or (isinstance(blocks, float) and np.isnan(blocks)) or not strand:
            skipped.append((val['Hugo_Symbol'], 'no CDS annotation')); continue
        seq = chrom_dict.get(val['Chromosome'])
        if seq is None:
            skipped.append((val['Hugo_Symbol'], 'chromosome not loaded')); continue

        fmap = bystander.cds_frame_map(blocks, strand)
        pos = int(val['Start_Position'])
        if pos not in fmap:
            skipped.append((val['Hugo_Symbol'], 'variant not in CDS')); continue
        codon_idx, _ = fmap[pos]

        codons = bystander.cds_codons(blocks, strand)
        if codon_idx >= len(codons):
            skipped.append((val['Hugo_Symbol'], 'codon index out of range')); continue
        cpos = codons[codon_idx]

        bases = [str(seq[p - 1]).upper() for p in cpos]
        if strand == '-':
            bases = [str(Bio.Seq.Seq(b).complement()) for b in bases]
        codon = ''.join(bases)

        syn = bystander.synonymous_codons(codon)
        if not syn:
            skipped.append((val['Hugo_Symbol'], 'no synonymous codon (Met/Trp)')); continue

        fewest = min(sum(a != b for a, b in zip(c, codon)) for c in syn)
        cands = [c for c in syn if sum(a != b for a, b in zip(c, codon)) == fewest]
        newcodon = cands[int(rng.integers(len(cands)))]

        fwd = sorted(cpos)
        ref_fwd = ''.join(str(seq[p - 1]).upper() for p in fwd)
        nb_ = list(newcodon)
        if strand == '-':
            nb_ = [str(Bio.Seq.Seq(b).complement()) for b in nb_][::-1]
        alt_fwd = ''.join(nb_)

        r = val.copy()
        r['Start_Position'], r['End_Position'] = fwd[0], fwd[-1]
        r['Reference_Allele'], r['Tumor_Seq_Allele2'] = ref_fwd, alt_fwd
        r['Variant_Type'] = 'ONP' if len(ref_fwd) == len(alt_fwd) else 'INDEL'
        r['silent_control_for'] = '%s:%s>%s@%d' % (
            val['Hugo_Symbol'], val['Reference_Allele'], val['Tumor_Seq_Allele2'], pos)
        r['silent_codon_change'] = '%s>%s' % (codon, newcodon)
        r['silent_aa'] = bystander.translate_codon(codon)
        rows.append(r)
    return pd.DataFrame(rows).reset_index(drop=True), skipped


def verify(silent_input, chrom_dict):
    """Cell 29's two assertions. Both are cheap and both catch a real failure.

    The protein must be unchanged -- that is the whole point of the control --
    and REF must match the genome, or pegg designs against the wrong sequence
    and the resulting pegRNAs edit nothing.
    """
    from pegg import bystander

    bad = [r['silent_codon_change'] for _, r in silent_input.iterrows()
           if bystander.translate_codon(r['silent_codon_change'].split('>')[0])
           != bystander.translate_codon(r['silent_codon_change'].split('>')[1])]
    assert not bad, 'codon change is not synonymous: %s' % bad

    for _, r in silent_input.iterrows():
        g = ''.join(str(chrom_dict[r['Chromosome']][p - 1]).upper()
                    for p in range(int(r['Start_Position']), int(r['End_Position']) + 1))
        assert g == r['Reference_Allele'], (r['Hugo_Symbol'], g, r['Reference_Allele'])
    print('checked: all synonymous, all REF alleles match the genome')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("variants", help="prepared variant table from prep_input.py "
                                     "(needs start_end_cds and transcript_strand)")
    ap.add_argument("-o", "--out", required=True,
                    help="control variant table to write; feed it to design.py "
                         "with the SAME parameters as the real library")
    ap.add_argument("--build", required=True, choices=sorted(prep_input.GENOMES))
    ap.add_argument("--seed", type=int, default=0,
                    help="ties between equally-close synonymous codons are "
                         "broken at random; fix this or the control library is "
                         "not reproducible")
    ap.add_argument("--all-variant-types", action="store_true",
                    help="do NOT restrict to single-amino-acid changes. Read "
                         "the docstring first: an indel's control is not "
                         "comparable to it as an activity control.")
    args = ap.parse_args()

    design_input = pd.read_csv(args.variants)
    for col in ("start_end_cds", "transcript_strand"):
        if col not in design_input.columns:
            raise SystemExit(
                "input has no %r column. Run prep_input.py with --tx-column "
                "first: these controls need the reading frame to find the "
                "codon to swap." % col)
    design_input = prep_input.normalize_chromosome(design_input)
    chrom_dict, _ = prep_input.load_genome(args.build)
    prep_input.assert_chroms_present(design_input, chrom_dict)

    silent_input, skipped = silent_controls(
        design_input, chrom_dict, seed=args.seed,
        single_aa_only=not args.all_variant_types)
    print('silent controls: %d of %d variants (single-amino-acid changes only)'
          % (len(silent_input), len(design_input)))
    for g, why in skipped:
        print('   skipped %-8s %s' % (g, why))

    verify(silent_input, chrom_dict)

    silent_input.to_csv(args.out, index=False)
    print("wrote %s (%d control variants)" % (args.out, len(silent_input)))
    if skipped:
        spath = os.path.splitext(args.out)[0] + "_skipped.csv"
        pd.DataFrame(skipped, columns=["gene", "reason"]).to_csv(spath, index=False)
        print("wrote %s (%d variants with no comparable control)"
              % (spath, len(skipped)))


if __name__ == "__main__":
    main()
