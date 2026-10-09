"""Verify that designed silent bystanders really are synonymous.

Lifted verbatim from early-leukemia-evolution/ch_aml_library_design.ipynb
cell 26: cds_cache, translate_raw, translate_cds, protein_after_edit and the
per-mutation audit loop. Do not edit here without editing there -- these two
must not drift.

Added as wrappers, not edits: argparse, the checkpoint loader, the
has_silent_bystander re-cast, and the seq_start precondition check.

HOW IT WORKS, and why not "at most one amino acid changed vs reference":
demanding that only holds for missense variants, and 20 of the 44 CH/AML
variants are frameshift, nonsense or in-frame indels that legitimately rewrite
many residues. Each bystander design is compared against the protein the PLAIN
designs for the SAME mutation agree on -- same intended edit, so any difference
is the bystander's doing, for every variant class.

PRECONDITION: seq_start exists only on the 'cBioPortal' input path (it is set
by prime.df_formatter). A library designed from 'WT_ALT' or 'PrimeDesign' input
cannot be audited by this script, because nothing maps a design back to genomic
coordinates. It fails by name rather than guessing.
"""

import argparse
import os
from collections import Counter

import Bio.Seq
import pandas as pd

import prep_input

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def cds_cache(chrom_seq, cds_blocks, strand):
    """CDS bases in transcription order, plus genomic position -> offset.

    Built ONCE per transcript. The earlier version copied the whole chromosome
    into a Python list for every pegRNA (list(str(chrom_seq)) on ~182 Mb), which
    is what made this cell take hours; only the CDS bases can affect the
    translation, so cache those and apply each RTT as a small override.
    """
    from pegg import bystander
    pos = bystander.cds_positions(cds_blocks, strand)
    raw = ''.join(str(chrom_seq[p - 1]) for p in pos).upper()
    return pos, raw, {p: k for k, p in enumerate(pos)}


def translate_raw(raw, strand):
    """Translate CDS bases already in transcription order."""
    return str(Bio.Seq.Seq(raw.translate(_COMP) if strand == '-' else raw).translate())


def translate_cds(chrom_seq, cds_blocks, strand):
    _, raw, _ = cds_cache(chrom_seq, cds_blocks, strand)
    return translate_raw(raw, strand)


def protein_after_edit(row, cache, strand):
    """Apply this pegRNA's RTT to the CDS and translate.

    The RTT replaces the genomic span from the nick to the end of the RHA. For a
    substitution that span is the same length as the RTT, but for an insertion or
    deletion it is not -- so the RTT is SPLICED in (replacing a span of a
    different length), not written base-over-base.

    Every base of the RTT is mapped to the genomic coordinate it occupies and
    substituted individually. Doing it per base rather than as one block is what
    lets an RTT that overhangs an exon edge still contribute the part of itself
    that is coding; blocks that did not fit entirely inside the CDS used to be
    discarded whole, which quietly removed the intended edit as well.

    Returns None when the RTT's arithmetic does not close, rather than guessing.
    """
    # RTT is stored reverse-complemented; undo that, then orient to the FORWARD
    # genome strand (the PAM strand says which way the RTT was written).
    rtt_fwd = str(Bio.Seq.Seq(row['RTT']).reverse_complement())
    if row['PAM_strand'] == '-':
        rtt_fwd = str(Bio.Seq.Seq(rtt_fwd).reverse_complement())

    ref = '' if pd.isna(row['REF']) else str(row['REF'])
    alt = '' if pd.isna(row['ALT']) else str(row['ALT'])
    d, rha = int(row['Distance_to_nick']), int(row['RHA_size'])
    if len(rtt_fwd) != d + len(alt) + rha:
        return None

    # genomic span the RTT replaces: d bases of LHA, the REF allele, then the RHA
    edit_g = int(row['seq_start']) + len(str(row['left_context']))
    g0 = edit_g - d if row['PAM_strand'] == '+' else edit_g - rha
    span = d + len(ref) + rha

    pos, raw, index = cache

    # Map each base of rtt_fwd to the genomic coordinate it occupies, so the
    # splice can be done per base rather than as one all-or-nothing block.
    # The RTT is LHA + ALT + RHA in PAM-strand order; on the - strand the arms
    # are swapped relative to the forward genome, so the seam sits at RHA.
    # Inserted bases have no reference coordinate of their own and are carried
    # as None -- they are emitted with whichever CDS base precedes them.
    seam = d if row['PAM_strand'] == '+' else rha
    rtt_coords = []
    for i in range(len(rtt_fwd)):
        if i < seam:                                # left homology arm
            rtt_coords.append(g0 + i)
        elif i < seam + len(alt):
            # An ALT base sits on a REF base for as far as the REF allele goes;
            # only ALT beyond len(REF) is a true insertion with no coordinate of
            # its own. Marking the whole ALT as inserted broke every design whose
            # edit sits at the RTT's first base (Distance_to_nick == 0): the SNP's
            # own coordinate was dropped and the protein came out a residue short.
            j = i - seam
            rtt_coords.append(g0 + seam + j if j < len(ref) else None)
        else:                                       # right homology arm
            rtt_coords.append(g0 + seam + len(ref) + (i - seam - len(alt)))

    # Group the RTT's bases by the coordinate they occupy. Inserted bases (None)
    # have no coordinate of their own, so they are appended to the coordinate of
    # the base that PRECEDES them in transcript order -- which, on a - strand
    # transcript, is the base that FOLLOWS them on the forward genome. Getting
    # this wrong puts an in-frame insertion a few residues off.
    override = {}
    carry = ''
    for i, c in enumerate(rtt_coords):
        if c is None:
            carry += rtt_fwd[i].upper()
        else:
            override[c] = override.get(c, '') + carry + rtt_fwd[i].upper()
            carry = ''
    if carry:
        # trailing insertion: attach to the last real coordinate seen
        for c in reversed(rtt_coords):
            if c is not None:
                override[c] = override[c] + carry
                break

    # Reference bases the edit removes: covered by the RTT's span but carried by
    # none of its bases.
    covered = {c for c in rtt_coords if c is not None}
    deleted = set(range(g0 + seam, g0 + seam + len(ref))) - covered

    # Walk the CDS in transcription order, substituting base by base. Only CDS
    # positions contribute, so an RTT that runs past an exon edge now keeps the
    # part of itself that is coding -- the previous version threw the WHOLE RTT
    # away in that case, which silently dropped the intended edit and made those
    # rows translate as wild type. On Dnmt3a R767* that was 45% of the rows, and
    # it is what produced the SNP "non-silent" calls in the earlier audit.
    #
    # On a - strand transcript `pos` runs backwards along the chromosome and
    # translate_raw() complements the result, so a multi-base override (an
    # insertion) has to be reversed to come out in transcript order.
    out = []
    for p in pos:
        if p in deleted:
            continue
        if p in override:
            piece = override[p]
            out.append(piece if strand == '+' else piece[::-1])
        else:
            out.append(raw[index[p]])
    return translate_raw(''.join(out), strand)


def audit(peg_df, gene_cds, chrom_dict, tx_column, gene_column="Hugo_Symbol",
          audit_per_mut=200, audit_seed=0):
    """Per-mutation bystander audit. Lifted from cell 26's loop.

    gene_cds from add_cds_to_variants(tx_column=...) is keyed by TRANSCRIPT id,
    not gene name -- grouping by Hugo_Symbol would look up a key that does not
    exist and skip every gene without saying so.

    With pegRNAs_per_mut='All' the table is ~665k rows and the check is
    per-row, which would take ~3 h. Whether a bystander is silent is a property
    of the DESIGN, and the failure is systematic per mutation (an affected
    mutation fails at a high rate, not at 1-2%), so a bounded RANDOM sample per
    mutation answers the same question. Random, not .head(): rows arrive
    ordered by PAM site then RTT/PBS length, so the first N would all share one
    corner of the design space.
    """
    rows = []
    for (gene, mi), sub in peg_df.groupby([gene_column, 'mutation_idx']):
        n_pool = len(sub)
        if n_pool > audit_per_mut:
            # sample each design type separately so both stay represented
            sub = (sub.groupby('has_silent_bystander', group_keys=False)
                   .apply(lambda g: g.sample(min(len(g), audit_per_mut // 2),
                                             random_state=audit_seed)))
        tx = sub[tx_column].iloc[0]
        ann = gene_cds.get(tx)
        if ann is None or not ann['valid']:
            continue
        raw = str(sub['Chromosome'].iloc[0]).replace('chr', '')
        chrom_key = raw if raw in ('X', 'Y') else int(raw)
        cache = cds_cache(str(chrom_dict[chrom_key]).upper(),
                          ann['cds'], ann['strand'])

        plain = [protein_after_edit(r, cache, ann['strand'])
                 for _, r in sub[~sub['has_silent_bystander']].iterrows()]
        plain = [x for x in plain if x]
        if not plain:
            continue
        # The plain designs for one mutation all carry the same edit, so they
        # must all translate to the same protein. Any disagreement means the
        # splice is wrong for some geometry, which would silently poison the
        # baseline the bystander rows are compared against -- so it is
        # recorded, not just assumed.
        intended = Counter(plain).most_common(1)[0][0]
        n_plain_disagree = sum(1 for x in plain if x != intended)

        bad = [r['pegRNA_id'] for _, r in sub[sub['has_silent_bystander']].iterrows()
               if (lambda q: q is not None and q != intended)(
                   protein_after_edit(r, cache, ann['strand']))]

        ref, alt = sub['REF'].iloc[0], sub['ALT'].iloc[0]
        rn, an = pd.isna(ref), pd.isna(alt)
        vtype = ('INS' if rn and not an else 'DEL' if an and not rn
                 else 'SNP' if len(str(ref)) == len(str(alt)) else 'INDEL')
        rows.append(dict(gene=gene, mutation_idx=mi, variant_type=vtype,
                         REF=ref, ALT=alt, frameshift=vtype in ('INS', 'DEL'),
                         pool=n_pool, audited=len(sub),
                         plain=len(plain), plain_disagree=n_plain_disagree,
                         bys=int(sub['has_silent_bystander'].sum()),
                         non_silent=len(bad), ids=';'.join(map(str, bad[:3]))))
    return pd.DataFrame(rows)


def report(aud, audit_per_mut, audit_seed):
    """Cell 26's printed summary."""
    print(aud[aud['non_silent'] > 0].to_string(index=False)
          if (aud['non_silent'] > 0).any()
          else 'every bystander design is silent')
    print()
    by_type = aud.groupby('variant_type').agg(
        mutations=('gene', 'size'), bystander=('bys', 'sum'),
        non_silent=('non_silent', 'sum'))
    by_type['rate_%'] = (100 * by_type['non_silent'] / by_type['bystander']).round(1)
    print(by_type.to_string())
    print()
    print('audited %s of %s pegRNAs (<= %d sampled per mutation, seed %d)'
          % ('{:,}'.format(aud['audited'].sum()),
             '{:,}'.format(aud['pool'].sum()), audit_per_mut, audit_seed))
    print('%d of %d bystander designs are NOT silent (%.1f%%), over %d of %d mutations'
          % (aud['non_silent'].sum(), aud['bys'].sum(),
             100 * aud['non_silent'].sum() / max(aud['bys'].sum(), 1),
             (aud['non_silent'] > 0).sum(), len(aud)))
    print()
    print('A frameshift shifts the reading frame, so a codon PEGG picked as synonymous')
    print('in the reference frame need not be synonymous in the edited one -- which is')
    print('why the indel variants carry nearly all of these. Filter with:')
    print("   peg_df = peg_df[~peg_df['pegRNA_id'].isin(non_silent_ids)]")
    print('or drop bystander designs for the frameshift mutations before ordering.')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("designs", help="design table from design.py (.csv or .csv.gz)")
    ap.add_argument("variants", help="the prepared variant table it was designed "
                                     "from, carrying start_end_cds / "
                                     "transcript_strand")
    ap.add_argument("-o", "--out", required=True, help="audit CSV to write")
    ap.add_argument("--build", required=True, choices=sorted(prep_input.GENOMES))
    ap.add_argument("--tx-column", required=True,
                    help="the transcript column the frame was annotated with "
                         "(e.g. tx_id_m). gene_cds is keyed by this, not by "
                         "gene name.")
    ap.add_argument("--gene-column", default="Hugo_Symbol")
    ap.add_argument("--per-mut", type=int, default=200,
                    help="max pegRNAs audited per mutation (random sample)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    peg_df = pd.read_csv(args.designs, low_memory=False)
    if 'seq_start' not in peg_df.columns:
        raise SystemExit(
            "this design table has no seq_start column, so no design can be "
            "mapped back to genomic coordinates. seq_start is set only on the "
            "'cBioPortal' input path (prime.df_formatter); a library designed "
            "from 'WT_ALT' or 'PrimeDesign' input cannot be audited here. "
            "Redesign from a MAF-format table if you need this check.")
    if 'has_silent_bystander' not in peg_df.columns:
        raise SystemExit("no has_silent_bystander column -- this table was "
                         "designed without silent_bystander=True, so there is "
                         "nothing to audit.")

    import design as _design
    peg_df = _design.recast_bystander_flag(peg_df)
    if args.tx_column not in peg_df.columns:
        raise SystemExit("%r is not in the design table; check df_formatter "
                         "kept it" % args.tx_column)

    variants = pd.read_csv(args.variants)
    chrom_dict, _ = prep_input.load_genome(args.build)
    db = prep_input.load_annotation(args.build)
    _v, gene_cds, _k = prep_input.annotate_frame(
        variants, db, tx_column=args.tx_column, gene_column=args.gene_column)

    aud = audit(peg_df, gene_cds, chrom_dict, args.tx_column,
                gene_column=args.gene_column,
                audit_per_mut=args.per_mut, audit_seed=args.seed)
    if not len(aud):
        raise SystemExit("nothing auditable: no mutation had both a valid CDS "
                         "and at least one plain design to compare against.")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    aud.to_csv(args.out, index=False)
    report(aud, args.per_mut, args.seed)
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
