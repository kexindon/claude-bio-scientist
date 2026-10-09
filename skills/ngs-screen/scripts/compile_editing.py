#!/usr/bin/env python3
"""STEP7 + STEP8, local: aggregated CRISPResso CSVs -> editing tables -> count matrix.

This is where the pipeline ends. LFC, z-scores, empirical p-values, FDR and hit
calling are the screen_stats skill, which reads exactly the two outputs written
here.

Lifted verbatim from:
  fsrlab-sensor-analysis/NGS-scripts/STEP7_crispresso_compiler/compiling_crispresso.ipynb
    cell 5  -> crispresso_compiler   (column list copied character for character;
                screen_stats' attach_editing_info reads these names)
    cell 7  -> MLE_merge_auto
    cell 9  -> merge_NGS_replicates
    (cell 11's MLE_merge is the legacy fixed-3-replicate version. Not lifted:
     MLE_merge_auto supersedes it and handles groups of any size. If you are
     reproducing a pre-2025 analysis that used it, note that it accepts exactly
     3 replicates and silently indexes k[0..2].)
  fsrlab-sensor-analysis/NGS-scripts/STEP8_counts_matrix/counts.ipynb
    cell 12 -> count_matrix
Do not edit the column lists or the percentage formulas here without editing
there -- these must not drift.

The notebooks drove these functions from hand-written parallel Python lists --
48 raw D25-xxxxx-N-7331E ids in one list and 48 friendly ids in another, kept in
register by eye. Here that mapping is a two-column TSV (--samples), so it is
reviewable and diffable.

ONE REAL BUG FIXED, in count_matrix:

    new_sample_ids = np.unique([x.rsplit('-',1)[0] for x in samp_id])

np.unique sorts alphabetically, while samps_to_merge is built in file order.
The loop then does join_out[new_sample_ids[i]] = join_out[samps_to_merge[i]].sum(...),
pairing the i-th lane pair with the i-th ALPHABETICAL name. With the real sample
set (lib, input, d5, d15-rep1..5, spleen1..5, bm1..5, men1..5) the two orders
differ, so lane-merged counts were written under the wrong sample column.
STEP7's merge_NGS_replicates already does this correctly with dict.fromkeys;
that is what is used here. Anything built with the notebook's count_matrix
before this fix should be rebuilt.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

# MLE_merge_auto is lifted verbatim, and `out = merged[col_order]` is a slice, so
# every percentage assignment raises SettingWithCopyWarning -- seven per group,
# which buries the real messages. The assignments are correct (out is the only
# reference and is returned), so the warning is silenced here rather than by
# editing lifted code that must stay comparable to the notebook.
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)


# ---------------------------------------------------------------------------
# STEP7 -- lifted verbatim
# ---------------------------------------------------------------------------

def crispresso_compiler(samp, sample_id, fp):
    """
    Takes in list of sample names (see cell above)
    Returns compiled dictionary of dataframes containing editing information for each sensor
    - key for each is the sample_id name

    Includes:
    1. corr_perc = pure correct editing
    2. target_base_edit_perc = target base editing perc (including edits with bystander editing)
    3. wt_perc
    4. byproduct information (indels, substitutions, ambiguous)
    """

    df_edits = []
    for k in samp:
        concated = pd.read_csv(os.path.join(fp, str(k) + '_crispresso_aggregated.csv'))
        concated = concated.fillna(0)

        # go through all of the samples and do it step by step
        sample_ids = []
        rii = []
        raaa = []
        r_lowqual = []
        ra_hdr = []
        ra_wt = []
        r_unaligned = []
        no_edit = []
        correct = []
        target_base_editing = []
        byproduct_all = []
        byproduct_indel = []
        byproduct_sub = []
        byproduct_ambig = []
        for i in np.unique(concated['Guide_ID']):
            sample_ids.append(i)

            subset = concated[concated['Guide_ID'] == i]

            wt = subset[subset['Amplicon'] == 'Reference']
            edit = subset[subset['Amplicon'] == 'HDR']

            rii1 = wt['Reads_in_input'].values[0]
            r1 = wt['Reads_aligned_all_amplicons'].values[0]
            rii.append(rii1)
            raaa.append(r1)
            r_lowqual.append(rii1 - r1)

            r2 = edit['Reads_aligned'].values[0]
            r3 = wt['Reads_aligned'].values[0]
            ra_hdr.append(r2)
            ra_wt.append(r3)
            r_unaligned.append(r1 - (r2 + r3))

            no_edit.append(wt['Unmodified'].values[0])
            correct.append(edit['Unmodified'].values[0])
            target_base_editing.append(edit['Modified'].values[0] + edit['Unmodified'].values[0])

            byprod_all = wt['Modified'].values[0] + edit['Modified'].values[0] + (r1 - (r2 + r3))  # add unaligned reads
            sub_all = wt['Only Substitutions'].values[0] + edit['Only Substitutions'].values[0]
            indel_all = wt['Only Deletions'].values[0] + wt['Only Insertions'].values[0] + wt['Insertions and Deletions'].values[0] + edit['Only Deletions'].values[0] + edit['Only Insertions'].values[0] + edit['Insertions and Deletions'].values[0]
            ambig_all = byprod_all - sub_all - indel_all

            byproduct_all.append(byprod_all)
            byproduct_indel.append(indel_all)
            byproduct_sub.append(sub_all)
            byproduct_ambig.append(ambig_all)

        cols = ["Guide_ID", "Reads_in_input", "Reads_lowqual", "Reads_aligned_all_amplicons", "Reads_aligned_WT", "Reads_aligned_HDR", "Reads_unaligned", "WT", "correct_edit", "target_base_edit", "byproduct_all", "byproduct_INDEL", "byproduct_sub", "byproduct_ambiguous"]
        col_vals = [sample_ids, rii, r_lowqual, raaa, ra_wt, ra_hdr, r_unaligned, no_edit, correct, target_base_editing, byproduct_all, byproduct_indel, byproduct_sub, byproduct_ambig]

        out = pd.DataFrame(dict(zip(cols, col_vals)))
        out['corr_perc'] = 100 * (out['correct_edit'] / out['Reads_aligned_all_amplicons'])
        out['target_base_edit_perc'] = 100 * (out['target_base_edit'] / out['Reads_aligned_all_amplicons'])
        out['WT_perc'] = 100 * (out['WT'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_all_perc'] = 100 * (out['byproduct_all'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_INDEL_perc'] = 100 * (out['byproduct_INDEL'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_sub_perc'] = 100 * (out['byproduct_sub'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_ambiguous_perc'] = 100 * (out['byproduct_ambiguous'] / out['Reads_aligned_all_amplicons'])
        out = out.fillna(0)
        out = out.sort_index()
        df_edits.append(out)

    edit_dict = dict(zip(sample_id, df_edits))

    return edit_dict


def MLE_merge_auto(edit_dict, samps_to_merge, new_name):
    """
    Combine replicates to generate best estimate of sensor editing.
    Automatically handles any number of replicates per group (>=1).
    """
    cols_to_add = [
        'Reads_in_input', 'Reads_lowqual', 'Reads_aligned_all_amplicons',
        'Reads_aligned_WT', 'Reads_aligned_HDR', 'Reads_unaligned', 'WT',
        'correct_edit', 'target_base_edit', 'byproduct_all',
        'byproduct_INDEL', 'byproduct_sub', 'byproduct_ambiguous',
    ]

    comb_holder = []

    for group in samps_to_merge:
        dfs = [edit_dict[sample] for sample in group]
        guide_check = [list(df['Guide_ID']) for df in dfs]
        assert all(g == guide_check[0] for g in guide_check), 'Guide_ID dont match'

        merged = dfs[0][cols_to_add].copy()
        for df in dfs[1:]:
            merged += df[cols_to_add]

        merged['Guide_ID'] = dfs[0]['Guide_ID']
        col_order = ['Guide_ID'] + cols_to_add
        out = merged[col_order]

        out['corr_perc'] = 100 * (out['correct_edit'] / out['Reads_aligned_all_amplicons'])
        out['target_base_edit_perc'] = 100 * (out['target_base_edit'] / out['Reads_aligned_all_amplicons'])
        out['WT_perc'] = 100 * (out['WT'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_all_perc'] = 100 * (out['byproduct_all'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_INDEL_perc'] = 100 * (out['byproduct_INDEL'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_sub_perc'] = 100 * (out['byproduct_sub'] / out['Reads_aligned_all_amplicons'])
        out['byproduct_ambiguous_perc'] = 100 * (out['byproduct_ambiguous'] / out['Reads_aligned_all_amplicons'])
        out = out.fillna(0)
        out = out.sort_index()   # keep original row order
        comb_holder.append(out)

    return dict(zip(new_name, comb_holder))


def merge_NGS_replicates(sample_ids, edit_dict):
    samps_to_merge = [
        [sample_ids[i], sample_ids[i + 1]]
        for i in range(0, len(sample_ids), 2)
    ]
    new_sample_ids = list(dict.fromkeys(x.rsplit("-", 1)[0] for x in sample_ids))
    comb_dict = MLE_merge_auto(edit_dict, samps_to_merge, new_sample_ids)
    return comb_dict


# ---------------------------------------------------------------------------
# STEP8 -- lifted, with the np.unique ordering bug fixed
# ---------------------------------------------------------------------------

def count_matrix(dictionary, samp, samp_id, fp, chosen_count_column,
                 guide_gene, class_gene, merge_lanes=True):
    df_holder = []
    for i, k in enumerate(samp):
        aa = pd.read_csv(os.path.join(fp, str(k) + '_count_df.csv'))
        aa = aa[['Guide_ID', chosen_count_column]]
        s = dictionary[k]
        aa = aa.rename(columns={chosen_count_column: str(s)})

        if i > 0:
            aa = aa.rename(columns={'Guide_ID': 'Guide_ID_x'})
        df_holder.append(aa)

    join_out = pd.concat(df_holder, axis=1)
    sample_id_expanded_to_add = []
    for x in samp_id:
        if x in join_out.columns:
            sample_id_expanded_to_add.append(x)
    join_out = join_out[['Guide_ID'] + sample_id_expanded_to_add]
    join_out = join_out.loc[:, ~join_out.columns.duplicated()]

    genes = []
    for i, val in join_out.iterrows():
        k = val['Guide_ID']
        gene = guide_gene[k]
        # make all the neutral guides have the same name
        class_ = class_gene[k]
        if class_ in ['non-targeting control', 'safe-targeting control']:
            gene = 'controlsgRNA'
        genes.append(gene)
    join_out['gene'] = genes
    join_out = join_out.rename(columns={'Guide_ID': 'sgRNA'})
    cols = ['sgRNA', 'gene'] + [c for c in join_out.columns if c not in ['sgRNA', 'gene']]
    join_out = join_out[cols]

    if merge_lanes:
        samps_to_merge = [
            [samp_id[i], samp_id[i + 1]]
            for i in range(0, len(samp_id), 2)
        ]
        # FIX: was np.unique(...), which sorts alphabetically and therefore
        # mislabelled the merged columns. dict.fromkeys keeps file order, which
        # is what samps_to_merge is in -- same as merge_NGS_replicates.
        new_sample_ids = list(dict.fromkeys(x.rsplit("-", 1)[0] for x in samp_id))
        if len(new_sample_ids) != len(samps_to_merge):
            raise SystemExit(
                'lane merge: ' + str(len(samp_id)) + ' lanes collapse to '
                + str(len(new_sample_ids)) + ' names but pair into '
                + str(len(samps_to_merge)) + ' groups. Every sample must have exactly two '
                'lanes named <sample>-1 and <sample>-2, consecutive in the samples file. '
                'Fix the samples file, or pass --no-merge-lanes and merge explicitly.')
        for i, x in enumerate(samps_to_merge):
            join_out[new_sample_ids[i]] = join_out[x].sum(axis=1)
            join_out = join_out.drop(columns=x)
    return join_out


# ---------------------------------------------------------------------------
# samples file + CLI -- wrappers around the lifted code
# ---------------------------------------------------------------------------

def read_samples(path):
    """Two columns, tab- or comma-separated, optional header: raw_folder, sample_id.

    raw_folder is the STEP4 folder_name, e.g.
    D25-13245-1-7331E_guide_split_BARCODES. sample_id is the friendly name with
    a lane suffix, e.g. input-1. Order matters: lanes of the same sample must be
    adjacent and named <sample>-1, <sample>-2.
    """
    sep = '\t' if path.endswith('.tsv') or path.endswith('.txt') else ','
    df = pd.read_csv(path, sep=sep, header=None, comment='#',
                     names=['raw_folder', 'sample_id'], dtype=str)
    if str(df.iloc[0]['raw_folder']).lower() in ('raw_folder', 'folder_name', 'folder'):
        df = df.iloc[1:].reset_index(drop=True)
    df['raw_folder'] = df['raw_folder'].str.strip()
    df['sample_id'] = df['sample_id'].str.strip()

    dup = df['sample_id'][df['sample_id'].duplicated()].tolist()
    if dup:
        raise SystemExit('duplicate sample_id in ' + path + ': ' + ', '.join(dup))

    # The merge pairs consecutive rows blindly. If any sample has one lane
    # instead of two the pairing shifts by one and merges reads from two
    # different biological samples with no error at all.
    stems = {}
    for sid in df['sample_id']:
        stems.setdefault(sid.rsplit('-', 1)[0], []).append(sid)
    odd = dict((k, v) for k, v in stems.items() if len(v) != 2)
    return df, odd


def cmd_samples(args):
    """Generate a samples file skeleton from a config file, for editing."""
    config = pd.read_csv(args.config, sep=' ')
    print('# raw_folder\tsample_id')
    print('# Fill in sample_id: a friendly name plus a lane suffix, e.g. input-1, input-2,')
    print('# d15-rep1-1, spleen1-2. Lanes of one sample must be adjacent and named -1/-2.')
    for folder in config['folder_name']:
        print(str(folder).lstrip('./') + '\t')
    return 0


def cmd_compile(args):
    samples, odd = read_samples(args.samples)
    if odd and not args.no_merge_lanes:
        for stem in sorted(odd):
            sys.stderr.write('WARNING: ' + stem + ' has lanes ' + repr(odd[stem])
                             + ', expected exactly two\n')
        raise SystemExit(
            'Lane merging pairs consecutive rows, so a sample with one lane shifts the '
            'pairing and merges two different biological samples silently. Fix the samples '
            'file, or pass --no-merge-lanes to write per-lane tables and merge explicitly '
            'with MLE_merge_auto.')

    edit_dict = crispresso_compiler(list(samples['raw_folder']),
                                    list(samples['sample_id']),
                                    args.crispresso_dir)

    if args.no_merge_lanes:
        out_dict = edit_dict
    else:
        out_dict = merge_NGS_replicates(list(samples['sample_id']), edit_dict)

    os.makedirs(args.out_dir, exist_ok=True)
    for sample_name, df in out_dict.items():
        path = os.path.join(args.out_dir, str(sample_name) + '_compact_unfiltered.csv')
        df.to_csv(path, index=False)
    print(str(len(out_dict)) + ' tables -> ' + args.out_dir + '/<sample>_compact_unfiltered.csv')

    neg = []
    for sample_name, df in out_dict.items():
        n = int((df['byproduct_ambiguous'] < 0).sum())
        if n:
            neg.append(sample_name + ' (' + str(n) + ')')
    if neg:
        sys.stderr.write(
            'NOTE: byproduct_ambiguous is negative for some guides in: ' + ', '.join(neg)
            + '. It is a residual (byproduct_all - sub_all - indel_all) and goes negative '
              'on low-read guides where CRISPResso double-counts a class. The lifted code '
              'only fillna(0)s, it does not clip. Filter on '
              'Reads_aligned_all_amplicons before interpreting byproduct splits.\n')
    return 0


def cmd_mle(args):
    """Merge biological replicates. Groups come from a spec file: name<TAB>s1,s2,s3."""
    loaded = {}
    for fname in sorted(os.listdir(args.compact_dir)):
        if fname.endswith('_compact_unfiltered.csv'):
            name = fname.replace('_compact_unfiltered.csv', '')
            loaded[name] = pd.read_csv(os.path.join(args.compact_dir, fname))
    if not loaded:
        raise SystemExit('no *_compact_unfiltered.csv under ' + args.compact_dir)

    groups = []
    names = []
    with open(args.groups) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.replace('\t', ',').split(',')
            names.append(parts[0].strip())
            members = [x.strip() for x in parts[1:] if x.strip()]
            missing = [m for m in members if m not in loaded]
            if missing:
                raise SystemExit('group ' + parts[0] + ' names samples that have no compact '
                                 'table: ' + ', '.join(missing))
            groups.append(members)

    out = MLE_merge_auto(loaded, groups, names)
    os.makedirs(args.out_dir, exist_ok=True)
    for condition, df in out.items():
        df.to_csv(os.path.join(args.out_dir, str(condition) + '.csv'), index=False)
    print(str(len(out)) + ' conditions -> ' + args.out_dir + '/<condition>.csv'
          + '  (this is the layout screen_stats reads)')
    return 0


def cmd_counts(args):
    samples, odd = read_samples(args.samples)
    lib = pd.read_csv(args.library)
    id_col = 'pegRNA_id' if 'pegRNA_id' in lib.columns else 'gRNA_id'
    for col in (id_col, 'gene_name_m', 'classification'):
        if col not in lib.columns:
            raise SystemExit('library has no ' + col + ' column; count_matrix needs it')

    guide_gene = dict(zip(lib[id_col], lib['gene_name_m']))
    class_gene = dict(zip(lib[id_col], lib['classification']))
    mapping = dict(zip(samples['raw_folder'], samples['sample_id']))

    out = count_matrix(mapping, list(samples['raw_folder']), list(samples['sample_id']),
                       args.counts_dir, args.count_col, guide_gene, class_gene,
                       merge_lanes=not args.no_merge_lanes)

    if args.editor:
        if 'Editor' not in lib.columns:
            raise SystemExit('--editor given but the library has no Editor column')
        keep = lib[lib['Editor'] == args.editor][id_col].tolist()
        before = len(out)
        out = out[out['sgRNA'].isin(keep)].reset_index(drop=True)
        print('Editor filter ' + args.editor + ': ' + str(before) + ' -> ' + str(len(out))
              + ' guides')
    else:
        sys.stderr.write(
            'WARNING: no --editor filter. An ABE count matrix must contain only ABE guides '
            'and a CBE matrix only CBE guides. Leaving the other editor in does not error; '
            'it changes every column\'s depth normalization and quietly destroys the LFC '
            'downstream.\n')

    n_ctrl = int((out['gene'] == 'controlsgRNA').sum())
    if n_ctrl == 0:
        sys.stderr.write(
            'WARNING: zero guides collapsed to gene "controlsgRNA". Only the exact strings '
            '"non-targeting control" and "safe-targeting control" in the classification '
            'column are collapsed. screen_stats builds its empirical null from these, so '
            'with none it has no null.\n')
    else:
        print(str(n_ctrl) + ' control guides under gene controlsgRNA')

    sample_cols = [c for c in out.columns if c not in ('sgRNA', 'gene')]
    zero = [c for c in sample_cols if out[c].sum() == 0]
    if zero:
        sys.stderr.write('WARNING: all-zero count columns: ' + ', '.join(zero) + '\n')

    out.to_csv(args.out, sep='\t', index=False)
    print(str(len(out)) + ' guides x ' + str(len(sample_cols)) + ' samples -> ' + args.out)
    print('Count column: ' + args.count_col
          + '. Next: screen_stats for LFC / empirical p / FDR.')
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('samples', help='print a samples-file skeleton from a config file')
    s.add_argument('config')
    s.set_defaults(func=cmd_samples)

    c = sub.add_parser('compile', help='STEP7a: aggregated CSVs -> per-sample editing tables')
    c.add_argument('--samples', required=True, help='TSV: raw_folder<TAB>sample_id')
    c.add_argument('--crispresso-dir', required=True,
                   help='local folder holding <raw_folder>_crispresso_aggregated.csv')
    c.add_argument('--out-dir', required=True)
    c.add_argument('--no-merge-lanes', action='store_true')
    c.set_defaults(func=cmd_compile)

    m = sub.add_parser('mle', help='STEP7b: merge biological replicates into conditions')
    m.add_argument('--compact-dir', required=True)
    m.add_argument('--groups', required=True,
                   help='one line per condition: name,sample1,sample2,...')
    m.add_argument('--out-dir', required=True)
    m.set_defaults(func=cmd_mle)

    n = sub.add_parser('counts', help='STEP8: per-guide count matrix')
    n.add_argument('--samples', required=True)
    n.add_argument('--counts-dir', required=True, help='local copy of the cluster counts/')
    n.add_argument('--library', required=True)
    n.add_argument('--out', required=True, help='output .txt (tab-separated)')
    n.add_argument('--count-col', default='bc_count',
                   choices=['bc_count', 'matched_guide_count', 'total_guide_count'])
    n.add_argument('--editor', default=None,
                   help="subpool filter on the library's Editor column, e.g. ABE or CBE")
    n.add_argument('--no-merge-lanes', action='store_true')
    n.set_defaults(func=cmd_counts)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
