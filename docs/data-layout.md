# Reference data layout

Every skill in this repo reaches reference data through `$LAB_DATA_ROOT`.
Nothing hardcodes a path, so the same skill runs on a personal laptop and on
the shared lab server, under any username.

```bash
# lab server, shared and read-only to everyone
export LAB_DATA_ROOT=/Users/Shared/lab-data/Database

# personal machine
export LAB_DATA_ROOT=~/Documents/GitHub/Database
```

Put the line in `~/.zshrc` so it survives new shells, then `source ~/.zshrc`.

## Expected tree

```
$LAB_DATA_ROOT/
├── RefGenome/
│   ├── ncbi-2023-09-12/        GRCh37/38 human, bgzipped .fna.gz
│   └── mouse-2023-09-13/       GRCm39 mouse
├── Genecode/
│   ├── release_19_GRCh37.p13/  gencode_v19_GRCh37.db
│   ├── release_44_GRCh38.p14/  gencode_v44_GRCh38.db
│   └── release_m33_GRCm39/     gencode_vm33_GRCm39.db
├── PublicDatabase/
│   ├── AACR-GENIE/             v13.1, v15.0, v17.0  ← controlled access
│   ├── COSMIC/                 Cancer Gene Census, mutant census
│   ├── ClinVar/
│   ├── MSK-IMPACT/
│   ├── MSK-HEME/               heme_msk_impact_2022
│   ├── genomAD/
│   └── MMRRC/
├── Pathogenecity/
│   ├── AlphaMIssense/          AlphaMissense_aa_substitutions.tsv
│   ├── SIFT/
│   └── dbNSFP/
├── Homology/                   HMD_HumanPhenotype.txt (MGI)
├── Encode/
├── H2M/
└── provean-1.1.5/
```

Roughly 72 GB in total. The GENCODE `.db` files are ~1.5 GB each; building
one from a GTF takes about 20 minutes, so they are worth keeping rather than
regenerating.

## Provisioning a new machine

Most of this tree is public data that is cleaner to re-download than to copy,
because you get known provenance instead of inheriting whatever state an old
copy is in:

| Source | Where from |
| --- | --- |
| Reference genomes | NCBI FTP (GRCh37, GRCh38, GRCm39), bgzipped |
| GENCODE annotations | gencodegenes.org — download the GTF, build the `.db` with `gffutils` |
| ClinVar | NCBI ClinVar FTP, VCF |
| COSMIC | cancer.sanger.ac.uk — needs a free academic account |
| AlphaMissense | Zenodo / the paper's supplement |
| MGI homology | informatics.jax.org |

What cannot simply be re-downloaded, and must be copied or re-requested:

- **AACR-GENIE** — controlled access under a data use agreement, obtained
  through Synapse with an approved account
- **MSK-IMPACT / MSK-HEME** — access terms depend on the release
- Anything derived or processed that is not in a public archive

## GENIE and controlled access

AACR-GENIE sits under `LAB_DATA_ROOT`, but unlike the rest of the tree it is
**not public**. It is distributed under a data use agreement that governs who
may access it and whether it may be shared.

Before pointing a shared machine, a second user account, or any automated
agent at `PublicDatabase/AACR-GENIE/`, confirm with whoever administers the
lab's GENIE access that this is covered. Putting it on a server other people
can read is a redistribution question, not a technical one.

Never commit any part of it to a repo.

## Making it read-only on a shared machine

On the lab server the data should be readable by everyone and writable only
by its owner, so another user's script cannot corrupt the reference data that
every analysis depends on:

```bash
sudo mkdir -p /Users/Shared/lab-data
sudo chown "$(whoami):staff" /Users/Shared/lab-data
chmod -R a+rX,go-w /Users/Shared/lab-data
```

## Verifying a copy

After moving the tree, check both ends match before deleting anything:

```bash
# file count and total size
find "$LAB_DATA_ROOT" -type f | wc -l
du -sh "$LAB_DATA_ROOT"

# checksum the large files that matter most
find "$LAB_DATA_ROOT" -name "*.db" -o -name "*.fna.gz" | sort | xargs shasum -a 256
```

Run the same two commands on the source and compare. Only then remove the
original.
