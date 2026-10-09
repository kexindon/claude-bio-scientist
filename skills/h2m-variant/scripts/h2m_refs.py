"""Reference-data resolver and loader for the `h2m_variant` skill.

The four heavy objects bioh2m needs are reached through $LAB_DATA_ROOT, never
through a hardcoded path. This module is a wrapper around the engine, not a
reimplementation of it: `h2m.genome_loader` and `h2m.anno_loader` do the work.

The load pattern itself is lifted from
  h2m-public/methods/run_h2m/h2m_db.py lines 17-21
and from h2m-public/h2m_package_tutorial.ipynb cell 13. Those files hardcode
bare filenames or /Users/kexindong/... paths; the only change here is that the
paths come from $LAB_DATA_ROOT. Do not edit the loader calls here without
editing there -- these must not drift.

Verified on disk 2026-10-09 against ~/Documents/GitHub/Database.
"""

import os

# bioh2m vendors its own module under the name `module`; the public surface is
# `import bioh2m as h2m`. methods/run_h2m/h2m_db.py still says `import h2m`,
# which is the pre-rename name and no longer resolves. Lift the logic, not that
# import line.
import bioh2m as h2m

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

# NOTE the spelling: the directory on disk is "Genecode", not "Gencode".
# docs/data-layout.md documents the .db files as living inside the per-release
# subdirectories, but on the machines checked they sit at the top of
# Genecode/ instead. Both layouts are searched, top level first.
_GENOME = {
    (37, "h"): ["RefGenome/ncbi-2023-09-12/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"],
    (38, "h"): ["RefGenome/ncbi-2023-09-12/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"],
    (None, "m"): ["RefGenome/mouse-2023-09-13/GCF_000001635.27_GRCm39_genomic.fna.gz"],
}

_ANNO = {
    (37, "h"): [
        "Genecode/gencode_v19_GRCh37.db",
        "Genecode/release_19_GRCh37.p13/gencode_v19_GRCh37.db",
    ],
    (38, "h"): [
        "Genecode/gencode_v44_GRCh38.db",
        "Genecode/release_44_GRCh38.p14/gencode_v44_GRCh38.db",
    ],
    (None, "m"): [
        "Genecode/gencode_vm33_GRCm39.db",
        "Genecode/release_m33_GRCm39/gencode_vm33_GRCm39.db",
    ],
}


def data_root():
    """$LAB_DATA_ROOT, or a one-line failure naming the doc that explains it."""
    root = os.environ.get("LAB_DATA_ROOT")
    if not root:
        raise SystemExit(
            "LAB_DATA_ROOT is not set -- see docs/data-layout.md. "
            "e.g. export LAB_DATA_ROOT=~/Documents/GitHub/Database"
        )
    root = os.path.expanduser(root)
    if not os.path.isdir(root):
        raise SystemExit("LAB_DATA_ROOT points at " + root + ", which is not a directory.")
    return root


def _pick(candidates, what):
    """First existing candidate, or a failure naming every path tried."""
    root = data_root()
    tried = []
    for rel in candidates:
        # Plain concatenation / os.path.join only. An f-string with a nested
        # double-quoted subscript is a SyntaxError before Python 3.11.
        full = os.path.join(root, rel)
        if os.path.exists(full):
            return full
        tried.append(full)
    raise SystemExit(
        "Missing reference file for " + what + ". Looked for:\n  "
        + "\n  ".join(tried)
        + "\nSee docs/data-layout.md for how to provision the tree."
    )


def genome_path(species, ver=None):
    """Path to the NCBI reference FASTA. species is 'h' or 'm'."""
    if species == "m":
        return _pick(_GENOME[(None, "m")], "mouse genome (GRCm39)")
    if ver not in (37, 38):
        raise SystemExit("ver must be 37 or 38 for human; got " + repr(ver))
    return _pick(_GENOME[(ver, "h")], "human genome (GRCh" + str(ver) + ")")


def anno_path(species, ver=None):
    """Path to the prebuilt gffutils .db. species is 'h' or 'm'."""
    if species == "m":
        return _pick(_ANNO[(None, "m")], "mouse GENCODE db (vM33)")
    if ver not in (37, 38):
        raise SystemExit("ver must be 37 or 38 for human; got " + repr(ver))
    return _pick(_ANNO[(ver, "h")], "human GENCODE db (GRCh" + str(ver) + ")")


# ---------------------------------------------------------------------------
# Preflight -- turn a mid-run KeyError into a setup error
# ---------------------------------------------------------------------------

def preflight(ver=37):
    """Check the .db files exist and are paired with the packaged tx dicts.

    `get_tx_batch` reads dict_tx_h.json / dict_tx_m.json frozen inside the
    installed bioh2m wheel. If the .db was built from a different GENCODE
    release, those transcript IDs are absent from it and the failure surfaces
    much later as a bare `KeyError` inside `GetTx` (db_input[tx_id]). Catch it
    here instead, in under a second, with no genome loaded.
    """
    import pandas as pd

    problems = []
    for species, probe in (("h", "TP53"), ("m", "Trp53")):
        path = anno_path(species, ver if species == "h" else None)
        db = h2m.anno_loader(path)
        col = "gene_name_h" if species == "h" else "gene_name_m"
        df = pd.DataFrame({col: [probe], "index": [0]})
        if species == "h":
            ok, _bad = h2m.get_tx_batch(df, "h", ver=ver)
        else:
            ok, _bad = h2m.get_tx_batch(df, "m")
        if len(ok) == 0:
            problems.append(probe + " is not in the packaged transcript dict")
            continue
        tx_id = ok.iloc[0]["tx_id_h" if species == "h" else "tx_id_m"]
        try:
            db[tx_id]
        except Exception:
            problems.append(
                "transcript " + str(tx_id) + " (from the bioh2m wheel) is not in "
                + os.path.basename(path)
                + " -- the GENCODE release and the wheel disagree"
            )
    if problems:
        raise SystemExit("Reference preflight failed:\n  - " + "\n  - ".join(problems))
    return True


# ---------------------------------------------------------------------------
# The heavy load, cached per process
# ---------------------------------------------------------------------------

_CACHE = {}


def load_refs(ver=37, quiet=False):
    """Return (records_h, index_list_h, records_m, index_list_m, db_h, db_m).

    ~1 minute and several GB on a local SSD -- `genome_loader` parses the whole
    FASTA into a list of Biopython records (the source comment says ~4 GB for
    human) and both species must be resident at once. Cached on `ver` so a
    notebook session pays it once. Never call this inside a loop.

    Nothing here is a quick question: if all you need is an ortholog or a
    canonical transcript, use h2m.query / h2m.get_tx_batch, which read packaged
    JSON and need no genome at all.
    """
    if ver in _CACHE:
        return _CACHE[ver]

    p_h = genome_path("h", ver)
    p_m = genome_path("m")
    a_h = anno_path("h", ver)
    a_m = anno_path("m")

    if not quiet:
        print("loading human genome  " + p_h)
    records_h, index_list_h = h2m.genome_loader(p_h)
    if not quiet:
        print("loading mouse genome  " + p_m)
    records_m, index_list_m = h2m.genome_loader(p_m)
    if not quiet:
        print("loading annotations   " + os.path.basename(a_h) + ", " + os.path.basename(a_m))
    db_h, db_m = h2m.anno_loader(a_h), h2m.anno_loader(a_m)

    _CACHE[ver] = (records_h, index_list_h, records_m, index_list_m, db_h, db_m)
    return _CACHE[ver]


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Resolve and check h2m reference data.")
    ap.add_argument("--ver", type=int, default=37, choices=(37, 38))
    ap.add_argument("--load", action="store_true", help="actually load the genomes (slow)")
    args = ap.parse_args()

    print("LAB_DATA_ROOT = " + data_root())
    for label, path in (
        ("human genome", genome_path("h", args.ver)),
        ("mouse genome", genome_path("m")),
        ("human anno  ", anno_path("h", args.ver)),
        ("mouse anno  ", anno_path("m")),
    ):
        print(label + "  " + path)
    preflight(args.ver)
    print("preflight OK (transcript dicts match the .db files)")
    if args.load:
        import time

        t0 = time.time()
        load_refs(args.ver)
        print("genomes loaded in " + str(round(time.time() - t0, 1)) + " s")


if __name__ == "__main__":
    main()
