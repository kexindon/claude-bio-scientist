"""pegRNA Golden Gate oligo design — the library behind `pegrna_oligo`.

Lifted verbatim from the validated notebook `pegrna_gga_oligo_order.ipynb` so the
two cannot drift: overhangs read off the SnapGene maps, BsmBI scanning, dedup,
IDT scale by length, and in-silico assembly.

Cloning cases (backbone + epeg decide the 3' extension oligos):

    backbone            epeg         up OH   down OH   3' ext top      3' ext bottom
    UPEmS-tevo/UPEC-tevo 'backbone'   CACC    CGCG      gtgc+ext        cgcg+rc(ext)
    UPEmS/UPEC           'oligo'      CACC    TTTT      gtgc+ext+tevo   aaaa+rc(tevo)+rc(ext)
    UPEmS/UPEC           'none'       CACC    TTTT      gtgc+ext        aaaa+rc(ext)

Impossible backbone/epeg combinations raise rather than yielding dead oligos.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import Iterable, Sequence

import pandas as pd

# --- sequence primitives -----------------------------------------------------

# str.translate is the fastest pure-Python reverse complement.
_COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def rc(seq: str) -> str:
    """Reverse complement."""
    return seq.translate(_COMP)[::-1]


# Fixed parts of the pegRNA cassette. SCAFFOLD is the sgRNA scaffold minus the
# 5' GTTTTAGAG and the 3' GTGC, which live on the flanking oligos as overhangs.
SCAFFOLD = "CTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG"
TEVOPREQ1 = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"

# Golden Gate junctions, as they appear on the plasmid top strand.
OH_VECTOR_UP = "CACC"      # U6 +1 -> protospacer
OH_PROTO_SCAF = "AGAG"     # protospacer -> scaffold (top strand)
OH_SCAF_EXT = "GTGC"       # scaffold -> 3' extension (top strand)

# BsmBI / Esp3I, CGTCTC(1/5). Both strands must be scanned: the site is not
# palindromic.
_BSMBI_FWD = re.compile("CGTCTC")
_BSMBI_REV = re.compile("GAGACG")
_POLYT = re.compile("T{5,}")

# --- backbones ---------------------------------------------------------------

@dataclass(frozen=True)
class Backbone:
    """A BsmBI-linearizable pegRNA vector.

    down_overhang is the 4-nt 5' top-strand overhang the vector presents
    downstream of the insert; it is what decides the 3' extension bottom oligo.
    """

    name: str
    up_overhang: str
    down_overhang: str
    has_tevopreQ1: bool

    def epeg_options(self) -> tuple[str, ...]:
        return ("backbone",) if self.has_tevopreQ1 else ("oligo", "none")


# Verified against the SnapGene maps (Addgene 196036/196037 and the -tevo derivatives).
BACKBONES: dict[str, Backbone] = {
    "UPEC": Backbone("UPEC", "CACC", "TTTT", False),
    "UPEmS": Backbone("UPEmS", "CACC", "TTTT", False),
    "UPEC-tevo": Backbone("UPEC-tevo", "CACC", "CGCG", True),
    "UPEmS-tevo": Backbone("UPEmS-tevo", "CACC", "CGCG", True),
}


def resolve_backbone(backbone: str | Backbone) -> Backbone:
    if isinstance(backbone, Backbone):
        return backbone
    key = str(backbone).strip()
    for candidate in (key, key.replace("_", "-"), key.replace("Cre", "C")):
        if candidate in BACKBONES:
            return BACKBONES[candidate]
    raise KeyError(
        f"Unknown backbone {backbone!r}. Known: {sorted(BACKBONES)}. "
        "For a new vector, build a Backbone(...) or use backbone_from_map()."
    )


def resolve_epeg(bb: Backbone, epeg: str | None) -> str:
    """Validate the epeg choice against the backbone; never guess silently."""
    allowed = bb.epeg_options()
    if epeg is None:
        if bb.has_tevopreQ1:
            return "backbone"
        raise ValueError(
            f"{bb.name} carries no tevopreQ1, so epeg= must be stated explicitly: "
            "'oligo' (3' extension oligos supply the motif) or 'none' (plain pegRNA)."
        )
    epeg = epeg.lower()
    if epeg not in allowed:
        raise ValueError(
            f"epeg={epeg!r} is impossible for {bb.name} "
            f"(tevopreQ1 in backbone: {bb.has_tevopreQ1}). Allowed: {allowed}."
        )
    return epeg


# --- column presets ----------------------------------------------------------

@dataclass(frozen=True)
class ColumnMap:
    """Where to find protospacer / 3' extension / design id in the input table.

    gene and mut are only used to label the oligos; None leaves that part of the
    name out.
    """

    protospacer: str
    extension: str
    design_id: str
    gene: str | None = None
    mut: str | None = None
    strip_prefix: int = 0  # PRIDICT ships oligo columns already carrying overhangs


COLUMN_PRESETS: dict[str, ColumnMap] = {
    "pegg": ColumnMap("Protospacer", "RTT_PBS", "pegRNA_id",
                      gene="Hugo_Symbol", mut="HGVSp_m"),
    "pridict": ColumnMap(
        "Protospacer-Oligo-FW", "Extension-Oligo-FW", "rank", strip_prefix=4
    ),
}

# --- oligo construction ------------------------------------------------------

@dataclass(frozen=True)
class Duplex:
    """An annealed oligo pair, named by its role in the assembly."""

    role: str
    top: str
    bot: str

    @property
    def left(self) -> str:
        """4-nt 5' overhang carried on the top strand."""
        return self.top[:4]

    @property
    def right(self) -> str:
        """4-nt 5' overhang carried on the bottom strand."""
        return self.bot[:4]

    @property
    def top_contribution(self) -> str:
        """What this duplex adds to the final top strand (overhang + duplexed core)."""
        return self.top

    def validate(self) -> None:
        if len(self.top) != len(self.bot):
            raise ValueError(
                f"{self.role}: top ({len(self.top)} nt) and bottom ({len(self.bot)} nt) "
                "must be equal length for 4-nt overhangs on both ends."
            )
        core = self.top[4:]
        if not rc(self.bot).startswith(core):
            raise ValueError(f"{self.role}: oligos do not anneal — core is not complementary.")


def normalize_protospacer(proto: str, g_start: str = "prepend") -> tuple[str, str]:
    """Force a G at position 1 for the U6 promoter. Returns (spacer, note)."""
    p = proto.strip().upper()
    if p.startswith("G"):
        return p, ""
    if g_start == "prepend":
        return "G" + p, f"G prepended ({len(p)}->{len(p) + 1} nt spacer)"
    if g_start == "replace":
        return "G" + p[1:], f"5' {p[0]}->G substituted"
    raise ValueError(f"Protospacer {p} does not start with G and g_start='require'.")


def protospacer_duplex(spacer: str) -> Duplex:
    # GTTTT/AAAAC is the 5' end of the sgRNA scaffold, split across the junction.
    return Duplex("protospacer", "cacc" + spacer + "gtttt", "ctctaaaac" + rc(spacer))


def scaffold_duplex() -> Duplex:
    return Duplex("scaffold", "agag" + SCAFFOLD, "gcac" + rc(SCAFFOLD))


def extension_duplex(ext: str, bb: Backbone, epeg: str) -> Duplex:
    """3' extension (RTT+PBS, 5'->3'), with the epeg motif where it belongs."""
    ext = ext.strip().upper()
    if epeg == "oligo":
        # Oligos carry tevopreQ1; the duplex ends on the vector's TTTT junction.
        return Duplex("extension", "gtgc" + ext + TEVOPREQ1,
                      "aaaa" + rc(TEVOPREQ1) + rc(ext))
    # epeg 'backbone' -> CGCG junction; epeg 'none' -> TTTT junction.
    return Duplex("extension", "gtgc" + ext, rc(bb.down_overhang).lower() + rc(ext))

# --- QC ----------------------------------------------------------------------

def bsmbi_sites(seq: str) -> list[tuple[int, str]]:
    s = seq.upper()
    hits = [(m.start(), "CGTCTC") for m in _BSMBI_FWD.finditer(s)]
    hits += [(m.start(), "GAGACG") for m in _BSMBI_REV.finditer(s)]
    return sorted(hits)


def qc_design(spacer: str, ext: str, insert_top: str) -> list[str]:
    """Design-level problems that would silently ruin the cloning or the pegRNA."""
    issues: list[str] = []
    for site, where in ((spacer, "protospacer"), (ext, "3' extension")):
        for pos, motif in bsmbi_sites(site):
            issues.append(f"BsmBI site {motif} at {where} position {pos + 1} — assembly will fail")
    # Sites can also be created across a junction, so scan the whole insert.
    extra = len(bsmbi_sites(insert_top))
    if extra:
        issues.append(f"{extra} BsmBI site(s) in the assembled insert (check junctions)")
    for m in _POLYT.finditer(ext.upper()):
        issues.append(
            f"poly-T run ({m.end() - m.start()} nt) at 3' extension position {m.start() + 1} "
            "— may terminate Pol III transcription early"
        )
    return issues


# --- IDT formatting ----------------------------------------------------------

def idt_scale(length: int) -> str:
    """Smallest IDT tube scale that covers this length (15-60/10-90/5-100 nt)."""
    if length <= 60:
        return "25nm"
    if length <= 90:
        return "100nm"
    if length <= 100:
        return "250nm"
    return "ULTRAMER"


# --- oligo naming ------------------------------------------------------------

IDT_NAME_LIMIT = 30


def name_token(text) -> str:
    """Turn a gene symbol or HGVS protein change into an oligo-name token.

    IDT only accepts letters, digits, '_' and '-', so 'p.R451*' becomes 'R451X'
    and 'D594_F595insNS' becomes 'D594F595insNS'. Nonsense '*' maps to 'X' rather
    than being dropped, which would silently make R451* read as R451.
    """
    s = str(text).strip()
    if not s or s.lower() in {"nan", "none"}:
        return ""
    if s[:2].lower() in {"p.", "c.", "g."}:
        s = s[2:]
    return re.sub(r"[^A-Za-z0-9]+", "", s.replace("*", "X"))


def name_stem(prefix: str, index: int) -> str:
    """`prefix_index`, or just `index` when prefix is empty — never a bare '_'."""
    prefix = str(prefix or "").strip().strip("_")
    return f"{prefix}_{index}" if prefix else str(index)


def build_oligo_name(prefix: str, index: int, gene: str, mut: str,
                     code: str, tag: str, limit: int = IDT_NAME_LIMIT) -> str:
    """prefix_index_gene_mut_CODE_TAG, trimmed to IDT's name limit.

    `prefix` is optional: NAME="" gives index_gene_mut_CODE_TAG with no leading
    underscore. The index alone already makes the name unique, so when gene+mut
    push past the limit the labels are shortened (mutation first, then gene)
    instead of raising.
    """
    head, suffix = name_stem(prefix, index), f"_{code}_{tag}"
    gene, mut = name_token(gene), name_token(mut)
    budget = limit - len(head) - len(suffix)
    if budget < 0:
        raise ValueError(
            f"{head}{suffix} is already {len(head) + len(suffix)} characters; "
            f"IDT allows {limit}. Shorten `name`."
        )
    for _ in range(2):
        label = "_".join(t for t in (gene, mut) if t)
        over = (len(label) + 1 if label else 0) - budget
        if over <= 0:
            break
        if len(mut) > over:          # trim the mutation, keep the gene whole
            mut = mut[:-over]
        elif mut:
            mut = ""
        else:
            gene = gene[:max(0, len(gene) - over)]
    label = "_".join(t for t in (gene, mut) if t)
    return f"{head}_{label}{suffix}" if label else f"{head}{suffix}"


def well_positions(n: int, order: str = "column") -> list[str]:
    """96-well labels, column-major by default (matches IDT plate upload)."""
    rows, cols = "ABCDEFGH", range(1, 13)
    if order == "column":
        wells = [f"{r}{c}" for c in cols for r in rows]
    else:
        wells = [f"{r}{c}" for r in rows for c in cols]
    if n > len(wells):
        raise ValueError(f"{n} oligos exceed one 96-well plate; split the order.")
    return wells[:n]

# --- main entry point --------------------------------------------------------

def build_oligos(
    designs: pd.DataFrame,
    name: str,
    backbone: str | Backbone,
    epeg: str | None = None,
    *,
    columns: str | ColumnMap = "pegg",
    gene_col: str | None = "__preset__",
    mut_col: str | None = "__preset__",
    include_scaffold: bool = False,
    g_start: str = "prepend",
    dedup: bool = True,
    plate: bool = False,
    strict: bool = True,
    start_index: int = 1,
    seen: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build IDT-ready oligos for a set of pegRNA designs.

    Returns (orders, report): `orders` is the IDT upload table, `report` is one
    row per design with the oligo names, assembled insert and QC findings.

    `start_index` offsets the oligo index, and passing a `seen` dict shares the
    dedup table between calls — together they let several backbones be built as
    one continuously numbered order.
    """
    bb = resolve_backbone(backbone)
    epeg = resolve_epeg(bb, epeg)
    cols = COLUMN_PRESETS[columns] if isinstance(columns, str) else columns
    if gene_col != "__preset__":
        cols = replace(cols, gene=gene_col)
    if mut_col != "__preset__":
        cols = replace(cols, mut=mut_col)

    for col in (cols.protospacer, cols.extension):
        if col not in designs.columns:
            raise KeyError(f"Column {col!r} not in input ({list(designs.columns)[:8]}...)")
    for col, what in ((cols.gene, "gene_col"), (cols.mut, "mut_col")):
        if col is not None and col not in designs.columns:
            raise KeyError(
                f"{what}={col!r} is not a column of the input "
                f"({list(designs.columns)[:8]}...). Set {what}=None to leave it "
                "out of the oligo names."
            )

    scaffold = scaffold_duplex()
    scaffold.validate()

    order_rows: list[dict] = []
    report_rows: list[dict] = []
    if seen is None:
        seen = {}  # sequence -> first oligo name, for dedup

    proto_series = designs[cols.protospacer].to_numpy()
    ext_series = designs[cols.extension].to_numpy()
    ids = (
        designs[cols.design_id].to_numpy()
        if cols.design_id in designs.columns
        else range(1, len(designs) + 1)
    )
    blanks = [""] * len(designs)
    genes = designs[cols.gene].to_numpy() if cols.gene else blanks
    muts = designs[cols.mut].to_numpy() if cols.mut else blanks

    for n, (raw_proto, raw_ext, design_id, gene, mut) in enumerate(
        zip(proto_series, ext_series, ids, genes, muts), start_index
    ):
        proto_in = str(raw_proto)[cols.strip_prefix:]
        ext_in = str(raw_ext)[cols.strip_prefix:]
        spacer, g_note = normalize_protospacer(proto_in, g_start)
        ext = ext_in.strip().upper()

        p = protospacer_duplex(spacer)
        e = extension_duplex(ext, bb, epeg)
        duplexes = [p, scaffold, e] if include_scaffold else [p, e]
        for d in (p, scaffold, e):
            d.validate()

        insert_top = (p.top + scaffold.top + e.top).upper()
        check_assembly_overhangs(bb, [p, scaffold, e], strict=strict)

        issues = qc_design(spacer, ext, insert_top)
        if issues and strict:
            raise ValueError(
                f"{name_stem(name, n)} ({design_id}) failed QC:\n  - "
                + "\n  - ".join(issues)
            )

        stem = build_oligo_name(name, n, gene, mut, "P", "T")[:-4]

        # entries: (oligo name, sequence, name of the identical oligo already ordered)
        entries: list[tuple[str, str, str | None]] = []
        for d in duplexes:
            code = {"protospacer": "P", "scaffold": "S", "extension": "E"}[d.role]
            for tag, seq in (("T", d.top), ("B", d.bot)):
                oligo_name = build_oligo_name(name, n, gene, mut, code, tag)
                key = seq.upper()
                duplicate_of = seen.get(key) if dedup else None
                if duplicate_of is None:
                    seen[key] = oligo_name
                    order_rows.append(
                        {
                            "Name": oligo_name,
                            "Sequence": seq,
                            "Scale": idt_scale(len(seq)),
                            "Purification": "STD",
                            "Length": len(seq),
                        }
                    )
                entries.append((oligo_name, seq, duplicate_of))

        report_rows.append(
            {
                "oligo_index": n,
                "design_id": design_id,
                "name_stem": stem,
                "gene": name_token(gene),
                "mut": name_token(mut),
                "spacer": spacer,
                "spacer_note": g_note,
                "extension": ext,
                "backbone": bb.name,
                "epeg": epeg,
                "oligos_ordered": ",".join(nm for nm, _, dup in entries if dup is None),
                "reused": ",".join(f"{nm}->{dup}" for nm, _, dup in entries if dup),
                "insert_bp": len(insert_top),
                "qc": "; ".join(issues) if issues else "ok",
            }
        )

    orders = pd.DataFrame(order_rows)
    if plate and not orders.empty:
        orders.insert(0, "Well Position", well_positions(len(orders)))
    return orders, pd.DataFrame(report_rows)


def check_assembly_overhangs(bb: Backbone, duplexes: Sequence[Duplex], strict: bool = True) -> None:
    """Every junction must join a complementary pair of 4-nt overhangs."""
    problems = []
    if duplexes[0].left.upper() != bb.up_overhang:
        problems.append(f"vector {bb.up_overhang} vs protospacer top {duplexes[0].left.upper()}")
    for a, b in zip(duplexes, duplexes[1:]):
        if rc(a.right.upper()) != b.left.upper():
            problems.append(f"{a.role}->{b.role}: {a.right.upper()} / {b.left.upper()}")
    last = duplexes[-1].right.upper()
    if rc(last) != bb.down_overhang:
        problems.append(f"{duplexes[-1].role}->vector: {last} vs {bb.down_overhang}")
    if problems and strict:
        raise ValueError("Overhang mismatch:\n  - " + "\n  - ".join(problems))

# --- in-silico Golden Gate against a real map --------------------------------

def read_snapgene(path: str) -> str:
    from Bio import SeqIO

    return str(SeqIO.read(path, "snapgene").seq).upper()


def digest_bsmbi(seq: str) -> tuple[str, str, str]:
    """Digest a circular plasmid with exactly two BsmBI sites.

    Returns (backbone_top, up_overhang, down_overhang) where backbone_top starts
    with the downstream 5' top overhang and ends flush before the upstream one.
    """
    fwd = [m.start() for m in _BSMBI_FWD.finditer(seq)]
    rev = [m.start() for m in _BSMBI_REV.finditer(seq)]
    if len(fwd) != 1 or len(rev) != 1:
        raise ValueError(f"Expected one site per strand, found {len(fwd)} CGTCTC / {len(rev)} GAGACG")
    i, j = fwd[0], rev[0]
    down = seq[i + 7:i + 11]           # 5' top overhang on the backbone
    up = seq[j - 5:j - 1]              # top-strand spelling of the bottom overhang
    backbone_top = seq[i + 7:] + seq[:j - 5] if j < i else seq[i + 7:j - 5]
    return backbone_top, up, down


def backbone_from_map(path: str, name: str | None = None) -> Backbone:
    """Derive a Backbone from a SnapGene map rather than trusting the registry."""
    seq = read_snapgene(path)
    _, up, down = digest_bsmbi(seq)
    return Backbone(name or path.rsplit("/", 1)[-1], up, down, TEVOPREQ1 in seq)


def assemble(vector_map: str, spacer: str, ext: str, epeg: str) -> dict:
    """Run the Golden Gate in silico and verify the product."""
    seq = read_snapgene(vector_map)
    backbone_top, up, down = digest_bsmbi(seq)
    bb = Backbone(vector_map.rsplit("/", 1)[-1], up, down, TEVOPREQ1 in seq)

    p, s, e = protospacer_duplex(spacer), scaffold_duplex(), extension_duplex(ext, bb, epeg)
    for d in (p, s, e):
        d.validate()
    check_assembly_overhangs(bb, [p, s, e])

    insert_top = (p.top + s.top + e.top).upper()
    product = backbone_top + insert_top
    expected_cassette = (
        OH_VECTOR_UP + spacer + "GTTTTAGAG" + SCAFFOLD + OH_SCAF_EXT + ext
        + (TEVOPREQ1 if epeg == "oligo" else "")
    )
    return {
        "backbone": bb,
        "product_bp": len(product),
        "product": product,
        "insert_bp": len(insert_top),
        "cassette_ok": expected_cassette.upper() in product,
        "residual_bsmbi": bsmbi_sites(product + product[:5]),
        "tevo_in_product": product.count(TEVOPREQ1),
    }

