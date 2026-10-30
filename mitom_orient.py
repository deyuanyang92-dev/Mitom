"""mitom_orient.py - re-orientation of annotated records: the orientation gene (cox1, else --fallback_genes) and its
aliases, rotation / reverse complement of a record with its features (partial ends kept), the check whether the
shifted round-1 annotation can be kept or round 2 is needed (--if_reannotation_2 auto, reannotation_decision.tsv),
and the rotation to the first exon after the intron step."""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from mito_common import (abbreviated_stop, feature_start_pos, merge_abutting_parts, rc_location as _rc_location,  # noqa: F401
                         move_position_qualifiers, reports_dir, reverse_complement_record, shift_location, write_gb,
                         write_nonempty_table, write_table)
from mito_names import ORIENTATION_ALIAS_GROUPS, gene_aliases, name_key

if TYPE_CHECKING:
    from mitom_inputs import SampleInfo


def gb_to_internal_fasta(gb_in: Path, fasta_out: Path, internal_id: str) -> None:
    """Convert a GenBank file to FASTA for Phase-2 re-annotation input: all records, the second contig as
    <internal_id>_2 and so on."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    out_recs = []
    for i, rec in enumerate(recs):
        topo   = rec.annotations.get("topology", "circular")
        seq_id = internal_id if i == 0 else f"{internal_id}_{i + 1}"
        out_recs.append(SeqRecord(rec.seq, id=seq_id, name=seq_id,
                                  description=f"topology={topo}"))
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with fasta_out.open("w", encoding="utf-8") as fh:
        SeqIO.write(out_recs, fh, "fasta")


normalize_gene_key = name_key
_GENE_ALIAS_GROUPS = ORIENTATION_ALIAS_GROUPS


# Records oriented by a fallback gene because the first gene was absent: {record id: gene}
ORIENTATION_FALLBACK: Dict[str, str] = {}


def orientation_genes(args: argparse.Namespace) -> str:
    """'cox1' plus --fallback_genes, as the comma list reorient_one_gb tries in order."""
    fb = [g.strip() for g in (getattr(args, "fallback_genes", "") or "").split(",")
          if g.strip() and g.strip().lower() != "none"]
    return ",".join([args.which_gene_first] + [g for g in fb if g != args.which_gene_first])


def find_orientation_feature(record: "SeqRecord", genes: str) -> Tuple[Optional["SeqFeature"], str]:
    """First gene of the comma list present in the record; when the gene is annotated in
    several pieces (e.g. exons around an intron), the 5'-most piece in its own orientation."""
    for gene in [g.strip() for g in genes.split(",") if g.strip()]:
        aliases = set(gene_aliases(gene))
        hits = [f for f in record.features if f.type in {"gene", "CDS", "tRNA", "rRNA"}
                and any(normalize_gene_key(n) in aliases for n in feature_gene_names(f))]
        if not hits:
            continue
        strand = hits[0].location.strand
        same = [f for f in hits if f.location.strand == strand]
        if strand == -1:
            return max(same, key=lambda f: int(f.location.end)), gene
        return min(same, key=lambda f: int(f.location.start)), gene
    return None, ""


def feature_gene_names(feature: "SeqFeature") -> List[str]:
    names: List[str] = []
    for k in ("gene", "product", "label"):
        for v in feature.qualifiers.get(k, []):
            names.append(str(v))
    return names


_LINEAR_SKIP_SEEN: set = set()      # sequences already reported by [LINEAR_REORIENT_SKIP] (run: step 3, then the final write)


def rotate_record(record: "SeqRecord", start_pos: int) -> "SeqRecord":
    seq_len = len(record.seq)
    topology = record.annotations.get("topology", "circular")
    if str(topology).lower() == "linear" and start_pos != 0:
        key = hash(str(record.seq))
        (logging.info if key in _LINEAR_SKIP_SEEN else logging.warning)(
            f"[LINEAR_REORIENT_SKIP] record={record.id}: target start={start_pos}; "
            "linear topology cannot be circularly rotated, sequence left unchanged.")
        _LINEAR_SKIP_SEEN.add(key)
        return SeqRecord(                     # moving a linear contig's tail before its head would corrupt it
            record.seq,
            id=record.id, name=record.name, description=record.description,
            dbxrefs=list(getattr(record, "dbxrefs", [])),
            annotations=dict(record.annotations),
            features=list(record.features),
        )
    new_seq = record.seq[start_pos:] + record.seq[:start_pos]
    new_rec = SeqRecord(
        Seq(str(new_seq)),
        id=record.id, name=record.name, description=record.description,
        dbxrefs=list(getattr(record, "dbxrefs", [])),
        annotations={**record.annotations, "topology": topology},
    )
    new_features: List["SeqFeature"] = []
    for feat in record.features:
        try:
            new_loc = shift_location(feat.location, start_pos, seq_len)
            new_features.append(SeqFeature(location=new_loc, type=feat.type,
                                            id=getattr(feat, "id", "<unknown>"),
                                            qualifiers=move_position_qualifiers(feat.qualifiers, seq_len, start_pos)))
        except Exception:
            new_features.append(feat)
    new_features.sort(key=lambda f: feature_start_pos(f))
    new_rec.features = new_features
    return new_rec


# --if_reannotation_2 auto (user 2026-10-21: re-annotating every record after reorientation wastes resources)
# ATP8 and ATP6 are NOT in the list: absent in whole phyla (nematodes, flatworms, many bivalves, anthozoans...),
# so their absence says nothing about a cut gene (review B finding 4).
STANDARD_GENES = ("COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                  "s-rRNA", "l-rRNA")
REANNOT_EDGE = 50          # nt from the ends of the round-1 sequence where a gene may have been cut
TRNA_LEN = (55, 95)        # complete tRNA genes are this long (nt); outside = cut or merged
RRNA_MIN = {"s": 600, "l": 900}   # shorter s-rRNA / l-rRNA (nt) at the sequence ends = cut


def edge_gene_truncated(f: "SeqFeature", seq: "Seq", total: int, code: int = 5) -> bool:
    """A CDS / rRNA / tRNA feature within REANNOT_EDGE nt of an end of the round-1 sequence (length `total`):
    does it look cut there? (A complete gene that merely starts at 1 or ends at the last base is not: NCBI
    vertebrate records start with trnF.) CDS: length not a multiple of 3 (except a T / TA stop tail, the
    mitochondrial incomplete stop codon), or no start codon where its 5' end is the sequence end, or no stop
    codon in the feature nor in frame within REANNOT_EDGE nt behind it where its 3' end is the sequence end
    (codon tables of `code`); tRNA: outside 55-95 nt; rRNA: s-rRNA < 600 nt, l-rRNA < 900 nt (unknown rRNA
    name: < 600). A join across the origin is an origin-spanning feature, judged like any interior one."""
    loc = f.location
    wraps = len(loc.parts) > 1 and any(int(p.end) == total for p in loc.parts) \
        and any(int(p.start) == 0 for p in loc.parts)       # join across the origin: annotated as such, its ends are inside
    left, right = (False, False) if wraps else (int(loc.start) < REANNOT_EDGE, int(loc.end) > total - REANNOT_EDGE)
    if not (left or right or wraps):
        return False
    n = len(loc)
    if f.type == "tRNA":
        return not TRNA_LEN[0] <= n <= TRNA_LEN[1]
    if f.type == "rRNA":
        name = " ".join(f.qualifiers.get("gene", []) + f.qualifiers.get("product", [])).lower()
        large = any(k in name for k in ("l-rrna", "rrnl", "16s", "large", "lsu"))
        return n < RRNA_MIN["l" if large else "s"]
    if left and right:                         # spans the whole sequence: cannot judge, re-annotate
        return True
    from Bio.Data import CodonTable
    table = CodonTable.unambiguous_dna_by_id.get(int(code), CodonTable.unambiguous_dna_by_id[5])
    nt = str(f.extract(seq)).upper()
    rem = len(nt) % 3
    tail = abbreviated_stop(nt)
    if rem and not tail:
        return True
    minus = loc.strand == -1
    five_at_edge = right if minus else left       # which end of the gene is a sequence end
    three_at_edge = left if minus else right
    if five_at_edge and nt[:3] not in table.start_codons:
        return True
    if three_at_edge and not tail and nt[-3:] not in table.stop_codons:
        # the core annotates a CDS up to its last aligned codon (the stop is added later): a stop in frame within
        # REANNOT_EDGE nt behind it, or a T / TA tail ending exactly at the sequence end, makes it complete
        if minus:
            beyond = str(seq[max(0, int(loc.start) - REANNOT_EDGE):int(loc.start)].reverse_complement()).upper()
            reaches_end = int(loc.start) <= REANNOT_EDGE
        else:
            beyond = str(seq[int(loc.end):int(loc.end) + REANNOT_EDGE]).upper()
            reaches_end = int(loc.end) + REANNOT_EDGE >= total
        codons = [beyond[i:i + 3] for i in range(0, len(beyond) - 2, 3)]
        if any(c in table.stop_codons for c in codons):
            return False
        return not (reaches_end and beyond[len(codons) * 3:] in ("T", "TA"))
    return False


def select_reannotation(samples: List["SampleInfo"], results1: Dict[str, Dict], out_root: Path,
                        code: int = 5) -> List["SampleInfo"]:
    """--if_reannotation_2 auto: the samples whose reoriented round-1 annotation fails reorient_check()
    (genetic code `code`); every decision with its reasons goes to reannotation_decision.tsv."""
    chosen, rows = [], []
    for si in samples:
        iid = si.internal_id
        reor = out_root / "01.anno_round1" / iid / "reoriented" / f"{iid}.reoriented.gbf"
        best = results1.get(iid, {}).get("best_gb")
        try:
            probs = reorient_check(Path(best), reor, code) if best else ["no round-1 GenBank"]
        except Exception as exc:                     # any doubt -> re-annotate
            probs = [f"check failed: {exc}"]
        if probs:
            chosen.append(si)
        rows.append(f"{si.sample_name}\t{'reannotate' if probs else 'shifted_round1'}\t{'; '.join(probs)}")
    write_table(reports_dir(out_root) / "reannotation_decision.tsv", "sample_name\tdecision\treasons", rows)
    logging.info(f"[PHASE2] auto: {len(chosen)} of {len(samples)} sample(s) re-annotated, "
                 f"{len(samples) - len(chosen)} use the shifted round-1 annotation (reannotation_decision.tsv)")
    return chosen


def reorient_check(round1_gb: Path, reor_gb: Path, code: int = 5) -> List[str]:
    """Reasons to re-annotate a reoriented record; [] = the shifted round-1 annotation can be used as is.
    Checks: (1) a gene within REANNOT_EDGE nt of the ends of the round-1 sequence that looks cut there
    (edge_gene_truncated, genetic code `code`; the core annotates the sequence as linear, so a gene across the
    original origin is cut or missed) - only when the sequence moved; (2) a standard protein gene / rRNA
    (STANDARD_GENES; not ATP8 / ATP6) missing in round 1; (3) every feature's sequence identical
    before and after the shift, same features; (4) no feature across the new origin. A linear record is never rotated
    (only turned to the plus strand): its end genes are cut by the real sequence ends and a missing gene is missing
    from the sequence, so (1) and (2) do not apply - a second round on the same sequence finds nothing new (2026-10
    benchmark, 15 linear annelid records: the same true tRNAs, 2 more false ones from the re-run)."""
    probs: List[str] = []
    r1 = {r.id: r for r in SeqIO.parse(str(round1_gb), "genbank")}
    for ro in SeqIO.parse(str(reor_gb), "genbank"):
        a = r1.get(ro.id)
        if a is None:
            probs.append(f"{ro.id}: no round-1 record")
            continue
        L = len(a.seq)
        genes = [f for f in a.features if f.type in ("CDS", "rRNA", "tRNA")]
        moved = str(ro.seq) != str(a.seq)
        linear = "linear" in (a.annotations.get("topology"), ro.annotations.get("topology"))
        if moved and not linear:
            edge = sorted({f.qualifiers.get("gene", [f.type])[0] for f in genes
                           if edge_gene_truncated(f, a.seq, L, code)})
            if edge:
                probs.append(f"{ro.id}: gene(s) at the end of the round-1 sequence, look cut: {', '.join(edge)}")
        have = {f.qualifiers.get("gene", [""])[0] for f in genes}
        miss = [g for g in STANDARD_GENES if g not in have]
        if miss and not linear:
            probs.append(f"{ro.id}: standard gene(s) missing in round 1: {', '.join(miss)}")

        def seqs(rec: "SeqRecord") -> Dict[Tuple[str, str], List[str]]:
            out: Dict[Tuple[str, str], List[str]] = {}
            for f in rec.features:
                if f.type in ("CDS", "rRNA", "tRNA"):
                    out.setdefault((f.type, f.qualifiers.get("gene", [""])[0]), []).append(str(f.extract(rec.seq)))
            return {k: sorted(v) for k, v in out.items()}
        before, after = seqs(a), seqs(ro)
        diff = sorted(k[1] or k[0] for k in set(before) | set(after) if before.get(k) != after.get(k))
        if diff:
            probs.append(f"{ro.id}: gene sequence differs after the shift: {', '.join(diff)}")
        n = len(ro.seq)
        cross = sorted({f.qualifiers.get("gene", [f.type])[0] for f in ro.features if f.type != "source"
                        and len(f.location.parts) > 1 and any(int(p.end) == n for p in f.location.parts)
                        and any(int(p.start) == 0 for p in f.location.parts)})
        if cross:
            probs.append(f"{ro.id}: feature(s) across the new origin: {', '.join(cross)}")
    return probs


def reorient_one_gb(
    gb_in: Path, gb_out: Path, gb_no_gene: Path,
    which_gene_first: str, force_forward: bool,
) -> Tuple[int, int, List[str]]:
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No GenBank records: {gb_in}")
    ok: List["SeqRecord"] = []
    no_gene: List["SeqRecord"] = []
    no_gene_ids: List[str] = []
    for rec in recs:
        try:
            feat, used = find_orientation_feature(rec, which_gene_first)
            if feat is None:
                no_gene.append(rec)
                no_gene_ids.append(rec.id)
                continue
            first = which_gene_first.split(",")[0].strip()
            if used != first:
                ORIENTATION_FALLBACK[rec.id] = used
                logging.info(f"[REORIENT_FALLBACK] record={rec.id}: '{first}' not found, "
                             f"sequence starts at '{used}'")
            work_rec = rec
            if force_forward and feat.location.strand == -1:
                work_rec = reverse_complement_record(rec)
                feat, _ = find_orientation_feature(work_rec, used)
                if feat is None:
                    no_gene.append(rec)
                    no_gene_ids.append(rec.id)
                    continue
            ok.append(rotate_record(work_rec, feature_start_pos(feat)))
        except Exception as exc:
            logging.warning(f"[REORIENT_REC_FAIL] record={rec.id}: {exc}")
            no_gene.append(rec)
            no_gene_ids.append(rec.id)
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    gb_no_gene.parent.mkdir(parents=True, exist_ok=True)
    if ok:
        SeqIO.write(ok, str(gb_out), "genbank")
    if no_gene:
        SeqIO.write(no_gene, str(gb_no_gene), "genbank")
    return len(ok), len(no_gene), no_gene_ids


def starts_with_gene(gb: Path, which_gene_first: str, force_forward: bool) -> bool:
    """Does a record of `gb` begin with its orientation gene (the first of which_gene_first present; on the plus
    strand when force_forward), i.e. no gene feature lies before it?"""
    for rec in SeqIO.parse(str(gb), "genbank"):
        feat, _ = find_orientation_feature(rec, which_gene_first)
        if feat is None or (force_forward and feat.location.strand == -1):
            continue
        genes = [f for f in rec.features if f.type in ("CDS", "rRNA", "tRNA")]
        if genes and feature_start_pos(feat) <= min(feature_start_pos(f) for f in genes):
            return True
    return False


def orientation_text(orient_ok: bool, si: "SampleInfo") -> str:
    """cox1_orientation_ok column: yes / no / no (linear, not rotated)."""
    return "yes" if orient_ok else "no (linear, not rotated)" if si.topology == "linear" else "no"


def write_orientation_fallback(out_root: Path, samples: List["SampleInfo"]) -> None:
    """oriented_by_fallback.txt: samples whose first gene was absent and that start at a fallback gene."""
    names = {si.internal_id: si.sample_name for si in samples}
    rows = sorted({f"{names.get(k.split('.')[0], k)}\t{g}" for k, g in ORIENTATION_FALLBACK.items()})
    write_nonempty_table(reports_dir(out_root) / "oriented_by_fallback.txt", "sample_name\tstart_gene", rows)


def add_exon_gap_introns(rec: "SeqRecord") -> int:
    """Intron features for the gaps between consecutive exons of a spliced CDS that lie inside the record and have no
    intron feature (introns across the origin are left out by mito_intron and restored here once the record is
    rotated): /gene and /number from the CDS, /note from its "intron across the origin" note when it has one.
    Returns the number added."""
    from Bio.SeqFeature import FeatureLocation
    from mito_intron import INTRON_NOTE
    have = [f for f in rec.features if f.type == "intron"]
    added = []
    for cds in [f for f in rec.features if f.type == "CDS" and len(f.location.parts) > 1]:
        parts = list(cds.location.parts)                      # transcript order
        gene = cds.qualifiers.get("gene", [""])[0]
        notes, ref_notes = {}, {}
        for x in cds.qualifiers.get("note", []):
            if x.startswith(INTRON_NOTE):              # mito_intron (de novo): type, E-value, splice score
                k, _, txt = x[len(INTRON_NOTE):].partition(":")
                notes[k.strip()] = txt.strip()
            m = re.match(r"intron (\d+) \([^)]*\) crosses the origin: given by the join only; (.*)", x)
            if m:                                      # mito_intron_ref: where the site is known from
                ref_notes[m.group(1)] = m.group(2)
        for k, (a, b) in enumerate(zip(parts, parts[1:]), 1):
            if a.strand != b.strand:
                continue
            s, e = (int(a.end), int(b.start)) if a.strand != -1 else (int(b.end), int(a.start))
            if not 0 < e - s < len(rec.seq) // 2 or s >= e:     # across the origin, abutting or not a gap
                continue
            if any(int(h.location.start) < e and int(h.location.end) > s for h in have + added):
                continue
            q = {"gene": [gene], "number": [str(k)]}
            q["note"] = ["; ".join(x for x in (notes.get(str(k)), ref_notes.get(str(k))) if x)
                         or "intron between the exons of the CDS"]
            added.append(SeqFeature(FeatureLocation(s, e, strand=a.strand), type="intron", qualifiers=q))
        if added:                                      # the introns are features now: the crossing notes are stale
            cds.qualifiers["note"] = [x for x in cds.qualifiers.get("note", []) if not x.startswith(INTRON_NOTE)
                                      and "crosses the origin: given by the join only" not in x]
    if added:
        rec.features = sorted(rec.features + added, key=lambda f: (int(f.location.start), f.type != "source",
                                                                   f.type != "gene"))
        logging.info(f"[REORIENT_INTRON] {rec.id}: {len(added)} intron feature(s) restored inside the record")
    return len(added)


def reorient_after_introns(args: argparse.Namespace, final_dir: Path) -> int:
    """The intron step can turn the orientation gene into exons across the origin: MitoM rotated the record to the
    core's piece of cox1, and when that piece is exon 2, exon 1 and the intron lie at the record's end (annelid cox1
    group II introns). Rotate such circular records again so the gene's 5'-most exon starts the sequence. Returns the
    number of records rotated."""
    genes, n = orientation_genes(args), 0
    for gbf in sorted(final_dir.rglob("*.gbf")):
        recs, changed = list(SeqIO.parse(str(gbf), "genbank")), False
        for i, rec in enumerate(recs):
            feat, _ = find_orientation_feature(rec, genes)
            if rec.annotations.get("topology") != "circular" or feat is None or feat.location.strand == -1:
                continue
            pos = feature_start_pos(feat)
            if pos:
                recs[i], changed = rotate_record(rec, pos), True
                n += 1
                add_exon_gap_introns(recs[i])
                logging.info(f"[REORIENT_INTRON] {rec.id}: {feat.qualifiers.get('gene', ['?'])[0]} starts at {pos + 1} "
                             "after the intron step: rotated to its first exon")
        if changed:
            write_gb(recs, gbf)
    return n
