"""mitom_metadata.py - the GenBank header and /organism of the original records restored in MitoM's output (`mitom
transfer_metadata`, `mitom replace_organism`, and the final write of `run`): accession matching (versions, RefSeq
prefixes, fuzzy candidates), header transfer at the text level, organism index and replacement."""
from __future__ import annotations

import difflib
import functools
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from mito_common import safe_name


def gb_record_accession(rec: "SeqRecord") -> str:
    ver  = rec.annotations.get("sequence_version", "")
    accs = rec.annotations.get("accessions", [])
    acc  = accs[0] if accs else ""
    if not acc:
        acc = rec.name if rec.name and rec.name not in (".", "<unknown>") else ""
    if not acc:
        acc = rec.id  if rec.id  and rec.id  not in (".", "<unknown>") else ""
    if ver and acc and not acc.endswith(f".{ver}"):
        versioned = f"{acc}.{ver}"
    else:
        versioned = acc
    return safe_name(versioned or rec.id or "unknown")


_REFSEQ_PREFIX_PAT = re.compile(r"^([A-Z]{2}_)")
# only accession-shaped IDs lose a trailing version: local sample names such as B64_2 are distinct samples
_NON_REFSEQ_VERSION_PAT = re.compile(r"^[A-Z]{1,4}\d{5,}[._]\d+$", re.I)
_LOCUS_BPCOUNT_RE  = re.compile(r"(LOCUS\s+\S+\s+)\d+(\s+[a-zA-Z]+)")      # LOCUS bp count


def _clean_acc(acc: Optional[str]) -> Optional[str]:
    if not acc:
        return None
    return acc.strip().rstrip(";").strip().replace(" ", "") or None


def _normalize_acc(acc: Optional[str]) -> Optional[str]:
    s = _clean_acc(acc)
    if not s:
        return None
    m = _REFSEQ_PREFIX_PAT.match(s)
    if m:
        prefix = m.group(1)
        body   = re.sub(r"[._]\d+$", "", s[len(prefix):])
        s      = prefix + body
    else:
        if _NON_REFSEQ_VERSION_PAT.match(s):
            s = re.sub(r"[._]\d+$", "", s)
    return s.rstrip(";").strip() or None


@functools.lru_cache(maxsize=4096)
def _acc_variants(acc: Optional[str]) -> Tuple[str, ...]:
    """All meaningful accession variants for lookup matching (cached; a tuple, so the cached value cannot be
    changed by a caller)."""
    clean = _clean_acc(acc)
    if not clean:
        return ()
    norm = _normalize_acc(clean)
    seen: List[str] = []

    def _add(x: Optional[str]) -> None:
        if x and x not in seen:
            seen.append(x)

    _add(clean)
    _add(norm)
    m = re.search(r"[._](\d+)$", clean)
    if m and norm and norm != clean:             # a version was removed (B64_2 must not gain B64 as an alias)
        ver = m.group(1)
        _add(f"{norm}.{ver}")
        _add(f"{norm}_{ver}")
    return tuple(seen)


def _split_gb_records(text: str) -> List[str]:
    """Split a multi-record GenBank text into individual record strings (no empty pseudo-record from newlines after
    the last //)."""
    records: List[str] = []
    buf: List[str] = []
    for line in text.splitlines(keepends=True):
        buf.append(line)
        if line.strip() == "//":
            chunk = "".join(buf)
            if chunk.strip():          # skip empty / whitespace-only chunks
                records.append(chunk)
            buf = []
    if buf:
        chunk = "".join(buf)
        if chunk.strip():
            records.append(chunk)
    return records


def _split_meta_body(record: str) -> Tuple[str, str]:
    lines = record.splitlines(keepends=True)
    meta: List[str] = []
    body: List[str] = []
    in_body = False
    for line in lines:
        if not in_body and (line.startswith("FEATURES") or line.startswith("ORIGIN")):
            in_body = True
        (body if in_body else meta).append(line)
    return "".join(meta), "".join(body)


def _get_record_accession(record: str) -> Optional[str]:
    """Accession from GB text: VERSION first, then ACCESSION, then LOCUS (Biopython truncates the LOCUS name to 16
    characters)."""
    for prefix in ("VERSION", "ACCESSION", "LOCUS"):
        for line in record.splitlines():
            if line.startswith(prefix):
                parts = line.split()
                if len(parts) >= 2:
                    return _clean_acc(parts[1])
    return None


def _meta_is_meaningful(meta: str) -> bool:
    if not meta.strip():
        return False
    lines = [ln for ln in meta.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("LOCUS"):
        return False
    useful = ("DEFINITION", "ACCESSION", "VERSION", "KEYWORDS",
              "SOURCE", "REFERENCE", "COMMENT", "DBLINK", "ORGANISM")
    return any(any(ln.startswith(p) for p in useful) for ln in lines[1:])


def build_metadata_index(src_text: str) -> Tuple[Dict, Dict, Dict]:
    exact_dict: Dict[str, str] = {}
    clean_dict: Dict[str, str] = {}
    norm_dict:  Dict[str, Tuple[str, str]] = {}
    for rec in _split_gb_records(src_text):
        raw_acc = _get_record_accession(rec)
        if not raw_acc:
            continue
        meta, _ = _split_meta_body(rec)
        if not _meta_is_meaningful(meta):
            continue
        for v in _acc_variants(raw_acc):
            exact_dict[v] = meta
        cl = _clean_acc(raw_acc)
        if cl:
            clean_dict.setdefault(cl, meta)
        norm = _normalize_acc(raw_acc)
        if norm:
            norm_dict.setdefault(norm, (raw_acc, meta))
    return exact_dict, clean_dict, norm_dict


def find_metadata_from_index(acc: str, exact_dict, clean_dict, norm_dict,
                              ) -> Tuple[Optional[str], Optional[str]]:
    for v in _acc_variants(acc):
        if v in exact_dict:
            mtype = "exact" if v == acc else f"clean '{acc}'->'{v}'"
            return exact_dict[v], mtype
    cl = _clean_acc(acc)
    if cl and cl in clean_dict:
        return clean_dict[cl], f"clean '{acc}'->'{cl}'"
    norm = _normalize_acc(acc)
    if norm and norm in norm_dict:
        orig, meta = norm_dict[norm]
        return meta, f"compat '{orig}'<->'{acc}'"
    return None, None


def _similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a.upper(), b.upper()).ratio()


def diagnose_unmatched(acc: str, exact_dict: Dict) -> Tuple[List[str], List[Tuple[int, str]]]:
    hints: List[str] = []
    candidates: List[Tuple[int, str]] = []
    acc_cl   = _clean_acc(acc)     or acc
    acc_norm = _normalize_acc(acc) or acc
    acc_up   = acc_cl.upper()
    seen: set = set()
    for k in exact_dict:
        ck = _clean_acc(k)
        if not ck or ck in seen:
            continue
        seen.add(ck)
        orig_norm = _normalize_acc(ck) or ck
        score = 0
        reason: Optional[str] = None
        if ck.upper() == acc_up and ck != acc_cl:
            score = 90
            reason = f"[case mismatch] orig='{k}' new='{acc}'"
        elif orig_norm == acc_norm and ck != acc_cl:
            score = 88
            reason = f"[version/style mismatch] orig='{k}' new='{acc}'"
        else:
            sim = _similarity(acc_norm, orig_norm)
            if sim > 0.78:
                score = int(sim * 55)
        if score > 0:
            if reason:
                hints.append(reason)
            candidates.append((score, k))
    if not hints:
        hints.append(f"[no diagnosis] no similar accession found for '{acc}'")
    candidates.sort(key=lambda x: -x[0])
    return hints[:5], candidates[:5]


def interactive_resolve(acc: str, candidates: List[Tuple[int, str]],
                        exact_dict: Dict, auto_accept: bool = False,
                        ) -> Tuple[Optional[str], Optional[str]]:
    if not candidates:
        return None, None
    best_score, best_acc = candidates[0]
    if auto_accept:
        if best_score >= 60:
            logging.info(f"[AUTO-ACCEPT] '{acc}' -> '{best_acc}' (score={best_score})")
            return exact_dict.get(best_acc), best_acc
        logging.warning(f"[AUTO-SKIP] '{acc}' top score {best_score} < 60")
        return None, None
    if not sys.stdin.isatty():
        logging.warning(
            f"[NON-INTERACTIVE] stdin is not a TTY; skipping interactive resolution "
            f"for '{acc}'. Top candidate: '{best_acc}' (score={best_score}). "
            f"Use --auto_accept_metadata yes to resolve automatically.")
        return None, None
    print(f"\n⚠ '{acc}' not matched. Candidates:")
    for i, (score, orig_acc) in enumerate(candidates, start=1):
        print(f"  [{i}] {orig_acc:<30} score={score}")
    print("  [s] skip   [q] quit")
    while True:
        try:
            choice = input("choose: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print()
            return None, None
        if choice == "q":
            sys.exit(0)
        if choice == "s":
            return None, None
        if choice.isdigit():
            idx = int(choice) - 1
            if 0 <= idx < len(candidates):
                chosen = candidates[idx][1]
                return exact_dict.get(chosen), chosen
        print("  invalid, retry.")


def replace_metadata_in_text(
    new_gb_text: str, exact_dict, clean_dict, norm_dict,
    smart: bool = False, auto_accept: bool = False,
) -> Tuple[str, Dict]:
    """Replace GB header sections in *new_gb_text* using the metadata index."""
    stats: Dict = {
        "total_new": 0, "no_accession": 0,
        "exact": 0, "compat": 0, "smart_resolved": 0, "not_found": 0,
        "exact_list": [], "compat_list": [], "smart_list": [], "not_found_list": [],
    }
    out_records: List[str] = []
    for rec in _split_gb_records(new_gb_text):
        stats["total_new"] += 1
        acc = _get_record_accession(rec)
        if not acc:
            stats["no_accession"] += 1
            out_records.append(rec)
            continue
        meta, match_type = find_metadata_from_index(acc, exact_dict, clean_dict, norm_dict)
        if meta is not None:
            _, body = _split_meta_body(rec)
            # the bp count of the source metadata's LOCUS line: the length of the re-annotated sequence
            new_locus_m = re.search(r"LOCUS\s+\S+\s+(\d+)\s+[a-zA-Z]+", rec)
            if new_locus_m:
                real_len = new_locus_m.group(1)
                meta = _LOCUS_BPCOUNT_RE.sub(rf"\g<1>{real_len}\g<2>", meta, count=1)
            out_records.append(meta + body if body else rec)
            if match_type == "exact":
                stats["exact"] += 1
                stats["exact_list"].append(acc)
            else:
                stats["compat"] += 1
                stats["compat_list"].append((acc, match_type))
            continue
        hints, candidates = diagnose_unmatched(acc, exact_dict)
        if smart and candidates:
            resolved_meta, chosen_acc = interactive_resolve(
                acc, candidates, exact_dict, auto_accept=auto_accept)
            if resolved_meta is not None:
                _, body = _split_meta_body(rec)
                out_records.append(resolved_meta + body if body else rec)
                stats["smart_resolved"] += 1
                stats["smart_list"].append((acc, chosen_acc))
                continue
        stats["not_found"] += 1
        stats["not_found_list"].append((acc, hints))
        out_records.append(rec)
    return "".join(out_records), stats


def transfer_gb_metadata_file(
    metadata_src_gb: Path, new_body_gb: Path, output_gb: Path,
    smart: bool = False, auto_accept: bool = False,
) -> Tuple[bool, int, int, List[str], str, List[str]]:
    if not metadata_src_gb.exists():
        return False, 0, 0, [], f"metadata source not found: {metadata_src_gb}", []
    if not new_body_gb.exists():
        return False, 0, 0, [], f"body source not found: {new_body_gb}", []
    src_text = metadata_src_gb.read_text(encoding="utf-8", errors="replace")
    new_text = new_body_gb.read_text(encoding="utf-8", errors="replace")
    exact_dict, clean_dict, norm_dict = build_metadata_index(src_text)
    if not exact_dict and not norm_dict:
        return False, 0, 0, [], "metadata source has no usable metadata", []
    replaced_text, stats = replace_metadata_in_text(
        new_text, exact_dict, clean_dict, norm_dict, smart=smart, auto_accept=auto_accept)
    output_gb.parent.mkdir(parents=True, exist_ok=True)
    output_gb.write_text(replaced_text, encoding="utf-8")
    n_replaced  = stats["exact"] + stats["compat"] + stats["smart_resolved"]
    unmatched   = [x[0] for x in stats["not_found_list"]]
    compat_keys = [x[0] for x in stats["compat_list"]]
    return True, n_replaced, stats["not_found"], unmatched, "ok", compat_keys


def build_organism_index(gb_files: List[Path]) -> Dict[str, Tuple[str, str]]:
    """Accession -> (organism_name, source_line_text) lookup, under every accession variant (_acc_variants, which
    include the normalised form)."""
    index: Dict[str, Tuple[str, str]] = {}
    for gb_path in gb_files:
        try:
            records = list(SeqIO.parse(str(gb_path), "genbank"))
        except Exception as exc:
            logging.warning(f"[ORG_INDEX] Cannot parse {gb_path}: {exc}")
            continue
        for rec in records:
            acc = gb_record_accession(rec)

            organism: str = rec.annotations.get("organism", "").strip()
            if not organism:
                for feat in rec.features:
                    if feat.type == "source":
                        vals = feat.qualifiers.get("organism", [])
                        if vals:
                            organism = str(vals[0]).strip()
                        break
            if not organism:
                organism = rec.annotations.get("source", "").strip()

            source_text: str = rec.annotations.get("source", "").strip()

            if not organism and not source_text:
                continue

            entry: Tuple[str, str] = (organism, source_text)

            for variant in _acc_variants(acc):
                if variant in index and index[variant] != entry:      # first wins, the conflict is logged
                    logging.warning(
                        f"[ORG_INDEX_CONFLICT] {variant}: keeping first "
                        f"{index[variant]!r}, ignoring {entry!r} from {gb_path.name}")
                else:
                    index.setdefault(variant, entry)

    return index


def lookup_organism(
    acc: str, org_index: Dict[str, Tuple[str, str]]
) -> Optional[Tuple[str, str]]:
    """Return (organism, source_text) for *acc* using exact → normalised fallback."""
    for v in _acc_variants(acc):
        if v in org_index:
            return org_index[v]
    norm = _normalize_acc(acc)
    if norm and norm in org_index:
        return org_index[norm]
    return None


def replace_organism_in_record(
    rec: "SeqRecord",
    new_organism: str,
    new_source: Optional[str] = None,
) -> "SeqRecord":
    """A copy of *rec* with the organism name replaced throughout: rec.annotations["organism"] when new_organism is
    non-empty, rec.annotations["source"] only when new_source is given and non-empty, /organism= on every source
    feature (added when absent)."""
    new_annotations = dict(rec.annotations)
    if new_organism:
        new_annotations["organism"] = new_organism
    if new_source:
        new_annotations["source"] = new_source

    new_features: List["SeqFeature"] = []
    for feat in rec.features:
        new_quals = dict(feat.qualifiers)
        if feat.type == "source" and new_organism:
            new_quals["organism"] = [new_organism]
        new_features.append(
            SeqFeature(
                location=feat.location,
                type=feat.type,
                id=getattr(feat, "id", "<unknown>"),
                qualifiers=new_quals,
            )
        )

    new_rec = SeqRecord(
        rec.seq,
        id=rec.id, name=rec.name, description=rec.description,
        dbxrefs=list(getattr(rec, "dbxrefs", [])),
        annotations=new_annotations,
        features=new_features,
    )
    return new_rec


def replace_organism_in_gb_file(
    gb_in:     Path,
    gb_out:    Path,
    org_index: Dict[str, Tuple[str, str]],
) -> Tuple[int, int, int]:
    """Replace /organism= (and optionally SOURCE) in every record of *gb_in*."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    out_recs: List["SeqRecord"] = []
    n_replaced = n_not_found = 0
    for rec in recs:
        acc    = gb_record_accession(rec)
        result = lookup_organism(acc, org_index)
        if result is not None:
            organism, source_text = result
            out_recs.append(
                replace_organism_in_record(rec, organism,
                                           new_source=source_text if source_text else None)
            )
            n_replaced += 1
            logging.info(
                f"[REPLACE_ORG] {acc}: /organism= → '{organism}'"
                + (f"  SOURCE → '{source_text}'" if source_text else "")
            )
        else:
            out_recs.append(rec)
            n_not_found += 1
            logging.warning(
                f"[REPLACE_ORG_MISS] {acc}: no entry in organism index — unchanged.")
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return n_replaced, n_not_found, len(recs)


def _sync_organism_qualifier(gb_in: Path, gb_out: Path) -> int:
    """Sync /organism= in FEATURES source to match rec.annotations["organism"] (added when absent)."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    out_recs: List["SeqRecord"] = []
    n_fixed = 0
    for rec in recs:
        organism = rec.annotations.get("organism", "").strip()
        if not organism:
            out_recs.append(rec)
            continue
        new_features: List["SeqFeature"] = []
        changed = False
        for feat in rec.features:
            new_quals = dict(feat.qualifiers)
            if feat.type == "source":
                current = str(new_quals["organism"][0]).strip() \
                          if new_quals.get("organism") else ""
                if current != organism:
                    new_quals["organism"] = [organism]
                    changed = True
            new_features.append(
                SeqFeature(
                    location=feat.location,
                    type=feat.type,
                    id=getattr(feat, "id", "<unknown>"),
                    qualifiers=new_quals,
                )
            )
        new_rec = SeqRecord(
            rec.seq,
            id=rec.id, name=rec.name, description=rec.description,
            dbxrefs=list(getattr(rec, "dbxrefs", [])),
            annotations=dict(rec.annotations),
            features=new_features,
        )
        out_recs.append(new_rec)
        if changed:
            n_fixed += 1
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return n_fixed


def gb_file_has_meaningful_metadata(gb_path: Path) -> bool:
    if not gb_path.is_file():
        return False
    try:
        txt = gb_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return any(_meta_is_meaningful(_split_meta_body(r)[0]) for r in _split_gb_records(txt))
