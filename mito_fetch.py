#!/usr/bin/env python3
"""mito_fetch.py - download complete mitochondrial genomes (GenBank flat files) from NCBI nuccore.

Records are selected by taxon (scientific name or taxid; the whole subtree) or by a raw Entrez query.
By default both RefSeq (NC_) and INSDC (GenBank) records are taken; `mitom dedup` later keeps one copy
of each genome (the RefSeq one when it exists). Output: one GenBank file with every record, a report
table, and the query that was used. A second run in the same directory only downloads records that
are not in the file yet (resume / update).

Usage:
  mitom fetch --taxon Annelida -o fetch_out                 # name or taxid, RefSeq + GenBank
  mitom fetch --taxon 6340 --refseq_only -o fetch_out       # RefSeq (NC_) only
  mitom fetch --query 'Nereididae[Organism] AND mitochondrion[Filter] AND "complete genome"[Title]' -o out
  mitom fetch --taxon Annelida -o fetch_out --fasta      # also fetch_out/Annelida.fasta
  then:  mitom dedup -i fetch_out/Annelida.gb -o dedup_out     or     mitom run -i fetch_out/Annelida.gb -o out

NCBI asks for an e-mail with every E-utilities request (--email, $NCBI_EMAIL, or "ncbi_email" in the
mito_env config file) and allows 10 requests/s with an API key (--api_key, $NCBI_API_KEY, "ncbi_api_key")
instead of 3/s without.

Outputs in -o:  <name>.gb  fetch_report.tsv  query.txt  fetch.log  (+ <name>.fasta with --fasta)
"""
from __future__ import annotations

import argparse
import http.client
import logging
import os
import re
import sys
import time
import urllib.error
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
from mito_common import gb_fasta_record, old_spelling, write_tsv  # noqa: E402
from mito_common import cli_guard, InputError, prepare_input, require_records  # noqa: E402

ESEARCH_PAGE = 10000          # NCBI maximum retmax of one esearch request
ACC_BATCH = 500               # UIDs per efetch rettype=acc request
GB_BATCH = 200                # records per efetch rettype=gb request (~15 kb each)
RETRIES = 5
TRANSIENT = (urllib.error.URLError, urllib.error.HTTPError, http.client.HTTPException,
             ConnectionError, TimeoutError, OSError)
REFSEQ_RE = re.compile(r"^[A-Z]{2}_\d")
MT_CORE = 'mitochondrion[Filter] AND "complete genome"[Title] NOT partial[Title]'


# ── query ────────────────────────────────────────────────────────────────────
def build_query(taxid: Optional[int], raw: Optional[str], refseq_only: bool,
                min_len: int = 0, max_len: int = 0) -> str:
    """Entrez term: the taxon subtree (txid[Organism:exp]) with the mitogenome filters, or the raw query as
    given; --refseq_only and the length window are appended to both."""
    if raw:
        term = f"({raw})"
    elif taxid is not None:
        term = f"txid{taxid}[Organism:exp] AND {MT_CORE}"
    else:
        raise ValueError("either a taxon or a query is needed")
    if refseq_only:
        term += " AND srcdb_refseq[PROP]"
    if min_len or max_len:
        term += f" AND {min_len or 1}:{max_len or 999999999}[SLEN]"
    return term


def _retry(call: Callable, what: str):
    """Run an E-utilities call; NCBI returns transient 429/5xx and connection resets under load."""
    delay = 2.0
    for attempt in range(1, RETRIES + 1):
        try:
            return call()
        except TRANSIENT as exc:
            if attempt == RETRIES:
                raise
            logging.warning(f"[FETCH] {what}: {exc} (retry {attempt}/{RETRIES - 1} in {delay:.0f}s)")
            time.sleep(delay)
            delay = min(delay * 2, 60)


def resolve_taxid(taxon: str) -> int:
    """Taxid of a scientific name (exact match in the taxonomy database); a numeric taxon is returned as is."""
    taxon = taxon.strip()
    if taxon.isdigit():
        return int(taxon)
    res = _retry(lambda: Entrez.read(Entrez.esearch(db="taxonomy", term=f'"{taxon}"[Scientific Name]')),
                 f"taxonomy lookup {taxon}")
    ids = res.get("IdList", [])
    if len(ids) != 1:
        if not ids:
            raise SystemExit(f"taxon '{taxon}' not found in the NCBI taxonomy (try the taxid)")
        raise SystemExit(f"taxon '{taxon}' matches {len(ids)} taxonomy records ({', '.join(ids)}); give the taxid")
    return int(ids[0])


def taxon_name(taxid: int) -> str:
    """Scientific name of a taxid (for the output file name); the taxid itself when the lookup fails."""
    try:
        res = _retry(lambda: Entrez.read(Entrez.esummary(db="taxonomy", id=str(taxid))), f"taxonomy {taxid}")
        return str(res[0]["ScientificName"])
    except Exception as exc:  # noqa: BLE001  (name is cosmetic)
        logging.warning(f"[FETCH] cannot read the name of taxid {taxid}: {exc}")
        return f"txid{taxid}"


def search_uids(term: str) -> List[str]:
    """Every nuccore UID matching term (paged esearch)."""
    uids: List[str] = []
    start = 0
    while True:
        res = _retry(lambda: Entrez.read(Entrez.esearch(db="nuccore", term=term, retstart=start,
                                                        retmax=ESEARCH_PAGE, idtype="acc")),
                     f"esearch retstart={start}")
        page = list(res.get("IdList", []))
        uids += page
        count = int(res.get("Count", 0))
        start += len(page)
        if not page or start >= count:
            break
    return uids


def existing_accessions(gb_path: Path) -> Dict[str, str]:
    """accession.version -> accession of the records already in gb_path (resume)."""
    out: Dict[str, str] = {}
    if not gb_path.is_file() or gb_path.stat().st_size == 0:
        return out
    text = gb_path.read_text(encoding="utf-8", errors="replace")
    end = text.rfind("\n//")
    tail = text[end + 3:] if end >= 0 else text
    if tail.strip():                                     # a cut-off last record (killed / out of disk): removed, so
        keep = text[:end + 3] + "\n" if end >= 0 else ""   # that the resume downloads it again (2026-10-31 review)
        logging.warning(f"[FETCH] {gb_path.name}: the last record is cut off ({len(tail)} bytes after the last //): "
                        f"removed, downloaded again")
        gb_path.write_text(keep, encoding="utf-8")
        if not keep:
            return out
    for rec in SeqIO.parse(str(gb_path), "genbank"):
        acc_ver = rec.id
        out[acc_ver] = acc_ver.split(".")[0]
    return out


def _chunks(items: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def download(uids: List[str], gb_path: Path, log_every: int = 1) -> int:
    """efetch the records (GenBank flat file) in batches, appended to gb_path. Returns records written."""
    written = 0
    gb_path.parent.mkdir(parents=True, exist_ok=True)
    with gb_path.open("a", encoding="utf-8") as out:
        for i, batch in enumerate(_chunks(uids, GB_BATCH), 1):
            # a batch is written only when every record in the response is complete (LOCUS ... //): a truncated
            # body or an error page from NCBI under load must not end up in the file, where a resume would keep it
            for attempt in range(1, RETRIES + 1):
                text = _retry(lambda: Entrez.efetch(db="nuccore", id=",".join(batch), rettype="gb",
                                                    retmode="text").read(), f"efetch batch {i}")
                if isinstance(text, bytes):
                    text = text.decode("utf-8", "replace")
                n = len(re.findall(r"^//\s*$", text, re.M))
                n_locus = len(re.findall(r"^LOCUS ", text, re.M))
                if n and n == n_locus:
                    break
                logging.warning(f"[FETCH] batch {i}: incomplete response ({n} complete of {n_locus} records; "
                                f"attempt {attempt}/{RETRIES})")
                time.sleep(2.0 * attempt)
            else:
                raise RuntimeError(f"[FETCH] batch {i}: no complete response after {RETRIES} attempts; nothing of this "
                                   f"batch was written - run the same command again to resume")
            if n < len(batch):
                logging.warning(f"[FETCH] batch {i}: asked {len(batch)} records, got {n} "
                                f"(the missing ones are tried again by the next run)")
            if not text.endswith("\n"):
                text += "\n"
            out.write(text)
            out.flush()
            written += n
            if i % log_every == 0:
                logging.info(f"[FETCH] {written}/{len(uids)} records")
    return written


def write_report(gb_path: Path, report: Path) -> int:
    """fetch_report.tsv: one row per record in gb_path."""
    rows = []
    for rec in (SeqIO.parse(str(gb_path), "genbank") if gb_path.is_file() else []):   # no file: nothing matched
        ann = rec.annotations
        lineage = ann.get("taxonomy", [])
        rows.append({
            "accession": rec.id,
            "refseq": "yes" if REFSEQ_RE.match(rec.id) else "no",
            "organism": ann.get("organism", ""),
            "length": len(rec.seq),
            "topology": ann.get("topology", ""),
            "date": ann.get("date", ""),
            "lineage": "; ".join(lineage),
            "definition": rec.description,
        })
    rows.sort(key=lambda r: (r["organism"], r["accession"]))
    write_tsv(report, list(rows[0].keys()) if rows else
              ["accession", "refseq", "organism", "length", "topology", "date", "lineage", "definition"], rows)
    return len(rows)


def setup_entrez(email: Optional[str], api_key: Optional[str]) -> None:
    cfg = mito_env.config()
    email = email or os.environ.get("NCBI_EMAIL") or cfg.get("ncbi_email")
    api_key = api_key or os.environ.get("NCBI_API_KEY") or cfg.get("ncbi_api_key")
    if not email:
        raise SystemExit("NCBI needs an e-mail address: --email you@example.org "
                         "(or NCBI_EMAIL, or \"ncbi_email\" in the mito_env config file)")
    Entrez.email = email
    Entrez.tool = "mitom"
    if api_key:
        Entrez.api_key = api_key


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--taxon", help="Scientific name or taxid; the whole subtree is taken (e.g. Annelida, 6340).")
    g.add_argument("--query", help="Raw Entrez query for nuccore (used instead of the taxon + mitogenome filters).")
    ap.add_argument("-o", "--out", dest="out_dir", metavar="OUT", default="fetch_out", help="Output folder (default: fetch_out).")
    old_spelling(ap, "--out_dir", "out_dir")
    ap.add_argument("--name", help="Base name of the GenBank file (default: taxon name, or 'ncbi' for --query).")
    ap.add_argument("--refseq_only", action="store_true", help="RefSeq (NC_) records only.")
    ap.add_argument("--min_len", type=int, default=0, help="Minimum sequence length (bp).")
    ap.add_argument("--max_len", type=int, default=0, help="Maximum sequence length (bp).")
    ap.add_argument("--max", type=int, default=0, help="Stop after this many new records (test runs).")
    ap.add_argument("--email", help="E-mail for NCBI E-utilities (or $NCBI_EMAIL / config ncbi_email).")
    ap.add_argument("--api_key", help="NCBI API key (or $NCBI_API_KEY / config ncbi_api_key).")
    ap.add_argument("--dry_run", action="store_true", help="Only print the query and the record count.")
    ap.add_argument("--fasta", action="store_true",
                    help="Also write <name>.fasta (every record of <name>.gb; headers >ACCESSION.V definition "
                         "topology=circular|linear, ready for mitom run). GenBank stays the main output: its "
                         "taxonomy sets the clade and genetic code, its header is copied into the results.")
    args = ap.parse_args(argv)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(out / "fetch.log", encoding="utf-8")])
    setup_entrez(args.email, args.api_key)

    taxid = None
    name = args.name
    if args.taxon:
        taxid = resolve_taxid(args.taxon)
        name = name or re.sub(r"[^A-Za-z0-9_.-]+", "_", taxon_name(taxid) if args.taxon.isdigit() else args.taxon)
        logging.info(f"[FETCH] taxon {args.taxon} -> txid{taxid}")
    name = name or "ncbi"
    term = build_query(taxid, args.query, args.refseq_only, args.min_len, args.max_len)
    (out / "query.txt").write_text(term + "\n", encoding="utf-8")
    logging.info(f"[FETCH] query: {term}")

    uids = search_uids(term)
    logging.info(f"[FETCH] {len(uids)} records match")
    if args.dry_run:
        print(f"query: {term}\nrecords: {len(uids)}")
        return 0

    gb_path = out / f"{name}.gb"
    have = existing_accessions(gb_path)
    bases = set(have.values())
    todo = [u for u in uids if u not in have and u.split(".")[0] not in bases]
    newer = [u for u in uids if u not in have and u.split(".")[0] in bases]
    logging.info(f"[FETCH] already in {gb_path.name}: {len(have)}  to download: {len(todo)}")
    if newer:
        logging.info(f"[FETCH] {len(newer)} record(s) have a newer version at NCBI than the copy in {gb_path.name} "
                     f"(kept as they are; delete the file to refresh): {', '.join(newer[:5])}"
                     + (" ..." if len(newer) > 5 else ""))
    if args.max and len(todo) > args.max:
        todo = todo[:args.max]
        logging.info(f"[FETCH] --max {args.max}: downloading the first {len(todo)}")
    n = download(todo, gb_path) if todo else 0
    total = write_report(gb_path, out / "fetch_report.tsv")
    if args.fasta and gb_path.exists():
        fa_path = out / f"{name}.fasta"
        SeqIO.write((gb_fasta_record(r) for r in SeqIO.parse(str(gb_path), "genbank")), str(fa_path), "fasta")
        logging.info(f"[FETCH] FASTA of the {total} records: {fa_path}")
    logging.info(f"[FETCH] done: {n} new, {total} records in {gb_path}, report {out / 'fetch_report.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
