# MitoM code review — Kimi K3, 2026-11-04

Reviewed by: **Kimi K3** (independent review; existing reports in `Mitom_review/` were not read).
Scope: production code of P:\Mitom (~40 top-level modules, Python). Read-only review: **no source code, tests or
configuration were changed**; this document is the only file written. Tests/, benchmark/, knowledge/ were only used
to confirm intended behaviour where cited.

Method: every production module was read; suspect lines were re-read before recording a problem. Items I could not
tie to specific lines were dropped. Biological statements I am not certain of are marked **(uncertain)** and kept at
low severity.

Findings are ordered by severity inside each section. Total: **36 problems** (4 high, 15 medium, 17 low).

---

## 1. Software engineering standards

### 1.1 HIGH — the output-folder lock kills the healthy run on Windows
`mitom_pipeline.py:1470-1489` (`lock_output`)

```python
other = int(lock.read_text().split()[0])
os.kill(other, 0)                      # still running?
```
On POSIX `os.kill(pid, 0)` probes a process. On Windows CPython implements `os.kill` via `TerminateProcess(handle,
sig)` for every signal except `CTRL_C_EVENT`/`CTRL_BREAK_EVENT` — so `os.kill(other, 0)` **terminates the other run's
process with exit code 0** when a second `mitom run` / `annotate` / `import` is pointed at an output folder that a
live run holds the lock of. The lock exists to prevent exactly this corruption (half-written outputs of two runs in
one folder); on Windows the probe kills the healthy run mid-write and then reports it as running.
Trigger: two `mitom run -o same_folder` on Windows.
Fix: probe with `if os.name != "nt": os.kill(other, 0) else:` use `psutil.Process(pid).is_running()` (already a
benchmark dependency) or ctypes `OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, ...)`; treat an unreadable pid as
stale only when the open genuinely fails with "no such process".

### 1.2 HIGH — multiprocessing hard-codes the `fork` context (breaks on Windows; plots fail silently)
Four production sites:
- `mito_intron.py:1612-1614` — `with mp.get_context("fork").Pool(a.jobs) as pool:` (used when `--jobs > 1`)
- `mito_intron_ref.py:445-447` (`build_library`) and `mito_intron_ref.py:2013-2016` (`main`, `--jobs > 1`)
- `mitom_steps.py:387-389` (`run_plots`) — created unconditionally, even for `--jobs 1`

`mp.get_context("fork")` raises `ValueError: cannot find context for 'fork'` on Windows. In the intron tools with
`--jobs > 1` this is an unhandled crash (traceback, exit 1). In `run_plots` the `except Exception` at
`mitom_steps.py:393` converts it into `ok=False` — **every genome map is silently missing on Windows** (the run
"succeeds"). Note also that `run_plots` forks a pool even when `max_tasks == 1`.
Trigger: `mitom intron --jobs 4` or any `mitom run` (plot step) on Windows.
Fix: pick the context by platform (`"fork" if os.name != "nt" else "spawn"`) — spawn requires picklable job state
(`_JOB` is a plain dict of picklable values in both intron tools; `_plot_file`'s job tuple is picklable too) — or
fall back to the sequential loop when fork is unavailable.

### 1.3 HIGH — nhmmer output redirected to a literal `/dev/null`
`mito_rrna_hmm.py:242`

```python
p = run_tool([exe, "--dna", "--cpu", "2", "-E", str(EVALUE), "--tblout", str(tbl), "-o", "/dev/null",
              str(hmm), str(fa)])
```
There is no `/dev/null` on Windows: nhmmer cannot open its `-o` file, exits non-zero, and `nhmmer()` raises
`RuntimeError` for every search. The core's rRNA step (`annotate_rrnas`) and the rRNA HMM fallback then fail for
every record on Windows (surfaced only as `failed:` rows in `rrna_hmm_report.tsv` — silently missing rRNAs).
Trigger: `mitom run` / `mito_rrna_hmm.py search` on Windows.
Fix: `"-o", os.devnull` (or `tempfile`/`"-"` stdout discard as the tool allows).

### 1.4 MEDIUM — `mito_intron.py` crashes when a directory input needs normalization
`mito_intron.py:1597-1598` and `mito_intron.py:1537`

```python
src_in = prepare_input(a.input)
files = iter_files([src_in], GB_EXT)
...
src_root = Path(a.input) if Path(a.input).is_dir() else f.parent
dest = out / f.relative_to(src_root)
```
`prepare_input()` unpacks `.gz` / converts BOM, Latin-1 and EMBL files into a temporary folder and returns that
folder. When the input is a *directory* needing such conversion, `files` live under the temp dir while `src_root` is
still the original `a.input`, so `f.relative_to(src_root)` raises `ValueError` and the whole run dies with a
traceback (cli_guard re-raises non-Biopython ValueErrors). `mito_intron_ref.py` does this correctly
(`root=src.parent if src.is_file() else src`, `mito_intron_ref.py:2010`).
Trigger: `python mito_intron.py -i folder_with_gb.gz_files/ -o out` (or a folder of EMBL/Latin-1 files).
Fix: pass the prepared `src_in` into `_JOB` and compute `src_root` from it, as in `mito_intron_ref.py`.

### 1.5 HIGH — value-less qualifiers (`/pseudo`) written back as `/pseudo=""` by most GenBank writers
Biopython parses a value-less `/pseudo` as `qualifiers["pseudo"] == [""]` and writes `[""]` back as `/pseudo=""`,
which is not valid INSDC (`/pseudo` is a flag; the project itself hit this: `mito_common.py:105-107`
`FLAG_QUALIFIERS` with the comment "2026-10-31 review: /pseudo=\"\" in the output"). `mito_common.write_gb()` fixes
these (and writes atomically), but many output paths still use plain `SeqIO.write(..., "genbank")`:

- `mito_curate.py:1376` — the curated output itself: `SeqIO.write(recs, str(dest), "genbank")`
- `mitom_orient.py:295,297` (`reorient_one_gb` — used by `mitom reorient` and by `finalize_genbank_output` in every
  run/import), `mitom_output.py:148` (`rewrite_genbank_record_ids`, step 1 of the final write)
- `mitom_metadata.py:470,513` (`replace_organism_in_gb_file`, `_sync_organism_qualifier`)
- `mitom_pipeline.py:346` (`prepare_import` — every `mitom import`), `mitom_pipeline.py:1035` (`workflow_dedup`)
- `mito_ncr.py:133`, `mito_trna.py:862`, `mito_engines.py:361`, `mitom_stage.py:141` (`write_input_annotation`)

Trigger: `mitom curate` / `import` / `reorient` / `dedup` on a GenBank record carrying `/pseudo`,
`/ribosomal_slippage` or `/trans_splicing` (e.g. an NCBI record with a pseudogene) → the deliverable GenBank contains
`/pseudo=""`. (In-pipeline, `mito_curate`'s own P10 writes `[None]` correctly, so the corruption needs a
parse→write round trip — exactly what these commands do.)
Fix: route every GenBank write through `mito_common.write_gb` (atomic + flag fix); that also removes the remaining
non-atomic writers (a killed run can currently leave a truncated `.gb`).

### 1.6 MEDIUM — `.tbl` writer gives `/pseudo` a value (`pseudo\tNone` / empty value)
`ncbi_validate.py:55-61` (`feature_lines`)

```python
for k, vals in f.qualifiers.items():
    ...
    for v in vals:
        ...
        out.append(f"\t\t\t{k}\t{' '.join(str(v).split())}")
```
No None/empty handling: in-memory features from `mito_curate` P10 carry `pseudo=[None]` → the 5-column table gets
`\t\t\tpseudo\tNone` (a qualifier value of the literal text "None"). Features re-read from disk carry `[""]` →
`\t\t\tpseudo\t` with a trailing tab and an empty 5th column instead of the bare value-less form. The None case
fires in `mito_curate`'s own before/after validation (`mito_curate.py:1372` validates the in-memory curated record);
the empty-value case fires in `run_submission` and `gb2tbl`.
Trigger: any pseudogene in records validated by table2asn through these paths.
Fix: `if v in (None, ""): out.append(f"\t\t\t{k}"); continue`.

### 1.7 MEDIUM — `transfer_metadata` swallows the exception text and always exits 0
`mitom_pipeline.py:945-955` (`workflow_transfer_metadata`)

```python
except Exception:
    if args.overwrite_new_body == "no":
        safe_copy(bf, out_file)
    n_fail += 1
    log_rows.append(f"{bf}\t{matched or 'NA'}\t{out_file}\tfail\t0\t0")
...
return 0
```
A bare `except Exception` with the exception never logged anywhere (compare `workflow_replace_organism`, which logs
`[FAIL] {bf.name}: {exc}`), and `n_fail` has no effect on the exit code: a run that failed for every input file
reports `[SUMMARY] transferred=0 failed=N` and **exits 0** — silent "success". `workflow_reorient` was fixed for
exactly this (`mitom_pipeline.py:504-508` returns 1).
Fix: `logging.error(f"[FAIL] {bf.name}: {exc!r}")`; `return 1 if n_fail else 0`.

### 1.8 MEDIUM — `mito_nuc`'s tool runner discards stderr of failed programs and bypasses the [CMD] log
`mito_nuc.py:142-143`

```python
def run(cmd: List[str]) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
```
On failure (e.g. makeblastdb/blastn/exonerate error) `CalledProcessError` propagates with `exc.stderr` captured but
never printed or logged — the traceback only says "returned non-zero exit status N", the program's own error text is
lost. It also bypasses `mito_common.run_tool`, so none of `mitom nuc`'s external calls appear as `[CMD]` lines in
`nuc.log`. (`mito_classify._blastx`, `mito_classify.py:125-132`, is also a raw `subprocess.run` — it at least embeds
the stderr tail in its RuntimeError, but its blastx calls are likewise invisible to the run's command log.)
Fix: use `mito_common.run_tool` / `run_or_raise` (logs `[CMD rc=...]` + stderr tail), or print `exc.stderr` on failure.

### 1.9 LOW-MEDIUM — per-sample log handler leaks on the error path (locks the file on Windows)
`mitom_stage.py:155-195` (`annotate_one`)

`get_sample_logger(log_name, workdir / "sample.log")` opens a `FileHandler`; it is closed only on the two happy
paths (`close_sample_logger` at lines 191/195). Exceptions in between — `engine_cmd(...)` construction,
`run_cmd_capture` raising `OSError` when the engine program is missing, `write_input_annotation` errors — propagate
to `run_annotation_stage`'s `except` with the handler still open. The named logger also stays in the logging
registry. On Windows the locked `sample.log` makes `cleanup_intermediates`' `shutil.rmtree(ignore_errors=True)`
silently leave work dirs behind (and `clear_previous_work` then raises `PermissionError` on the next run).
Fix: `try: ... finally: close_sample_logger(log_name)` around the whole body.

### 1.10 LOW — locale-dependent text I/O in output files (non-ASCII corrupted on Windows)
- `mitom_output.py:72,75` (`merge_gb_files`): `gb.read_text()` / `merged_out.write_text(...)` without
  `encoding="utf-8"` — on a cp1252 Windows locale, UTF-8 GenBank files with non-ASCII content (e.g. an organism
  "Müller") are mis-decoded and written back mojibake into `final_annotation_all.gb`.
- `mito_core.py:783`: `(out / f"{prefix}.core.tsv").open("w")` without encoding.
- `trna_ml.py:171`: `json.load(open(cdir / "models.json"))` without encoding (and an unclosed handle).
Most other writers in the project correctly pass `encoding="utf-8"`.
Fix: add `encoding="utf-8"` (and `with` blocks) to these.

### 1.11 LOW — a feature that cannot be transformed is kept with its old coordinates, silently
`mitom_orient.py:112-118` (`rotate_record`) and `mito_common.py:294-301` (`reverse_complement_record`)

```python
except Exception:
    new_features.append(feat)
```
If `shift_location` / `rc_location` raises for a feature (unusual position type, malformed location), the *original*
feature is appended unchanged — it then sits at **wrong coordinates** on the rotated / reverse-complemented
sequence, with nothing logged. Dropping it (or failing loudly) is safer than a silently misplaced gene.
Fix: at minimum `log.warning` with record id and feature; preferably drop or abort the record.

### 1.12 LOW — executable discovery does not probe Windows suffixes
`mito_env.py:153` (`_exe`: only `p / name` with `os.access(..., X_OK)`) and `mitom:97-99` (the launcher probes
`b / "python3"` / `b / "python"` in conda env bins). On Windows the binaries are `blastx.exe`, `python.exe`, etc.;
the config/env-bin scan never finds them (only the `shutil.which` PATH fallback works), and the launcher fails to
locate a conda Python with Biopython on Windows.
Fix: try `name + ext` for each `PATHEXT` suffix (or use `shutil.which(name, path=dir)`).

### 1.13 LOW — output of failed programs is still parsed
`mito_intron.py:199-213` (the module's `run_tool` wrapper): on `p.returncode != 0` it warns once per program but
**returns `p.stdout` anyway** — a truncated exonerate/blastp stream is then parsed as if complete, which can yield
partial HSPs and wrong models. (Timeout/OSError correctly return `""`.)
Fix: return `""` on non-zero exit after warning, or have callers treat non-empty stderr/exit as "no result".

### 1.14 LOW — group-I subtype CMs: non-atomic download, cached forever when partial
`mito_intron.py:551-564` (`subtype_cms`): `urllib.request.urlretrieve(it["download_url"], SUBTYPE_CM_DIR /
it["name"])` writes straight to the final cache path with no timeout, checksum or size check. A killed/interrupted
download leaves a truncated `.cm` that every later run picks up (`glob("*.cm")`), making cmsearch fail or type
introns wrongly.
Fix: download to a temp file, verify (size/`cmpress`/manifest), then `os.replace`.

### 1.15 LOW — input normalization reads each whole file into memory twice
`mitom_inputs.py:129-131` (`p.read_bytes().decode("utf-8")` to test decodability) and again at
`mitom_inputs.py:149`. For very large FASTA inputs this doubles peak memory of `mitom run`'s first step.
Fix: decode streaming in chunks (e.g. `codecs.iterdecode`) or check a bounded prefix plus incremental validation.

### 1.16 LOW — lock cleanup at exit can itself raise
`mitom_pipeline.py:1486-1487`: the `atexit` lambda does `lock.is_file()` then `lock.read_text()` — a TOCTOU race
(the file can vanish in between) raises `FileNotFoundError` inside the atexit handler, printing a traceback at
interpreter exit and potentially masking the run's real exit path.
Fix: wrap the body in `try/except OSError`.

### 1.17 LOW — duplicated logic and nominal module boundaries (maintainability)
- `mito_compare._lift_except` (lines 118-133) re-implements `mito_common.move_positions` (including the
  origin-split join) — two copies of the same `/transl_except` coordinate math to keep in sync.
- "Most common CDS `/transl_table`" is implemented three times: `mito_curate.record_code` (1289-1292),
  `ncbi_validate.record_gencode` (77-79), `mito_check.record_code` (90-93) — with the same tie-handling quirk
  (`max(set(codes), key=codes.count)` is order-unstable on ties).
- Same-named helpers with different semantics: `mito_ncr.span`/`extent` vs `mito_compare.span` vs
  `mito_rrna_hmm._positions` — easy to cross-wire.
- `mitom_pipeline.py:62-110` re-exports ~100 private names (`_LOCUS_BPCOUNT_RE`, `_similarity`, ...) from its
  parts for `batch_mitoz` compatibility — the "module boundary" is nominal; any script can (and tests do) reach into
  privates.
Fix: consolidate on the `mito_common` helpers; keep alias shims where compatibility requires.

---

## 2. Annotation algorithms and logic

### 2.1 MEDIUM — 3'-end rule silently falls back to a stop deep inside the next feature
`mito_intron.py:890-900` (`extend_to_stop`)

```python
if s[c0:c0 + 3] in stops:
    t = abbreviated() if c0 + 3 > limit + STOP_OVERLAP_MAX else None
    return (exons[:-1] + [(a, t or c0 + 3)]), False
```
The docstring promises: "A first complete stop more than STOP_OVERLAP_MAX nt inside the next feature gives way to
an incomplete stop right before it." But when `abbreviated()` finds no in-frame T/TA ending exactly at the limit,
`abbreviated()` returns `None` and `t or c0 + 3` uses the **deep stop anyway** — the written CDS runs arbitrarily far
into the downstream feature. `mito_curate` P3 only curates single-interval CDSs (`simple(f)`), so the overlap is
never repaired downstream and lands in the final GenBank.
Trigger: an intron model whose first in-frame stop lies > 30 nt inside the next feature with no T/TA before it
(e.g. the Glycera-type cox1→cox2 case the comment cites).
Fix: when `abbreviated()` is `None` and the stop is beyond `limit + STOP_OVERLAP_MAX`, end at `limit` as 3'-partial
(stop not determined) instead of `c0 + 3`.

### 2.2 MEDIUM — curation uses the record's genetic code but never writes `/transl_table`
`mito_curate.py:1289-1292` (`record_code`), `mito_curate.py:1364` (`code = a.genetic_code or record_code(rec) or 5`)

The curator translates with the record's majority `/transl_table` (or 5) but the curated output gets no
`/transl_table` on CDS features that lacked one. Everything downstream then treats the record as code 1:
`ncbi_validate.record_gencode` finds no code → `fasta_defline` writes no `[mgcode=]` → table2asn validates with the
standard code → false `InternalStop` ERRORs on TGA-containing genes, and the `.tbl`/final GenBank carry no code at
all (also contradicting `mito_check`'s own `no_transl_table` warning). `mito_core`, `mito_rescue`,
`mito_intron(_ref)` all write `transl_table` per CDS — curation is the gap (standalone `mitom curate` / `import` of
annotators that omit it, e.g. MitoZ-flavoured GenBank).
Trigger: `mitom curate` on a GenBank file whose CDSs have no `/transl_table`.
Fix: when the determined code differs from none, set `f.qualifiers["transl_table"] = [str(code)]` on every CDS that
lacks it (and note it in the report).

### 2.3 LOW-MEDIUM — `mito_trna.py` treats linear FASTA records as circular
`mito_trna.py:805-808`

```python
L = len(rec.seq)
global CIRC_L
CIRC_L = L if rec.annotations.get("topology", "circular") == "circular" else 0
preds = unwrap(preds, L)
```
For FASTA input there is no `topology` annotation, so the default `"circular"` applies: `overlap()` compares
predictions across a non-existent origin (a tRNA predicted at each end of a linear record can be merged into one
group), `unwrap()` extends any `start > end` prediction past the end, and `feature()` then writes an origin-crossing
`join(...)` on a linear record. `topology_from_header` (mito_common) is not consulted.
Trigger: `python mito_trna.py -i linear_contigs.fasta ...` with a tRNA prediction touching the ends.
Fix: use the FASTA `topology=` tag when present, else default to linear (`CIRC_L = 0`) and skip `unwrap`.

### 2.4 LOW — library CDS without `/transl_table` translated with code 1
`mito_intron_ref.py:222`

```python
code = int(f.qualifiers.get("transl_table", ["1" if f.type == "CDS" else "0"])[0] or 0)
```
For the intron-site library (`ref_genes`), a CDS without `/transl_table` is sanity-translated with the **standard
code (1)** — a mitochondrial CDS without the qualifier will mismatch its `/translation` (TGA = stop) and the whole
gene (with its junction sites) is silently dropped from the library. `mito_intron.load_ref_gb` defaults to the
command's genetic code instead, which is the better heuristic.
Fix: default to the record's majority CDS code, else the CLI code — and log dropped entries.

### 2.5 LOW — HSP split threshold off by one vs the documented 30 nt
`mito_intron_ref.py:567` (`SPLIT_GAP = 30  # a subject insertion of >= 30 nt inside one HSP is an intron`) vs
`mito_intron_ref.py:583-584`:

```python
if last_t is not None and t - last_t > SPLIT_GAP + step:
    pieces.append({})
```
The inserted-columns count between two consecutive aligned columns is `t - last_t - step`; the condition fires only
when that is `> SPLIT_GAP`, i.e. an insertion of **>= 31 nt** (blastn) — a 30-nt insertion does not split, contrary
to the comment.
Fix: decide the intended boundary and align code/comment (`>= SPLIT_GAP + step` splits at 30). (Impact is small:
introns of exactly 30 nt are below `MIN_QUERY_INTRON` anyway, so this is a documentation/logic mismatch, not a
calling change.) **(uncertain about intended semantics — hence low)**

### 2.6 LOW — MITOS2 CDS fragments chained in letter order may invert a minus-strand gene **(uncertain)**
`mito_engines.py:140-150` (`_chains`): CDS fragments `cox1-a`, `cox1-b` are chained in letter order. If MITOS2 names
fragments in *genome* order, then on the minus strand the letter order is the reverse of transcript order, and the
resulting `CompoundLocation` yields a reversed `extract()`/translation (5' and 3' exons swapped) for
`--engine mitos2` / `mitom import` of such genes. Whether this fires depends on MITOS2's fragment naming order on the
minus strand, which I could not verify offline **(uncertain)**.
Fix: order chained parts by strand (ascending genome for +, descending for -) or by protein-coordinate evidence;
add a minus-strand split-gene test (e.g. hexacoral nad5).

### 2.7 LOW — `mito_ncr.py` standalone driver: no input validation, case-sensitive suffixes, no cli_guard
`mito_ncr.py:119-136`: a missing input folder or a folder without GenBank files prints
`[mito_ncr] 0 non-coding regions ... -> out` and **exits 0** (silent "success"); the suffix filter
`x.suffix in (".gb", ...)` is case-sensitive (`.GB` missed on Linux, inconsistent with
`mito_common.iter_gb_files`' `.suffix.lower()`); parse errors surface as raw tracebacks (no `cli_guard`, unlike the
sibling tools).
Fix: `sys.exit(cli_guard(main))`, `require_records`, case-insensitive suffixes.

### 2.8 LOW — the core's own outputs are not atomically written (and one without encoding)
`mito_core.py:677-684` (`write_genbank` → plain `SeqIO.write`) and `mito_core.py:783`
(`(out / f"{prefix}.core.tsv").open("w")`, no encoding). A core killed mid-write leaves a truncated `.gbf` (the
pipeline is protected by the exit code; standalone `mito_core.py annotate` is not) — the project has
`mito_common.write_gb` for exactly this.
Fix: use `write_gb` for the `.gbf` and `encoding="utf-8"` for the TSV.

---

## 3. Biological / NCBI conventions

### 3.1 MEDIUM — tunicates get the vertebrate code (2) instead of 13 in the `import --clade auto` path
`mitom_inputs.py:89` (`auto_clade`): `code = genetic_code or (next(iter(codes)) if codes else
PHYLUM_CODE.get(phylum, 5))` with `mito_taxa.py:30` (`PHYLUM_CODE`) and `PHYLUM_SYNONYMS` mapping Tunicata →
Chordata. For a GenBank record with taxonomy but **no CDS `/transl_table`**, `mitom import --clade auto` selects
code 2 for a tunicate; NCBI uses code 13 (Ascidian Mitochondrial) there — AGA/AGG become stops and proteins are
truncated/fragmented in annotation. Two sibling paths handle this correctly: `mitom_groups.decide` uses
`mito_classify.code_of_lineage` (ref_taxa votes, Tunicata → 13), and `mito_tbl.record_code` has
`SUBPHYLUM_CODE = {"Tunicata": 13, ...}`. Three different lineage→code implementations disagree. (Related
robustness note: `code_of_lineage`→`ref_taxa()` raises `FileNotFoundError` if `models/ref_taxa.tsv.gz` is absent —
`decide()` has no guard there, so a broken install crashes `mitom run --clade auto` on taxonomy-only input.)
Fix: one shared helper (taxonomy → code, subphylum-aware, with a missing-model fallback to `PHYLUM_CODE`) used by
`auto_clade`, `decide` and `mito_tbl`.

### 3.2 MEDIUM — atp9 is silently dropped from the MITOS2 engine and `mitom import`
`mito_engines.py:237-243` (`mitos2_record`)

```python
elif f.type in ("CDS", "rRNA") and name in mito_names.CORE_NAMES:
    ...
else:
    continue
```
MITOS2 annotates `atp9` in sponges and some cnidarians; `CORE_NAMES` lacks it (`mito_names.PRODUCTS` has
`"atp9": "ATP synthase F0 subunit 9"`, and `mito_intron_ref` can write atp9 models). With `--engine mitos2` or when
importing MITOS2 results, atp9 features are discarded without a note. For a phylum where it is a genuine, often
annotated gene this is a silent false negative in the deliverable.
Fix: include `atp9` in the accepted names (product string exists) and pass it through as a CDS.

### 3.3 LOW — `/translation` of intron-added CDS does not start with M for alternative start codons
`mito_intron.py:1189-1194` (`apply_model`: `prot = str(nt[:...].translate(table=gencode, to_stop=True))`) and
`mito_intron_ref.py:1455-1460,1883` (`q["translation"] = [m.protein]`). A CDS starting with GTG/TTG/ATT is written
with its literal first residue (V/L/I) while NCBI's convention (and `mito_curate.translation`, lines 946-948) writes
M for any valid start codon of the table. In `mitom run` the curation step rewrites translations, so only standalone
`mitom intron` outputs are affected.
Fix: apply the same rule (`M` + rest) when the first codon is a start codon and the CDS is not 5'-partial.

### 3.4 LOW — a zero-length site at the sequence boundary becomes an inverted interval
`mito_tbl.py:178-181` (`clean_record`): the `100^101`-style site fix requires `0 < int(loc.start) < n`; a
between-bases site exactly at position 0 or n is left as zero-length, and `feature_lines` then emits
`start > end` (e.g. `16570..16569`) — a Bad interval for table2asn.
Fix: clamp such sites into the sequence (`max(0, start-1)`, `min(n, start+1)`).

### 3.5 LOW — single-position `/transl_except` values are dropped by curation
`mito_curate.py:195` (`_EXCEPT_RE = re.compile(r"pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)")`) matches only
ranged `pos:N..M`. A legal single-position exception (`pos:123,aa:Trp`) does not match, so `internal_excepts`
silently discards it (single-position `aa:TERM` is fine — it is re-added by the caller — but a non-TERM
single-position exception, e.g. a one-codon read-through, is lost from the output).
Fix: make the range part optional (`(\d+)(?:\.\.(\d+))?`), as `mito_compare._lift_except` already does.

---

## 4. Debug functionality

Verification summary of the required checks (what I traced end-to-end):
- `mitom debug` parsing of run.log (steps, `[CMD]`, `[RECORD]`, `[DONE]`, `[ERROR]`, old 1.0.0 format) — works as
  documented; exit codes 0/1/2 correct.
- `--debug` propagation: `mitom` launcher sets `MITOM_DEBUG=1` before `execv`, `mitom_pipeline.main` re-sets it
  (line 1565-1569) before any child runs → annotation core subprocesses, intron subprocesses, Snakemake workers and
  fork workers all inherit it; `setup_root_logger(debug=True)` puts DEBUG lines of every module in run.log. OK.
- Child-process stderr: the annotation core's stdout/stderr are captured to `logs/core.std{out,err}.log`
  (`mitom_stage.run_cmd_capture`), the last stderr line goes into the `[RECORD]` event, and `run_tool` logs
  `[CMD rc=...]` + stderr tail to run.log (also in non-debug mode). OK — with the gaps below.
- `reports/failed_logs`: written by `cleanup_intermediates` — but only for a subset of failure kinds (4.1).

### 4.1 MEDIUM — logs of topology-phase and round-2 failures are deleted before `mitom debug` can read them
`mitom_output.py:325-334` (`cleanup_intermediates` copies logs only for `summary.all_problem_samples()` =
annotation_failed + final_write_failed + no_target_gene + metadata_transfer_failed + metadata_key_mismatch +
reorient_failed). Two failure kinds emit `[RECORD]` lines pointing at their sample log but are **not** in that set:

- `mitom_topology.py:170-172` — `record_event(si.sample_name, f"annotation as {topology} failed: ...")` for the
  `02.anno_topology` re-annotation (its workdir is deleted by the same cleanup);
- `mitom_pipeline.py:668` — `record_event(si.sample_name, "round-2 annotation failed: ...")` for
  `03.anno_round2_reoriented` failures (round-1 annotation is kept, so the sample is not a "problem sample").

With a default run (`--keep_intermediate no`), `mitom debug`'s `record_logs` then finds neither a `failed_logs`
copy nor the live workdir — the actual stderr of the failed attempt is unrecoverable, contradicting
"core logs of problem samples are always kept (reports/failed_logs/)" (`mitom_pipeline.py:1236` help text).
Fix: add samples named in any `[RECORD]` "annotation failed" event (or the 02/03 failure lists) to the
`problem_samples` set passed to `cleanup_intermediates`.

### 4.2 MEDIUM — `mitom debug` cannot resolve sample-log paths on Windows, nor the `02` phase anywhere
`mito_debug.py:216` and `mito_debug.py:230`:

```python
m = re.search(r"sample log: (\S+?/(0[13]\.[\w.]+)/([^/\s]+))/sample\.log", text)
...
m = re.search(r"/0[13]\.[\w.]+/([^/\s]+)/sample\.log", "\n".join(events))
```
(a) `annotate_one` logs `sample log: {workdir / 'sample.log'}` — a `pathlib` rendering, i.e. **backslashes on
Windows**, which these forward-slash regexes never match; (b) the phase class `0[13]` excludes `02.anno_topology`,
so even on POSIX a topology-phase sample log is not located; (c) `relocated()` (`mito_debug.py:263-274`) rewrites
work paths only in forward-slash form, so the "run the failed program by hand" path is stale on Windows too.
Fix: accept `[\\/]` as separators, include `02`, and normalize with `PurePath`-style parts instead of regex on raw
strings.

### 4.3 MEDIUM — failed records of the curate / intron / rescue / rrna_hmm steps never reach `mitom debug`'s records section
The debug "records" section is driven by `[RECORD]` lines plus the four `FAILED_LISTS` (`mito_debug.py:50-51`).
These per-record failures are written only into their own TSV reports and are invisible to `mitom debug`:

- `mito_curate.py:1368-1370` — a record that crashes curation becomes a `curation_report.tsv` row
  ("left unchanged") with **nothing logged** to run.log;
- `mito_intron.py:1533-1535` — `error, record left unchanged: ...` rows (same for `mito_intron_ref._annotate_file`);
- `mito_rescue.py:201-204` — `status: failed: ...` rows in `rescue_report.tsv`;
- `mito_rrna_hmm.py:370-372,413-414` — `failed: ...` rows in `rrna_hmm_report.tsv`.

A user asking "`mitom debug`: which records failed?" gets an incomplete answer whenever a record fails in these
steps. Fix: emit `record_event(...)` for these (or add the four reports to `FAILED_LISTS`-style parsing — but the
`[RECORD]` line is the mechanism the parser already reads everywhere).

### 4.4 MEDIUM — group runs (`--clade auto` over several phyla): no per-record log enrichment, and group failed_logs overwrite each other
`mito_debug.py:390-400`: group logs (`reports/groups/<group>.run.log`) are parsed with `parse(g, log.parent)` but
`add_output(gr)` is never called for them, so the "first error in reports/failed_logs/..." lines that the single-run
report produces per failed record are missing for groups.
Additionally `mitom_groups.py:300-303` merges each group's `reports/failed_logs/*` into one
`reports/failed_logs/` with unchanged names (`dst = t`): identical sample names across groups (possible with
FASTA input, where names repeat) silently overwrite — the last group copied wins.
Fix: call `add_output(gr)` per group (with the group's own `out/`) and prefix merged failed_logs with the group name.

### 4.5 LOW — some external calls stay invisible in the `[CMD]` stream
`mito_classify.py:125-132` (`_blastx`, raw `subprocess.run`) and `mito_nuc.py:142-143` (see 1.8): these programs
never produce `[CMD rc=...]` lines, so `mitom debug`'s "failed programs" section cannot show their exit code /
stderr tail; a classify blastx failure appears only as a note string.
Fix: route them through `run_tool` (or log an equivalent `[CMD]` line).

### 4.6 LOW — sample logs capped at INFO even in `--debug`; logger objects accumulate
`mitom_log.py:167-182`: `get_sample_logger` hard-sets `lg.setLevel(logging.INFO)` — with `--debug` the per-sample
`sample.log` still has no DEBUG lines (the core's own stderr log does, so impact is limited). Also each
`get_sample_logger` creates a new named logger that stays in the logging registry after
`close_sample_logger` (only handlers are removed) — thousands of samples grow the registry for the process's
lifetime.
Fix: honour `debug_mode()` for the level; `logging.Logger.manager.loggerDict.pop(name, None)` on close (or reuse one
logger with changing handlers).

---

## Fix checklist (ordered)

1. **`lock_output` (1.1)**: replace `os.kill(pid, 0)` with a Windows-safe liveness probe — *stops the lock from
   killing healthy runs on Windows*.
2. **fork→spawn/fallback (1.2)**: platform-choose the multiprocessing context in `mito_intron.py:1612`,
   `mito_intron_ref.py:445/2014`, `mitom_steps.py:388` — *fixes `--jobs > 1` crashes and the silently missing maps on
   Windows*.
3. **`os.devnull` (1.3)**: `mito_rrna_hmm.py:242` — *restores the rRNA step on Windows*.
4. **`write_gb` everywhere (1.5)**: route the listed `SeqIO.write(..., "genbank")` sites through
   `mito_common.write_gb` — *fixes `/pseudo=""` and makes the writes atomic*.
5. **`.tbl` value-less qualifiers (1.6)**: bare `\t\t\tpseudo` for `None`/`""` in `ncbi_validate.feature_lines`.
6. **`mito_intron` prepared-input root (1.4)**: compute `src_root` from `src_in`, as `mito_intron_ref` does.
7. **`extend_to_stop` 3' fallback (2.1)**: end 3'-partial at the limit when no T/TA, never the deep stop.
8. **write `/transl_table` in curation (2.2)**: set the determined code on CDSs lacking it.
9. **`transfer_metadata` (1.7)**: log the exception; exit 1 when any file failed.
10. **failed_logs coverage (4.1)** + **debug regexes (4.2)**: keep 02/03 failure logs; accept Windows separators and
    the `02` phase in `mito_debug`.
11. **record_event for curate/intron/rescue/rrna_hmm record failures (4.3)**; per-group `add_output` and
    group-prefixed failed_logs (4.4).
12. **`mito_nuc.run` (1.8)** and **`mito_classify._blastx` (4.5)**: use `run_tool`/`run_or_raise` so stderr and
    `[CMD]` lines are kept.
13. **Genetic-code helper (3.1)**: one subphylum-aware taxonomy→code function for `auto_clade`, `decide`,
    `mito_tbl`; guard `ref_taxa()` when the model file is missing.
14. **atp9 passthrough (3.2)** in `mito_engines.mitos2_record`.
15. **Sample-logger `try/finally` (1.9)**; **`mito_trna` linear default (2.3)**; **single-position `/transl_except`
    (3.5)**; **`mito_ncr` driver hardening (2.7)**; encoding/atomic-write cleanups (1.10, 2.8, 1.16).

Test gaps worth adding after the fixes: a Windows CI job (would have caught 1.1-1.3 and 4.2), an end-to-end test
that `mitom debug` finds a 02/round-2 failure log (4.1), a `/pseudo` round-trip test for every GenBank-writing
command (1.5), and a table2asn golden test on a record with a pseudogene (1.6, 3.3).
