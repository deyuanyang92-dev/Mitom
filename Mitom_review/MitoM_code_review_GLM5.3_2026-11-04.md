# MitoM code review

- Reviewed by: **GLM 5.3** (independent review, no other review document was read or reused)
- Date: **2026-11-04**
- Mode: **read-only** — no source code, tests, or configuration was changed; the only file written is this document.
- Scope: production code in `P:\Mitom` (36 top-level `.py` files, ~19.5k lines), per the requested focus:
  engineering standards; annotation algorithms (`mito_core.py`, `mito_curate.py`, `mito_trna.py`, `trna_ml.py`,
  `mito_rrna_hmm.py`, `mito_intron.py`, `mito_intron_ref.py`, `mito_nuc.py`, `mito_ncr.py`, `mito_engines.py`,
  `mito_classify.py`); biological / NCBI conventions (`mito_tbl.py`, `ncbi_validate.py`, `mito_names.py`,
  genetic codes); and a dedicated end-to-end review of the debug functionality (`mito_debug.py`, `mitom_log.py`,
  `--debug` flags, pipeline error handling). `tests/`, `benchmark/`, `knowledge/` were not reviewed except where a
  production question required looking at a test fixture.
- Conventions used below: every problem has **severity (high/medium/low)**, **file + line numbers**, a **short code
  quote**, **why it is wrong**, a **realistic triggering input**, and a **concrete potential fix (not implemented)**.
  Items whose biological or platform consequence I could not fully verify are marked **不确定**; per instructions
  they are not upgraded in severity.

---

## 1. Software engineering standards

### 1.1 HIGH — `multiprocessing.get_context("fork")` in production code breaks on Windows

- Files / lines:
  - `mitom_steps.py:388` — `with mp.get_context("fork").Pool(max(1, min(len(jobs), ...))) as pool:`
  - `mito_intron.py:1613` — `with mp.get_context("fork").Pool(a.jobs) as pool:      # workers inherit _JOB (no pickling)`
  - `mito_intron_ref.py:446` — `with mp.get_context("fork").Pool(jobs) as pool:`
  - `mito_intron_ref.py:2015` — `with mp.get_context("fork").Pool(a.jobs) as pool:`
- Why wrong: on Windows `multiprocessing.get_context("fork")` raises
  `ValueError: cannot find context for 'fork'`. Consequences per site:
  - `run_plots` (`mitom_steps.py:387-392`): the `except Exception` swallows it, so **genome maps are never produced on
    Windows** (default `--plot yes`); the step is logged as `[PLOT] failed` and the run "succeeds" without maps.
  - `mito_intron.py` / `mito_intron_ref.py` with `--jobs > 1`: the exception aborts the whole script, so
    `mitom intron` and the pipeline's intron step produce **no introns at all** on Windows whenever more than one job
    is requested (pipeline default `--jobs` is `max(1, CPUs // threads)`, i.e. > 1 on multi-core machines).
  - The design also depends on fork semantics: workers read the module globals `_JOB` set in the parent
    (`mito_intron.py:1518`, `mito_intron_ref.py:1909`), so simply switching to `"spawn"` would break state passing.
- Realistic trigger: `mitom run -i x.gb -o out` on Windows (maps lost); `mitom run` on Windows with
  `--intron_fix auto` and default `--jobs` (no introns, exit non-zero for the step).
- Potential fix: choose the context per platform (`mp.get_context("fork" if hasattr(os, "fork") else "spawn")`)
  and pass the shared state explicitly via `Pool(initializer=..., initargs=...)` or make the worker arguments fully
  picklable (`functools.partial` / small payload dicts), so both fork and spawn work. Alternatively document Windows
  as unsupported for these steps and fail loudly instead of silently.

### 1.2 HIGH — nhmmer is sent to the POSIX-only path `/dev/null`

- File / line: `mito_rrna_hmm.py:242`
  ```python
  p = run_tool([exe, "--dna", "--cpu", "2", "-E", str(EVALUE), "--tblout", str(tbl), "-o", "/dev/null",
                str(hmm), str(fa)])
  ```
- Why wrong: `models/rrna_hmm` ships with the install (64 files present), and `--rrna_method hmm` is the default, so
  `annotate_rrnas` calls `nhmmer` for every record whose phylum has a model. On Windows there is no `/dev/null`
  (that is `NUL`), so a native HMMER build fails to open the output file, `run_tool` returns a non-zero exit and
  `nhmmer()` raises `RuntimeError(f"nhmmer failed (exit {p.returncode})...")`
  (`mito_rrna_hmm.py:246-248`). In the annotation core this exception is **not caught per record**: it propagates out
  of `mito_core.annotate` → the engine subprocess exits non-zero → **every record fails annotation** on Windows
  whenever a phylum model matches. (In the later rescue step `mito_rrna_hmm.search_record` it is caught and becomes
  a "failed: ..." report row — inconsistent handling of the same root cause.)
- Realistic trigger: `mitom run -i genome.fa --clade Mollusca` on Windows (or any run where `--phylum`/taxonomy
  matches a model in `models/rrna_hmm`).
- Potential fix: write nhmmer's human-readable output to a temp file (it is discarded anyway) or use
  `os.devnull`; e.g. `"-o", str(tmp / "nhmmer.out")` or `"-o", os.devnull`.

### 1.3 HIGH — stale-lock probe `os.kill(pid, 0)` kills the other run on Windows

- File / line: `mitom_pipeline.py:1471-1479`
  ```python
  other = int(lock.read_text().split()[0])
  os.kill(other, 0)                      # still running?
  ```
- Why wrong: on POSIX `os.kill(pid, 0)` is a harmless liveness probe. On Windows, `os.kill` has no "signal 0"
  semantics: any value other than the console-event signals is passed to `TerminateProcess` as the exit code — the
  documented Windows behavior is "any other value for sig will cause the process to be unconditionally killed".
  So on Windows, when run B finds the lock of a still-running run A, it **terminates run A** (exit code 0), then takes
  over the lock and writes into the same output folder. Run A dies silently mid-annotation.
- Realistic trigger: two `mitom run` commands (or a re-run started while the previous one is still going) with the
  same `-o` folder on Windows.
- Potential fix: probe liveness portably — e.g. on Windows use `ctypes` (`OpenProcess`/`GetExitCodeProcess`), or
  simply trust the lock with an age check (`mtime` older than N hours → stale), or use `os.kill` only when
  `hasattr(os, "fork")`.

### 1.4 MEDIUM — locale-dependent text encoding on all GenBank/FASTA writes (and one reader)

- Files / lines (representative):
  - `mito_common.py:126` — `with os.fdopen(fd, "w") as h:` in the atomic `write_gb` (no `encoding=`, default newline)
  - `mito_common.py:141` — `with open(path) as fh:` in `read_fasta` (no `encoding=`)
  - Biopython `SeqIO.write(records, str(path), "genbank")` opens the file with the **locale** encoding; production
    sites include `mito_core.py:782`, `mito_curate.py:1376`, `mito_trna.py:862`, `mito_ncr.py:133`,
    `mito_engines.py:361`, `mitom_stage.py:141`, `mitom_orient.py:295/297`, `mitom_metadata.py:470/513`,
    `mitom_output.py:148`, `mitom_pipeline.py:346/1036/1067/1069`, `mito_nuc.py:1718/1732`.
- Why wrong: on Windows the default encoding is the ANSI code page (cp1252/cp936/...), not UTF-8. MitoM deliberately
  supports Latin-1 / BOM inputs (`mitom_inputs.normalize_input` converts them), so records with non-ASCII
  organism names or reference titles flow through the pipeline. On a cp936 (Chinese) locale, writing the final
  GenBank with an organism like `Bactrilum caf\u00e9ensis` raises `UnicodeEncodeError` → the record fails at the
  "final GenBank not written" stage; on cp1252 it silently writes wrong bytes (mojibake). `read_fasta`
  (`mito_common.py:141`) reading a UTF-8 FASTA under cp1252/cp936 likewise garbles or crashes on headers, while
  `count_seq_records` (`mito_common.py:570`) reads the same files as UTF-8 — inconsistent. As a secondary effect,
  text-mode writes without `newline=""` produce CRLF GenBank files on Windows (Biopython re-parses them fine; NCBI
  tools are tolerant, but it is inconsistent with the LF-locked TSV writers).
- Realistic trigger: a Geneious/ENA GenBank record with an accented organism or author name, on a Windows machine
  with a non-UTF-8 locale.
- Potential fix: pass `encoding="utf-8"` (or `"utf-8", newline="\n"`) at every GenBank/FASTA write — at minimum in
  `write_gb` (`os.fdopen(fd, "w", encoding="utf-8", newline="\n")`) and add `encoding="utf-8", errors="replace"` to
  `read_fasta`; or centralize all GenBank writes through `write_gb`.

### 1.5 MEDIUM — `tarfile.extractall(filter="data")` needs a newer Python than the declared 3.8 minimum

- File / line: `trna_ml.py:109` — `t.extractall(tmp, filter="data")`
- Why wrong: `mitom_pipeline.check_environment` (`mitom_pipeline.py:124-125`) only requires Python 3.8
  (`if (major, minor) < (3, 8): raise RuntimeError(...)`). The `filter=` keyword of `tarfile.extractall` exists only
  in Python ≥ 3.12 (and the back-ported security releases 3.11.4 / 3.10.12 / 3.9.17 / 3.8.17). On an older 3.8/3.9
  interpreter the call raises `TypeError: extractall() got an unexpected keyword argument 'filter'`, so MitoM-ML's
  clade models can never unpack and every record falls back to MiTFi tRNAs — silently different results on different
  Python patch versions.
- Realistic trigger: a conda environment with Python 3.9.5 and the shipped `models/trna_ml/clades.tar.gz` (first
  MitoM-ML call: `clade_dir()` → `extractall`).
- Potential fix: guard the kwarg —
  `try: t.extractall(tmp, filter="data")\nexcept TypeError: t.extractall(tmp)` — and/or raise the declared minimum
  Python to a version that has the parameter.

### 1.6 MEDIUM — annotation-stage subprocesses have no timeout; a hung engine stalls the run forever

- Files / lines:
  - `mitom_stage.py:36-46` (`run_cmd_capture`): `p = run_tool(cmd, cwd=cwd)` — `run_tool`'s `timeout` default is
    `None`.
  - `mitom_snakemake.py` worker (`render_snakemake_worker_script`, ~line 88): `p = subprocess.run(cmd, ...,
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)` — no timeout, and no `run_tool` logging either.
- Why wrong: `run_or_raise` already defines `TOOL_TIMEOUT = 3600` for the core's own programs
  (`mito_common.py:466-475`), but the *engine* subprocess (a whole `mito_core.py` / MitoZ / MITOS2 run) is started
  without any bound. MitoZ and MITOS2 are third-party pipelines that can hang (network, java, blast deadlocks); with
  the threads backend the `as_completed` loop waits forever, and the "a hung program never stalls a run" promise is
  broken at the stage level. The Snakemake worker has the same gap (and no `[CMD]` entry at all, see 4.2).
- Realistic trigger: `--engine mitoz` where MitoZ's java step hangs on one sample.
- Potential fix: pass a generous stage-level timeout (e.g. 24 * TOOL_TIMEOUT) to `run_tool` in `run_cmd_capture`,
  mark the sample failed with a `record_event` on `TimeoutExpired`, and give the Snakemake worker the same
  `subprocess.run(..., timeout=...)`.

### 1.7 MEDIUM — `workflow_transfer_metadata` swallows the exception and its message completely

- File / line: `mitom_pipeline.py:946-950`
  ```python
          except Exception:
              if args.overwrite_new_body == "no":
                  safe_copy(bf, out_file)
              n_fail += 1
              log_rows.append(f"{bf}\t{matched or 'NA'}\t{out_file}\tfail\t0\t0")
  ```
- Why wrong: this is the only failure path for a body file in `transfer_metadata`, and the exception text is
  discarded — not logged, not in the TSV. The user sees status `fail` with no way to know why (parse error? bad
  metadata record? permission?), and there is no `[RECORD]`/`[ERROR]` line for `mitom debug` to pick up. Note the
  sibling workflow `workflow_replace_organism` right below does it correctly: `except Exception as exc: ...
  logging.error(f"[FAIL] {bf.name}: {exc}")` (`mitom_pipeline.py:1004-1010`).
- Realistic trigger: any exception while transferring headers into one of many files (e.g. a GenBank parsing error
  in the metadata file).
- Potential fix: `except Exception as exc:` plus `logging.error(f"[FAIL] {bf}: {exc}")`, and append the reason to the
  TSV row.

### 1.8 LOW — sample-logger `FileHandler` leaks on worker exception

- File / line: `mitom_stage.py:147-177` (`annotate_one`): `lg = get_sample_logger(log_name, ...)` at the top,
  `close_sample_logger(log_name)` only on the two normal returns (lines 195 and 199).
- Why wrong: if `run_cmd_capture`/`engine_cmd` raises (e.g. `OSError` for a missing engine, or a future
  `TimeoutExpired` per 1.6), the exception propagates to `run_annotation_stage`'s `except` (line 254-258) which
  records the event — but `close_sample_logger` never runs for that sample, leaving the `FileHandler` (and its open
  file) attached to a process-global named logger until process exit. A batch with hundreds of failing samples holds
  hundreds of open handles.
- Realistic trigger: `--engine mitoz` with `mitoz` missing from PATH mid-run, or any per-sample `OSError`.
- Potential fix: wrap the body of `annotate_one` in `try/finally`, or close the logger in the `except` branch of
  `run_annotation_stage`.

### 1.9 LOW — non-atomic GenBank writes in modules that already import the atomic writer

- Files / lines:
  - `mito_curate.py:1376` — `SeqIO.write(recs, str(dest), "genbank")      # always: N1 names ...`
  - `mito_trna.py:862` — `SeqIO.write(recs_out, str(out / sub / f"{f.stem}.gbf"), "genbank")`
  - `mito_ncr.py:133` — `SeqIO.write(recs, str(dest), "genbank")`
  - `mito_engines.py:361`, `mitom_stage.py:141`, `mitom_orient.py:295/297`, `mitom_pipeline.py:346/1036`.
- Why wrong: `mito_common.write_gb` exists precisely so "a failure (or a killed run) never leaves a half-written
  file" (`mito_common.py:114-134`), and `mito_intron.py`/`mito_intron_ref.py`/`mito_rrna_hmm.py` do use it. The
  modules above write final/curated GenBank files directly with `SeqIO.write`, so a kill or disk-full during these
  steps leaves a truncated `.gbf` that downstream steps or the user may read as a complete record (the codebase
  itself acknowledges this failure mode: `truncated_genbank` in `mito_common.py:580-600` exists to detect it).
- Realistic trigger: Ctrl-C or out-of-space during `mitom curate` on a large folder.
- Potential fix: route these writes through `write_gb` (it also keeps the value-less-qualifier fix, 3.1).

### 1.10 LOW — `mitom` launcher relies on `os.execv`, which has Windows quirks (不确定 on exact effect)

- File / line: `mitom` (launcher), last block of `main()`: `os.execv(py, [py, str(HERE / script), *pre, *args])`.
- Why wrong: on Windows the `_execv` family is implemented by spawning a new process and terminating the caller;
  exit-status propagation to the waiting shell is a known weakness of this family on Windows (the new program runs
  under a new PID), so `%ERRORLEVEL%`/`$LASTEXITCODE` after `mitom ...` may not reflect the real exit code of the
  pipeline — which matters because MitoM's contract is "0 done, 1 failure, 2 usage" and `mitom debug`'s own exit
  codes. I could not execute this to verify the exact behavior on this machine (不确定); flagging because the whole
  CLI's exit-code contract rides on it.
- Realistic trigger: `mitom run ...` in a Windows script that checks the exit code.
- Potential fix: on Windows use `subprocess.run` + `sys.exit(rc)` instead of `os.execv` (keep `execv` on POSIX), or
  verify and document the behavior.

---

## 2. Annotation algorithms and logic

### 2.1 MEDIUM — `mito_curate` rule O2 (`origin_cut`) is inverted for minus-strand origin-spanning joins

- File / line: `mito_curate.py:655-663`
  ```python
          head, tail = f.location.parts                  # transcript order: the part before the origin first
          if len(tail) > MAX_ORIGIN_EXT or len(head) < 60:
              return
          ...
          f.location = FeatureLocation(head.start, head.end, head.strand)
  ```
- Why wrong: the comment assumes parts[0] is "the part before the origin". That is only true for **plus**-strand
  joins. The core writes minus-strand origin-spanning joins in transcript order as
  `CompoundLocation([last, first])` = `[[0, tail), [s, L)]` (`mito_core.py:668-671`), so for a minus-strand gene
  parts[0] is the **origin-side 5' piece** and parts[1] is the main body that ends at the sequence end.
  Consequences for a minus-strand join:
  1. Common case (main body > 300 nt): `len(tail) > MAX_ORIGIN_EXT` is true because `tail` is the *main body* → O2
     silently does nothing; the documented NoStop repair never runs for minus-strand genes, and the join is later
     only marked 3'-partial by `origin_join_check` (P6) although its 3' end is not at the origin.
  2. Dangerous case (main body ≤ 300 nt and origin-side piece ≥ 60 nt — small genes such as atp8/nad4l running into
     the cox1-first rotation): the rule **keeps `head` (the origin-side fragment) and deletes the gene's main body**,
     replacing a real ~250-nt CDS with a ~70-nt fragment at position 1.
- Realistic trigger: `mitom run` output where the gene preceding cox1 across the origin lies on the minus strand
  and the core joins it across the origin without an in-frame stop (the exact scenario the rule was written for,
  "annelid nad2 running 8-23 nt into cox1", can occur on either strand).
- Potential fix: select the parts by geometry instead of by list order, which is strand-correct for both:
  keep the part whose `int(p.end) == self.L` (the pre-origin body), drop the part whose `int(p.start) == 0`, and
  apply the guard to the dropped part (`len(dropped) <= MAX_ORIGIN_EXT and len(kept) >= 60`).

### 2.2 LOW — `mito_trna.group()` compares each prediction only with the first member of a group

- File / line: `mito_trna.py:264-269`
  ```python
      for p in sorted(preds, key=lambda x: (x.start, x.end)):
          for g in groups:
              q = g[0]
              ov = overlap(p, q)
              ... if p.strand == q.strand and same and ov >= 0.5 * min(p.length, q.length):
                  g.append(p)
  ```
- Why wrong: `q = g[0]` is the earliest-starting prediction of the group. A prediction that overlaps a *later*
  member of the group by ≥ 50 % but not the first member fails the test and starts a new group. The same tRNA locus
  then appears as two groups, both of which can be selected (`--select mitos2_or_support`: both get
  `support >= min_support`), producing duplicate tRNA features at one locus; with `--select mitom` the same family
  can be written twice. Chain-grouping (compare against all members, or merge groups when a member matches) is the
  usual fix.
- Realistic trigger: two tools call the same tRNA with starts > 50 % of the length apart (common for degenerate
  mitochondrial tRNAs), and a third tool agrees with the second boundary — three predictions, one true locus, two
  groups.
- Potential fix: test `p` against every member of the group (or against the union span of the group), and when `p`
  matches two existing groups, merge them.

### 2.3 LOW — `mito_trna.from_input()` takes the coordinate hull of a multi-part input tRNA

- File / line: `mito_trna.py:232-237`
  ```python
      out.append(Pred("input", int(f.location.start) + 1, int(f.location.end), f.location.strand or 1, aa, ...))
  ```
- Why wrong: `location.start/end` is the hull of a `CompoundLocation`. An input tRNA written as a join across the
  origin (or a trans-spliced piece) becomes a `Pred` spanning (nearly) `1..L`, which can then overlap every group,
  poison `conflicts()` (`mito_trna.py:276-284` flags overlapping groups as conflicts), and, if selected, write a
  whole-genome tRNA feature. `parse_laslett`/`parse_mitos_bed` handle origin-crossing calls explicitly
  (start > end → `unwrap`); the input-annotation reader does not.
- Realistic trigger: `mitom trna -i annotated.gb` where the input GenBank has a tRNA as
  `join(16001..16560,1..60)` (annotated by another tool on a circular record).
- Potential fix: iterate `f.location.parts` per strand (as `mito_trna.py:812` does for output) and emit one `Pred`
  per part with the wrap marker, or extend the end by `+ L` for a two-part origin join, mirroring `unwrap()`.

---

## 3. Biological / NCBI conventions

### 3.1 MEDIUM — value-less qualifiers (`/pseudo`, `/trans_splicing`) are written into the .tbl with a value

- File / line: `ncbi_validate.py:54-63` (`feature_lines`, used by both `mitom_tbl` and the pipeline's submission step)
  ```python
              for v in vals:
                  if k == "note" and INTERNAL_NOTE.search(str(v)):
                      continue
                  out.append(f"\t\t\t{k}\t{' '.join(str(v).split())}")
  ```
- Why wrong: `mito_curate` P10 sets `x.qualifiers["pseudo"] = [None]` (`mito_curate.py:886, 898`) and
  `mito_intron_ref._apply_model` sets `trans_splicing` to `[None]` (`mito_intron_ref.py:1887-1888`); Biopython
  correctly writes them value-less (`/pseudo`) in GenBank, and reads them back as `[""]`
  (`mito_common.py:107-109` documents this). But in the 5-column feature table the value column is emitted
  unconditionally:
  - records parsed from GenBank (the normal pipeline path, `05.curated/` → `run_submission`) give `v = ""` → the row
    is `pseudo\t` (trailing tab, empty value);
  - in-memory records in `mito_curate --validate` give `v = None` → the literal text `pseudo\tNone`.
  The NCBI 5-column spec has no value for these qualifiers; a stray `None`/empty fifth column can make table2asn
  misparse the qualifier (exact behavior 不确定 — it may also be tolerated; severity kept at medium because the
  same code path writes every submission `.tbl`).
- Realistic trigger: any record where P10 marked a CDS `/pseudo` (internal stop) — the submission `.tbl` and the
  before/after curation validation both hit it.
- Potential fix: `if v in (None, ""): out.append(f"\t\t\t{k}") else: out.append(f"\t\t\t{k}\t{...}")`.

### 3.2 MEDIUM — an explicit `--clade Chordata` forces genetic code 2 on tunicates and lancelets

- Files / lines:
  - `mito_taxa.py:14` — `CLADE_CODES = {"Chordata": 2, "Echinodermata": 9, "Porifera-sponges": 4}`
  - `mitom_pipeline.py:307-309` — `args.genetic_code = CLADE_CODES.get(args.clade, 5)` when `--genetic_code` is absent.
- Why wrong: the Chordata clade includes Tunicata (mitochondrial code **13**) and Cephalochordata (code **5**) — the
  codebase itself encodes this knowledge in `mito_tbl.SUBPHYLUM_CODE = {"Tunicata": 13, "Cephalochordata": 5}`
  (`mito_tbl.py:141`). With an explicit `--clade Chordata`, all records are annotated and curated under code 2, in
  which AGA/AGG are stop codons: tunicate CDS get internal stops and are written `/pseudo` by P10
  (`mito_curate.py:880-901`), starts/stops are misplaced, and `mito_curate` runs with the wrong codon table. Only
  GenBank input with `/transl_table` gets a `[CODE] ... differs from --genetic_code 2 — check` warning
  (`mitom_inputs.py:845-849`); FASTA input of a tunicate gets no warning at all.
- Realistic trigger: `mitom run -i tunicate_mitos.gb --clade Chordata -o out` (a natural thing to do for a
  non-vertebrate chordate), or the same with FASTA input.
- Potential fix: when `--clade` is Chordata and the record lineage (or `--clade auto`) contains Tunicata /
  Cephalochordata / Hemichordata, use the subphylum code (as `mito_tbl.record_code` already does for
  submission files), or make code 2 require an explicit `--genetic_code 2` for non-vertebrate chordates and warn
  otherwise.

### 3.3 LOW — `ncbi_validate.fasta_defline` does not ASCII-clean the organism, unlike `mito_tbl.defline`

- File / line: `ncbi_validate.py:66-75`
  ```python
      org = rec.annotations.get("organism", "unknown organism")
      ...
      mods = [f"[organism={org}]"] + ...
      (d / "rec.fsa").write_text(fasta_defline(...) + "\n" + str(rec.seq) + "\n")
  ```
- Why wrong: `mito_tbl.ascii_text` exists because NCBI accepts ASCII only
  ("text is made ASCII (NCBI rejects anything else, SEQ_INST.NonAsciiAsn)", `mito_tbl.py:45-47, 166-168`), and the
  `mito_tbl` defline applies it. `ncbi_validate.validate_record` writes the raw organism into `rec.fsa` and opens
  the file with the locale encoding (compound of 1.4): a non-ASCII organism either crashes the write (cp936) or
  produces a `.fsa` that table2asn flags / mangles. Also, when a record has no header organism (core output),
  `[organism=unknown organism]` is submitted, which `mito_tbl` explicitly treats as a placeholder requiring user
  action; here it silently passes validation.
- Realistic trigger: `mitom validate -i final.gbf` on a record whose GenBank header carries an accented organism
  name (transferred from the original record in step "write").
- Potential fix: reuse `mito_tbl.ascii_text` (or inline the NFKD fold) on `org` in `fasta_defline`, and write the
  `.fsa` with `encoding="utf-8"` / plain ASCII.

---

## 4. Debug functionality (end-to-end review of `mito_debug.py`, `mitom_log.py`, `--debug`, error capture)

Verification results for the questions asked (verified against the code, not the docs):

- **Does `mitom debug` locate the failing step/record/program from run.log?** Yes. `parse()`
  (`mito_debug.py:120-200`) reads `[STEP n/N key]`, `[CMD rc=...]` (+ stderr tail lines after `    | `), `[RECORD ...]`,
  `[ERROR]`, `[DONE] exit=...`, `[ENV]`/missing-program lines, the `reports/*_failed*.txt` lists and
  `ncbi_submission/validation_messages.tsv` (`add_output`, lines 237-262); group runs (`--clade auto`) are read from
  `reports/groups/*.run.log` (lines 390-404). Reproduction commands are emitted per failed program and per record.
- **Do `--debug` flags raise the log level everywhere?** Yes: the launcher strips `--debug` and sets `MITOM_DEBUG=1`
  (`mitom`, `main()`), `mitom_pipeline.main` re-sets it for children (`mitom_pipeline.py:1565-1569`),
  `setup_root_logger(..., debug=True)` sets root to DEBUG (`mitom_log.py:47-67`), the `mitom.tools` logger is always
  DEBUG with propagation so `[CMD ...]` lines reach `run.log` even in normal runs, and `mito_common` adds a stderr
  handler for `mitom.tools` in debug mode (`mito_common.py:421-428`). The step scripts (`mito_core`, `mito_curate`,
  `mito_intron`, ...) inherit the env. Exceptions: findings 4.2 and 4.6 below.
- **Are child-process stderr captured and surfaced?** Yes for everything run through `run_tool` — command line,
  exit code, seconds and a stderr tail (20 lines, 200 with `--debug`) are logged as one `[CMD ...]` block that
  `mitom debug` parses (`mito_common.py:440-461`, `mito_debug.py:103-118`), and each sample's full engine stderr is
  kept in `<workdir>/logs/core.stderr.log` with the tail repeated in `sample.log` and in the `[RECORD]` line
  (`mitom_stage.py:178-192`). Exceptions: 4.2, 4.4, 4.7 below.
- **Do pipeline steps log the record id on failure?** Yes for annotation (round 1 and round 2), reorientation, final
  write, rescue, NCBI ERRORs (`record_event` in `mitom_pipeline.py:627, 668, 762, 846` etc.), and per-record curation
  errors (`curation_report.tsv` "error" rows, `mito_curate.py:1367-1369`). Gaps: 4.1 and 4.3 below.
- **Are failed records written to `reports/failed_logs`?** Yes for problem samples — `cleanup_intermediates` copies
  `sample.log` + `logs/*.log` of every `all_problem_samples()` sample into `reports/failed_logs/` before deleting
  intermediates (`mitom_output.py:322-334`). Gap: 4.1 below.

Problems found:

### 4.1 MEDIUM — round-2 failed records: their logs are deleted, and `mitom debug`'s "sample log" pointer dangles

- Files / lines:
  - `mitom_pipeline.py:664-670` — round-2 failures only `record_event(...)` (the event includes the workdir sample-log
    path written by `annotate_one`), they are added to **no** summary list;
  - `mitom_output.py:38-44` — `RunSummary.all_problem_samples()` does not include any round-2 category;
  - `mitom_pipeline.py:832-833` — `cleanup_intermediates(out_root, samples, summary.all_problem_samples())` then
    deletes `03.anno_round2_reoriented/<id>/logs/core.stderr.log` **without** copying it to `reports/failed_logs`.
- Why wrong: the run-level story is "failed records' logs are always kept (reports/failed_logs/)"
  (`--keep_intermediate` help, `mitom_pipeline.py:1236`). For a round-2 failure the `[RECORD]` line survives, but
  `mitom_debug.record_logs` (`mito_debug.py:213-223`) finds neither a `failed_logs` copy nor the (deleted) workdir,
  so the record's own error is lost; only the parent `[CMD rc=...]` block (last 20 stderr lines) remains. A user
  debugging "record X kept the shifted round-1 annotation" cannot see the core's actual stderr.
- Realistic trigger: `--if_reannotation_2 yes` (or `auto`) where round 2 fails for a few records (timeout/blast
  crash), default `--keep_intermediate no`.
- Potential fix: add a `reannotated_failed` list to `RunSummary`, append round-2 failures to it, and include it in
  `all_problem_samples()` (it is what `cleanup_intermediates` uses to decide the `failed_logs` copies).

### 4.2 MEDIUM — several production subprocesses bypass `run_tool`, so they leave no `[CMD ...]` trail for `mitom debug`

- Files / lines:
  - `mito_classify.py:128-131` — `_blastx` uses raw `subprocess.run(...)`
  - `mito_nuc.py:146` — `run(cmd)` uses `subprocess.run(cmd, capture_output=True, text=True, check=True)`
  - `mito_nuc.py:1277-1282` — `_exonerate_all` runs exonerate with raw `subprocess.run` **and ignores the exit code**
    (`out = subprocess.run([...]).stdout`; a crashing exonerate silently yields no histone models);
  - `mitom_snakemake.py` worker (~line 88) — `subprocess.run(cmd, ...)` for the engine, no logging.
- Why wrong: `run_tool` is documented as "the one way the MitoM modules run programs" and is the source of the
  `[CMD rc=...]` blocks that `mitom debug` lists under "Failed programs". Programs started through these four raw
  call sites never appear there: a missing/crashing blastx during `--clade auto` classification shows up only as a
  record note "blastx failed: ...", with no command line, exit code or stderr tail; a crashing exonerate in
  `mitom nuc` is fully silent. The Snakemake worker additionally hides the engine's exit code from the log trail
  (it is only in the per-sample result json / `core.stderr.log`).
- Realistic trigger: `mitom run -i contigs.fasta` (`--clade auto` needs classify's blastx) with a broken BLAST
  install; `mitom nuc annotate` with a crashing exonerate.
- Potential fix: use `run_tool(...)` (it already handles timeouts and logging) at these sites; at minimum capture
  and log failures in `_exonerate_all` (`if p.returncode: logging.warning(...)`) and in the Snakemake worker.

### 4.3 LOW — table2asn non-timeout failures are treated as "no messages" / "no NCBI ERROR"

- Files / lines:
  - `ncbi_validate.py:119-126` — `validate_record` catches only `subprocess.TimeoutExpired`; the return code of
    `run_tool([table2asn, ...])` is never checked, and a missing `.val` simply yields zero rows;
  - `mito_tbl.py:558-566` — same pattern in `validate()` (only `TimeoutExpired` handled).
- Why wrong: when table2asn crashes (segfault, bad `.fsa`, licensing banner on a non-zero exit), the record is
  reported as "0 ERROR / 0 WARNING" and the summary line says the record is clean; in the pipeline, no `.sqn` is
  copied but `submission_summary.tsv` shows it without NCBI ERROR (`mitom_steps.py:150-167`). The run-level
  diagnosis "table2asn not found: skipped" also cannot be distinguished from a per-record crash.
- Realistic trigger: a record whose `.fsa`/`.tbl` triggers a table2asn bug (the code itself documents that
  table2asn 1.28 "hangs for ever" on some inputs and needs `-N` renaming — fragile enough that hard crashes occur).
- Potential fix: check `p.returncode`; on non-zero, emit a synthetic row
  `{"code": "MITOM.Table2asnFailed", "severity": "ERROR", ...}` (mirroring the existing
  `MITOM.Table2asnTimeout` row) and skip the record in the ERROR counts rather than counting it as clean.

### 4.4 LOW — Snakemake backend points `mitom debug` at a `sample.log` that is never written

- File / line: `mitom_snakemake.py:231-233`
  ```python
                record_event(obj.get("sample_name", ...), f"annotation failed ({phase_name}, exit "
                             f"{obj.get('returncode', '?')}): {last}  (sample log: {obj.get('workdir')}/sample.log)")
  ```
- Why wrong: the Snakemake worker (`render_snakemake_worker_script`) writes only `logs/core.stdout.log`,
  `logs/core.stderr.log`, `core.returncode.txt` and `core.runtime_sec.txt` — there is no `sample.log` (that file
  exists only in the threads backend, `mitom_stage.annotate_one`). The `[RECORD]` line therefore names a
  non-existent file; `mitom_debug.record_logs` still finds `core.stderr.log` next to it, so the damage is a wrong
  pointer/message rather than lost data — but anyone following the message by hand hits a missing file.
- Realistic trigger: `--scheduler_backend snakemake`, any failed sample.
- Potential fix: point the message at `logs/core.stderr.log`, or have the worker also write a `sample.log`.

### 4.5 LOW — `mito_intron.run_tool` reports only the first failure of each program; later failures are silent

- File / line: `mito_intron.py:203-208`
  ```python
          p = run_program(cmd, timeout=TIMEOUT)
          if p.returncode != 0 and what not in _FAILED:     # (2026-10-31 review: failures were silent: no introns)
              _FAILED.add(what)
              print(f"[mito_intron] WARNING: {what} failed (exit {p.returncode}): ...")
  ```
- Why wrong: `what` is the program name ("exonerate", "blastp", "cmsearch"). After the first failure of exonerate is
  reported, every later exonerate failure — possibly with a *different* stderr, on specific records — is swallowed
  completely (the function returns `""` and the record silently gets no introns). The one-per-program de-dup hides
  per-record diagnostic information that `mitom debug` would otherwise show via the intron step's captured log.
- Realistic trigger: exonerate works for most records but crashes on a few (e.g. a record with a huge repeat) —
  only the first crash is ever mentioned.
- Potential fix: warn per (record, program) — pass the record id into `what` — or print every failure but keep a
  per-program counter for the summary line.

### 4.6 LOW — in non-debug runs the step scripts' internal `[CMD ...]` lines are lost (no logging handlers)

- Files / lines: `mito_intron.py` / `mito_trna.py` / `mito_curate.py` / `mito_intron_ref.py` `main()` never call
  `logging.basicConfig`/`setup_root_logger`; `mito_common` attaches a `mitom.tools` stderr handler only when
  `MITOM_DEBUG` is set (`mito_common.py:421-428`).
- Why wrong: in the pipeline the parent's `run_tool` logs `[CMD rc=1 ...] python mito_intron.py ...` with the last
  20 lines of the script's stderr, which is usually enough — but the exact command lines of the programs run
  *inside* the step (exonerate args, blastp args) are logged at DEBUG on `mitom.tools` and, without any handler,
  vanish (`logging.lastResort` only prints WARNING+). With `--debug` they land in `intron_fix.log` and the step
  logs. The practical effect: for a non-debug run, `mitom debug`'s "Reproduce → the failed program by hand" section
  can reproduce the failed *step script*, but not the specific inner exonerate invocation that failed.
- Realistic trigger: non-debug `mitom run` where the intron step fails for one record; user asked to reproduce the
  failing exonerate call.
- Potential fix: give the standalone `main()`s a minimal `logging.basicConfig(level=..., stream=sys.stderr)` (or
  always attach the `mitom.tools` stderr handler, downgrading its level when not in debug mode).

### 4.7 LOW — `mito_debug.rerun()` quotes reproduce commands for POSIX shells only

- File / line: `mito_debug.py:299` — `return " ".join(shlex.quote(w) for w in res)`
- Why wrong: on Windows the printed commands (`mitom run -i 'P:\x\my file.gb' ...`) are single-quoted, which neither
  `cmd.exe` nor PowerShell accepts; a Windows user copying the "Reproduce" lines gets argument-splitting errors.
- Realistic trigger: `mitom debug mitom_out` on Windows with a space in a path.
- Potential fix: use `subprocess.list2cmdline` on Windows (`os.name == "nt"`) and `shlex.quote` elsewhere.

---

## Fix checklist

High priority (Windows blockers; the platform this installation runs on):

1. Replace `mp.get_context("fork")` with a portable context + explicit worker state
   (`mitom_steps.py:388`, `mito_intron.py:1613`, `mito_intron_ref.py:446`, `mito_intron_ref.py:2015`). (1.1)
2. Replace `-o /dev/null` for nhmmer with a temp file / `os.devnull` (`mito_rrna_hmm.py:242`). (1.2)
3. Make the stale-lock probe non-lethal on Windows (`mitom_pipeline.py:1478`). (1.3)

Medium priority:

4. Add `encoding="utf-8"` (and `newline="\n"`) to all GenBank/FASTA writes, at least `write_gb`
   (`mito_common.py:126`) and `read_fasta` (`mito_common.py:141`). (1.4)
5. Guard `tarfile.extractall(filter=...)` for older Pythons (`trna_ml.py:109`). (1.5)
6. Add stage-level timeouts in `run_cmd_capture` and the Snakemake worker (`mitom_stage.py:41`,
   `mitom_snakemake.py` worker). (1.6, 4.2)
7. Log the exception in `workflow_transfer_metadata`'s bare `except` (`mitom_pipeline.py:947`). (1.7)
8. Fix O2 `origin_cut` part selection for minus-strand origin joins — keep the part ending at `L`, drop the part
   starting at 0 (`mito_curate.py:655-663`). (2.1)
9. Emit no value column for value-less qualifiers in `.tbl` (`ncbi_validate.py:62`). (3.1)
10. Use the subphylum code for Tunicata/Cephalochordata under `--clade Chordata`, or warn loudly
    (`mitom_taxa.py:14`, `mitom_pipeline.py:307-309`). (3.2)
11. Track round-2 failures in `RunSummary` so their logs land in `reports/failed_logs`
    (`mitom_pipeline.py:664-670`, `mitom_output.py:38-44`). (4.1)
12. Route `mito_classify._blastx`, `mito_nuc.run`/`_exonerate_all`, and the Snakemake worker through `run_tool`
    (and check exonerate's exit code). (4.2)

Lower priority:

13. Close the sample logger on worker exception (`mitom_stage.py:147-199`). (1.8)
14. Use `write_gb` for the remaining GenBank writers (`mito_curate.py:1376`, `mito_trna.py:862`,
    `mito_ncr.py:133`, `mito_engines.py:361`). (1.9)
15. Verify/replace `os.execv` in the `mitom` launcher on Windows. (1.10)
16. Group predictions against all members in `mito_trna.group()`; handle multi-part input tRNAs in
    `from_input()`. (2.2, 2.3)
17. ASCII-clean the organism in `ncbi_validate.fasta_defline`. (3.3)
18. Surface table2asn non-zero exits as validation rows instead of "clean" (`ncbi_validate.py:119-126`,
    `mito_tbl.py:558-566`). (4.3)
19. Fix the Snakemake "sample log:" pointer (`mitom_snakemake.py:232`). (4.4)
20. Report per-record program failures in `mito_intron.run_tool` (`mito_intron.py:203-208`). (4.5)
21. Add a basic logging config to the standalone step `main()`s (`mito_intron.py`, `mito_trna.py`, ...). (4.6)
22. Windows-compatible quoting in `mito_debug.rerun()` (`mito_debug.py:299`). (4.7)
