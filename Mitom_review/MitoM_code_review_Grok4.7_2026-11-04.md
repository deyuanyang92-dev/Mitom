# MitoM code review

Reviewed by **Grok 4.7**. Date **2026-11-04**. Read-only. No source code, tests, or config were changed. The only file written is this review.

Line numbers are from the tree as read on that date. Each item is a problem plus a potential fix. Confirm the cited lines still match before editing. Do not treat an omitted behavior as approval of it.

Severity:

- **high** — a realistic input produces a wrong gene, wrong coordinates, a silently wrong translation, or drops a record that should be kept
- **medium** — wrong on a real but narrower input, or a failure the user cannot see from the exit code
- **low** — off-by-one at a boundary, a wrong counter, or a rare qualifier

NCBI codon differences below are from the genetic-code tables at <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>. Re-check that page before changing a translation table.

`mito_classify.code_of_lineage` (`mito_classify.py` 104–116) can already vote a code from reference genomes that carry `/transl_table`, including 13, 14, 21, 24, and 33. Do not replace that vote. The bugs below are the **hardcoded** tables and the paths that never consult the record's own `/transl_table`.

---

## 1. Genetic code and NCBI conventions

### 1.1 high — Trematodes are translated with code 9; NCBI code 21 has ATA = Met

**Where:** `mito_taxa.py` 30–31 and 41; `mito_tbl.py` 219–232; `mitom_pipeline.py` 307–310.

```python
PHYLUM_CODE = {"Chordata": 2, "Echinodermata": 9, "Hemichordata": 9, "Platyhelminthes": 9, ...}
**dict.fromkeys(("Turbellaria", "Trematoda", "Cestoda", "Monogenea"), "Platyhelminthes"),
```

```python
code = next((SUBPHYLUM_CODE[t] for t in tax if t in SUBPHYLUM_CODE), PHYLUM_CODE.get(phylum, 5))
```

`Trematoda` is not a phylum token, so it is folded to `Platyhelminthes` and then to code 9. `SUBPHYLUM_CODE` has no trematode entry, so a full lineage `Platyhelminthes; Trematoda; ...` still gets 9. The CLI text at `mitom_pipeline.py` 309 tells the user to pass `9` for Platyhelminthes. Exonerate can use code 21 (`mito_core.py` 108 lists 21 in `EXONERATE_CODES`); the code is never selected.

**Why:** NCBI table 21 (Trematode Mitochondrial Code) matches table 9 except **ATA = Met in 21 and Ile in 9**. There is no new internal stop, so table2asn accepts the protein and the methionine/isoleucine difference is silent.

**Input:** A trematode GenBank record with taxonomy containing `Trematoda` and no CDS `/transl_table`, run through `mitom gb2tbl` (or curated/annotated after the user follows the "9 Platyhelminthes" hint). Any ATA in a CDS is written as isoleucine.

**Fix:** If the lineage contains `Trematoda`, use 21. Do this before the synonym fold, in `PHYLUM_CODE` / `record_code` and in the CLI hint. Leave an existing `/transl_table=21` untouched. Keep code 9 only for flatworm lineages that are not trematodes, and see 1.6 for code 14.

### 1.2 high — `--clade Chordata` and the phylum table use vertebrate code 2

**Where:** `mito_taxa.py` 14 and 30–35; `mitom_pipeline.py` 307–311; `mito_tbl.py` 141 and 227–230.

```python
CLADE_CODES = {"Chordata": 2, "Echinodermata": 9, "Porifera-sponges": 4}
args.genetic_code = CLADE_CODES.get(args.clade, 5)
```

```python
SUBPHYLUM_CODE = {"Tunicata": 13, "Cephalochordata": 5}
**dict.fromkeys(("Vertebrata", "Craniata", "Tunicata", "Cephalochordata"), "Chordata"),
```

**Why:** NCBI table 2 is Vertebrata, not Chordata.

- Tunicates are table 13: AGA and AGG are glycine, not stops. Table 2 cuts the CDS at the first AGR.
- Lancelets are table 5 in this codebase (`SUBPHYLUM_CODE`): AGA and AGG are serine. Table 2 again cuts at AGR.

`mito_tbl.record_code` does check the literal tokens `Tunicata` and `Cephalochordata` before the phylum code. That path is fine and should stay. It is skipped when:

- the user passes `--clade Chordata` and omits `--genetic_code` (`mitom_pipeline.py` 311 sets 2 for the whole run)
- the lineage says only `Chordata`
- the lineage uses the old name `Urochordata` (not in `SUBPHYLUM_CODE` or `PHYLUM_SYNONYMS`). With `Chordata; Urochordata`, the phylum match is Chordata and the subphylum match fails, so the code is 2. With `Urochordata` alone, the phylum lookup fails and the fallback is 5, which is also not 13.

**Input:** A tunicate FASTA or GenBank run as `mitom run --clade Chordata` without `--genetic_code`. Or a GenBank record whose taxonomy says `Urochordata` and whose CDS have no `/transl_table`, passed to `mitom gb2tbl`.

**Fix:** Map tunicates / `Urochordata` to 13 and lancelets to 5 before applying Chordata → 2. Do not let `--clade Chordata` override a record whose own taxonomy or `/transl_table` says 13 or 5. Vertebrata / Craniata stay on 2.

### 1.3 high — curation falls back to code 5 and never writes `/transl_table`

**Where:** `mito_curate.py` 1288–1292 and 1364; `translation` at 932–949. A search of `mito_curate.py` shows `/transl_table` is only read, never assigned.

```python
def record_code(rec) -> Optional[int]:
    """The genetic code of a record: the most common CDS /transl_table (None without one)."""
...
code = a.genetic_code or record_code(rec) or 5
```

`translation()` calls `Seq(nt).translate(table=self.code)` and stores `/translation` only. The help text at 1301–1302 says the old bug was "code 5 was used for every record". The remaining fallback is still 5 whenever the CDS has no `/transl_table` and the user did not pass `--genetic_code`. This function does not look at taxonomy, so a vertebrate record is not rescued by `PHYLUM_CODE`.

**Why:** Table 5 reads AGA/AGG as serine and ATA as methionine. Table 2 reads AGA/AGG as stops. A vertebrate CDS with no `/transl_table` is extended through every AGR until a TAA/TAG. The boundaries are then stored, but the qualifier that table2asn will trust is still missing. `mito_tbl` later guesses from taxonomy (1.1, 1.2) and drops `/translation` (`ncbi_validate.py` 33, `SKIP_QUALS`). The submitted protein is retranslated with a different code from the one that set the ends.

The same split happens when `--genetic_code` disagrees with an existing `/transl_table`: ends follow the new code, the qualifier stays on the old one.

**Input:** A vertebrate mitogenome whose CDS lines have no `/transl_table`, run through `mitom curate` without `--genetic_code`.

**Fix:** Choose the code in one place (explicit flag, else existing `/transl_table`, else the taxonomy rules in 1.1 and 1.2, else stop and ask). Write that integer onto every CDS as `/transl_table` before computing `/translation`. Do not default a missing code to 5.

### 1.4 medium — the curation validator's FASTA header omits the mitochondrial code

**Where:** `ncbi_validate.py` 65–79 and 113.

```python
if gencode:
    mods.append(f"[{'gcode' if nuclear else 'mgcode'}={gencode}]")
...
def record_gencode(rec) -> Optional[int]:
    codes = [int(f.qualifiers["transl_table"][0]) for f in rec.features
             if f.type == "CDS" and f.qualifiers.get("transl_table", [""])[0].isdigit()]
    return max(set(codes), key=codes.count) if codes else None
```

```python
(d / "rec.fsa").write_text(fasta_defline(rec, record_gencode(rec), nuclear) + "\n" + str(rec.seq) + "\n")
```

**Why:** `validate_record` asks only for `/transl_table`. After 1.3, a record curated with an assumed code still has no qualifier, so the `.fsa` has `[location=mitochondrion]` and no `[mgcode=]`. `mito_tbl.defline` (264–279) does call `record_code`, which can fill a code from taxonomy. The validator path does not. table2asn is then checking a different translation from the one curation just computed. Mitochondrial TGA is a stop in the standard code.

**Input:** The same vertebrate record as 1.3, when `mitom curate` runs table2asn (`mito_curate.py` 1360–1373 calls `nv.validate_record`).

**Fix:** Pass the code actually used for translation into `fasta_defline`. After 1.3 writes `/transl_table`, `record_gencode` will see it; until then, do not validate with an empty code.

### 1.5 high — `/pseudo` is submitted as the text `None`

**Where:** `mito_curate.py` 896–900; `mito_tbl.py` 191; `ncbi_validate.py` 55–61.

```python
x.qualifiers["pseudo"] = [None]    # value-less /pseudo
```

```python
vals = [str(v) for v in f.qualifiers[k]]
```

```python
out.append(f"\t\t\t{k}\t{' '.join(str(v).split())}")
```

**Why:** INSDC `/pseudo` is a flag with no value. Biopython will write a GenBank `/pseudo` from a `None` value. Both table writers then call `str(None)`, which is the four characters `None`. The feature table line is `pseudo	None`. table2asn does not treat that as the flag, so a CDS that was marked pseudo because of internal stops still fails `InternalStop` and can still be translated. The note at `mito_curate.py` 898 shows the author intended a valueless qualifier.

**Input:** Any CDS whose translation, after curation, still contains an internal `*`. Rule P10 sets `/pseudo` on the gene and the CDS. `mitom gb2tbl` or `validate_record` emits `pseudo	None`.

**Fix:** In both writers, emit a valueless qualifier with an empty value column (`pseudo` and no text). Never call `str` on `None`. Keep Biopython's `None` in the SeqFeature so the GenBank writer still emits `/pseudo`.

### 1.6 medium — the static table cannot emit code 14 (TAA = Tyr)

**Where:** `mito_taxa.py` 30 and 41; `mito_tbl.py` 230. No assignment in this repository sets 14 except whatever `code_of_lineage` copies from a reference `/transl_table`.

**Why:** NCBI table 14 (Alternative Flatworm Mitochondrial Code) matches table 9 except **TAA = Tyr, not a stop**. TAG stays a stop. Translating those records with table 9 cuts the protein at the first TAA. NCBI's own note on table 14 says the TAA=Tyr evidence is disputed, so do not move the whole phylum to 14. What the code does today is stronger than that dispute: `Cestoda` is always folded to Platyhelminthes → 9, and a record that needs 14 cannot get it from the fallback.

**Input:** A cestode or other flatworm GenBank file with no `/transl_table`, passed to `mitom gb2tbl`. An existing `/transl_table=14` is kept by `mito_curate.record_code` (it only reads the qualifier) and should stay kept.

**Fix:** Do not guess 9 over an existing 14. If the taxonomy fallback is used, prefer the code on the closest reference genomes (`code_of_lineage`) over `PHYLUM_CODE`. Add 14 as a possible result for the taxa whose references use it. Do not hard-code every flatworm as 14.

### 1.7 medium — every hemichordate is code 9; two families are 24 and 33

**Where:** `mito_taxa.py` 30, `"Hemichordata": 9`. `mito_core.py` 108–109 already knows 24 and 33 are real tables (they are not built into exonerate and are passed as an explicit 64-codon table). Nothing in `mito_taxa.py` or `mito_tbl.py` selects them.

**Why:** Confirm the codon lists on the NCBI page before editing translations. The assignment itself is not in doubt:

- Table 24, Rhabdopleuridae: not table 9. AAA and AGG are not the table-9 amino acids (table 9 has AAA = Asn and AGG = Ser).
- Table 33, Cephalodiscidae: TAA = Tyr, so table 9 truncates at the first TAA; AGG also differs.

Enteropneusts were not re-checked here. Do not change the whole phylum off 9.

**Input:** A *Rhabdopleura* or *Cephalodiscus* record with no `/transl_table`, annotated with the phylum default or submitted through `gb2tbl`.

**Fix:** If the lineage contains `Rhabdopleuridae`, use 24. If it contains `Cephalodiscidae`, use 33. Leave other hemichordates on 9 until a reference `/transl_table` says otherwise.

### 1.8 medium — trnL1/trnL2 and trnS1/trnS2 follow MITOS2, not NCBI (intentional)

**Where:** `mito_names.py` 126–135, used by curation at `mito_curate.py` 1199 (`trna_symbol`).

```python
def trna_symbol(aa: str, ac: str) -> str:
    """... NCBI's vertebrate records (human NC_012920.1) number Leu / Ser the other way;
    MitoM follows MITOS2 (user's decision 2026-10-31)."""
    if aa == "L":
        return "trnL1" if ac.endswith("ag") else "trnL2" if ac.endswith("aa") else "trnL"
    if aa == "S":
        return "trnS1" if ac.endswith("ct") else "trnS2" if ac.endswith("ga") else "trnS"
```

**Why:** This is a documented choice, not an accident. It still renames genes on the way to NCBI. Anticodon `tag` (CUN) becomes `trnL1` and `taa` (UUR) becomes `trnL2`. Anticodons ending in `ct` (AGN / AGY) become `trnS1` and `tga` (UCN) becomes `trnS2`. Human RefSeq NC_012920 is the other way: trnL1 = UUR (`taa`), trnL2 = CUN (`tag`), trnS1 = UCN (`tga`), trnS2 = AGY (`gct`). Curating a RefSeq record that already has anticodons swaps all four names.

**Input:** `mitom curate` on NC_012920.1, or any vertebrate record whose tRNAs carry `/anticodon`.

**Fix:** Write the NCBI numbers on features that will be submitted. Keep the MITOS2 numbers only inside the caller that needs them. Do not change both conventions in one symbol without a flag; the 2026-10-31 comment records a deliberate choice, so the submission path is the part to split out.

### 1.9 medium — a bare `atp` is renamed to ATP6 during curation

**Where:** `mito_names.py` 39 and 71–79; `mito_curate.py` 1174–1176.

```python
"ATP6": ("atp6", "atpase6", "atp", "atpsynthasef0subunit6"),
```

```python
if g not in pcg and g in PCG_SYN and g != "atp":   # a bare "atp" stays ambiguous
```

```python
name = PCG_SYN.get(key(cur)) or PCG_SYN.get(key(f.qualifiers.get("product", [""])[0]))
```

**Why:** `gene_name()` refuses the bare synonym `atp`. `Curator.names` does not call `gene_name()`. It uses `PCG_SYN`, which still maps `atp` → `ATP6`. A feature that only says `atp` (an atp8 or atp9 that lost its digit, or a generic ATPase note stored in `/gene`) becomes ATP6, with the ATP6 product string.

**Input:** A CDS whose `/gene` is `atp` or `ATP`, run through `mitom curate`.

**Fix:** Use the same exclusion as `gene_name()`: if the normalized token is `atp`, leave it and log it. Do not add `atp` to the ATP6 synonym tuple.

### 1.10 medium — `/transl_except` must contain `..`, so single-base and `join` exceptions are deleted

**Where:** `mito_curate.py` 194–210. The writer of the one-base form is `mito_common.py` 370–377.

```python
_EXCEPT_RE = re.compile(r"pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)")
```

```python
pos = f"{end - tail + 1}..{end}" if tail > 1 else f"{end}"
return f"(pos:{pos},aa:TERM)"
```

**Why:** NCBI writes a one-base exception as `(pos:15432,aa:TERM)` or `(pos:complement(5),aa:Trp)`, and a stop split by the origin as `pos:join(...)`. The regex requires a two-number range. `internal_excepts` drops every non-matching qualifier. TERM exceptions are supposed to be rewritten by the caller, but a non-TERM exception (an ambiguous codon translated as Trp, the case named in the comment on line 198, NC_001636.1 ATP8) is removed. The residue becomes `X` on the next translation. A one-base TERM written by `transl_except_term` itself does not match this regex on a later pass.

**Input:** A CDS with `/transl_except="(pos:15432,aa:Trp)"` or `(pos:complement(5),aa:Trp)`, or `(pos:join(16000..16000,1..1),aa:TERM)`.

**Fix:** Accept a single number, `complement(single number)`, and `join(...)`. Keep any exception whose positions still lie on the final CDS and whose amino acid is not the stale TERM that this pass is about to rewrite. Make `transl_except_term` and the parser agree on the one-base form.

### 1.11 medium — start/stop checks ignore `/codon_start`

**Where:** `mito_curate.py` 888–902 versus 937.

```python
nt = s[a - 1:b]
...
if nt[:3] not in self.starts:
```

```python
nt = str(f.extract(self.rec.seq)).upper()[int(f.qualifiers.get("codon_start", ["1"])[0]) - 1:]
```

**Why:** `translation()` skips `codon_start - 1` bases (the qualifier is 1-based: 1, 2, or 3). `pcg_check` always reads the first three bases of the feature and the tail of that same untrimmed string. With `/codon_start=2` or `3`, the start-codon test (P5), the internal-stop test (P4/P10), and the T/TA tail test (P6) run in a different frame from the translation. A legal partial CDS can be marked 5′-partial, 3′-partial, or `/pseudo`.

Nothing in this repository writes `codon_start` 0. The bug is 2 and 3 on incoming records. `tests/test_mito_nuc.py` 122–133 is a separate nuclear ORF path and is not this function.

**Input:** An imported CDS of 100 nt with `/codon_start=2`, a start codon in frame at base 2, and a stop in that same frame. P5 looks at bases 1–3, which are not that start codon.

**Fix:** Apply the same `codon_start` trim before every codon test in `pcg_check`. After a 5′ edit, recompute `codon_start` as 1, 2, or 3 only.

### 1.12 medium — any organelle string containing `plastid` is written as a chloroplast

**Where:** `mito_tbl.py` 244–252.

```python
def location(rec) -> Tuple[str, str]:
    org = " ".join(source_of(rec).get("organelle", [])).lower()
    if "mitochondr" in org:
        return "mitochondrion", "mgcode"
    if "chloroplast" in org or "plastid" in org:
        return "chloroplast", "gcode"
```

**Why:** `/organelle="plastid:apicoplast"` contains the substring `plastid`. The FASTA modifier becomes `[location=chloroplast]`. An apicoplast is not a chloroplast. The same branch forces `gcode` rather than leaving the original location string.

**Input:** A GenBank record with `/organelle="plastid:apicoplast"` passed to `mitom gb2tbl`.

**Fix:** Return `chloroplast` only for `plastid:chloroplast` or a string that is exactly `chloroplast`. Pass any other `plastid:...` value through as the location modifier.

### 1.13 medium — the download query requires the adjacent words "complete genome"

**Where:** `mito_fetch.py` 50.

```python
MT_CORE = 'mitochondrion[Filter] AND "complete genome"[Title] NOT partial[Title]'
```

**Why:** Entrez double quotes are a phrase search. This matches a title containing the two words `complete genome` next to each other (`mitochondrion, complete genome`). It does not match the common INSDC title `complete mitochondrial genome`. The same repository's fetch helpers in `benchmark/mito_sample.py` and `scripts/fetch_rrna_refs.py` do not quote this phrase. The miss is silent: those records are simply absent from `mitom fetch`.

**Input:** `mitom fetch` for a taxon whose mitogenome titles use `complete mitochondrial genome`.

**Fix:** Drop the quotes, or combine the real title patterns (`complete genome[Title]`, `complete mitochondrial genome[Title]`, `complete mtDNA[Title]`) and still exclude `partial[Title]`.

### 1.14 low — *Radopholus* is code 14, and the nematode default is 5

**Where:** `mito_tbl.py` 230, `PHYLUM_CODE.get(phylum, 5)`. Nematoda is not in `PHYLUM_CODE`, so the default is 5.

**Why:** NCBI table 14 lists *Radopholus similis* and *R. arabocoffeae* with TAA = Tyr. Table 5 treats TAA as a stop and truncates those genes. *C. elegans* and ascarids really are table 5, so the phylum default must stay 5.

**Input:** A *Radopholus* record with no `/transl_table`.

**Fix:** If the organism or lineage contains `Radopholus`, do not assume 5. Use 14, or require an explicit `/transl_table`.

### 1.15 low — `subunit 1` is not a COX synonym in the table curation uses

**Where:** `mito_names.py` 36; `mito_curate.py` 1175.

```python
"COX1": ("cox1", "coi", "co1", "coxi", "cytochromecoxidasesubuniti"),
```

**Why:** The synonym is the Roman numeral `subuniti`. `Curator.names` looks up `PCG_SYN` and does not call `gene_name()`, whose product regex does accept an Arabic numeral. A CDS with an empty `/gene` and `/product="cytochrome c oxidase subunit 1"` is not renamed to COX1. `/gene=cox1` or `COI` still works.

**Input:** A CDS with only that product string.

**Fix:** Point N1 at `gene_name()`, or add the Arabic `subunit1` / `subunit2` / `subunit3` keys next to the Roman ones.

---

## 2. Annotation algorithms

### 2.1 high — two stop codons in the aligned frame replace the reading frame

**Where:** `mito_core.py` 233–241, applied at 385–386, then consumed by `place_cds` (244–291).

```python
def coding_frame(window: str, tab: int, tae: int, table: CodonTable.CodonTable) -> int:
    stops = set(table.stop_codons)
    n = [sum(1 for c in range(tab + d, tae - 2, 3) if window[c:c + 3] in stops) for d in range(3)]
    if n[0] < 2:
        return tab
    return tab + min(range(3), key=lambda d: n[d])
```

**Why:** The aligned frame is kept only when it has 0 or 1 stop. At 2 stops the function moves the CDS start by 1 or 2 nt to whichever frame has the fewest stops. `place_cds` then cuts every codon on that new frame, so both coordinates and the translation change. The comment on lines 239–240 says this guard exists because a short gene's other frames may contain no stops. The guard does not cover the short gene once the true frame itself has two stops.

Animal mitochondrial code has two stop codons (TAA, TAG) out of 64. An atp8 of about 50 codons therefore expects roughly 1.6 stops in a random frame. A correct atp8 frame that contains two stops (two sequencing errors, or two real premature stops that should be reported) loses to a frame with 0 or 1. The shift of the whole gene is then written as the CDS.

The intended case is real and must keep working: the comment cites NC_086642.1 nad5, 78 internal stops, both ends 5 nt off. `min(range(3), key=...)` prefers the lowest index on a tie, so an equal count does not leave frame 0.

**Input:** Build a plus-strand window whose exonerate protein alignment is in-frame and covers atp8, and put two in-frame TAA codons inside `[tab, tae)` while another frame has fewer than two. `coding_frame` returns `tab+1` or `tab+2`.

**Fix:** Shift only when the new frame has strictly fewer stops and the translated new frame still matches the aligned reference protein over the aligned columns. Do not shift a CDS shorter than a stated codon count (atp8-sized). Leave the 78-stop nad5 case on the shifted frame.

### 2.2 high — intron-end extension treats a compound location as one interval

**Where:** `mito_intron.py` 992–1023 (`prev_feature_end`, `next_feature_start`), the same pair at 914 (`core_ends`), and `overlap` at 1050–1051.

```python
fs, fe = int(ft.location.start) + 1, int(ft.location.end)
if overlap((fs, fe), g_span) >= 30:
    continue
end_tx = fe if m.strand == 1 else L + 1 - fs
```

```python
def overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
```

**Why:** One `start` and one `end` are not the exons. `overlap` is a linear interval: if the first number is greater than the second, the overlap is 0 and the neighbor is invisible.

- Plus strand, origin-spanning `join`, parts in transcript order. The first part starts near the end of the record and the last part ends near base 1. `fs > fe`, overlap is 0, and `extend_to_stop` is not stopped by that neighbor. Its `limit` stays at `L` (line 1013).
- Plus strand, ordinary multi-exon gene. `start` of the first exon and `end` of the last exon include the intron. A different gene sitting in that intron (cox1 inside a coral nad5 intron) overlaps the outer span by ≥ 30 and is skipped by the `continue`. The outer exons do not block the extension either, because they were not tested as separate parts.
- Minus strand, several exons in transcript order (5′ exon at the higher coordinate). `start` of the first part is greater than `end` of the last part, overlap is 0, and the same extension runs to the end of the sequence.

`core_ends` uses the same pair, so an existing origin-spanning CDS also fails the overlap test and its curated end is not reused.

**Input:** A circular record whose nad5 is `join(15000..16000,1..400)` and whose next gene is a simple CDS. Extend nad5. Or a plus-strand nad5 with an intron that contains cox1. The cox1 feature is skipped and the 3′ walk is not limited by it.

**Fix:** Walk `location.parts`. Compute overlap and transcript coordinates per exon. A neighbor that shares ≥ 30 nt with the gene being rebuilt can be ignored; a neighbor that only falls inside an intron of that gene must still limit the extension. For a circular record, test the wrapped interval instead of `min/max` on a reversed pair.

### 2.3 high — the reference-intron 3′ end keeps an internal stop or swallows the next codon

**Where:** `mito_intron_ref.py` 1158–1187. `gene_model` uses the returned coordinate as the last exon end at 1412–1416. There is no next-gene argument.

```python
p = x
while p < e0 + 300:
    ...
    if G.sub(strand, p, p + 3) in stops:
        if p >= e0 - 90:
            return p + 3, False
        break
    p += 3
q = x + 3 * ((e0 - x) // 3)
tail = G.sub(strand, q, q + 2)
if tail == "TA":
    return q + 2, False
if tail[:1] == "T":
    return q + 1, False
return (p + 3, False) if p < e0 + 300 else (q, True)
```

**Why:** Three different cuts:

1. Any in-frame stop from 90 nt before the reference end to 300 nt after it is accepted (`p >= e0 - 90` and the loop bound `e0 + 300`). The CDS can run into the next gene. Nothing in `gene_end` consults the next feature.
2. A stop more than 90 nt before the reference end executes `break`, leaving `p` on that stop. Control then falls through to the T/TA test at the reference end `q`. If that codon starts with T (about a quarter of codons) the function returns an incomplete stop at `q` and the earlier stop stays inside the CDS. The last `return` is the line that would have ended at the early stop, and it runs only when the reference-end codon does not start with T. The comment on that line ("an early stop ends the CDS there") does not match this order.
3. `TA` is tested before the codon is known to be a stop. `TAT`/`TAC` (tyrosine) are cut to `TA`. `TTG` and any other `Tnn` sense codon are cut to `T`. The real incomplete stop `T` followed by the next gene's `ATG` is the string `TATG`, which is taken as `TA` and consumes the `A`.

**Input:** Reference exon end `e0` on a frame where an in-frame TAA sits 120 nt earlier and the codon at `e0` is `TTG`. The CDS keeps the TAA and ends on the `T` of `TTG`. Separately, a reference end whose next three bases are `TATG` (incomplete T, then ATG) returns `q+2`.

**Fix:** On an in-frame stop before the reference end, end at that stop (`p+3`) and do not look at `q`. Past the reference end, accept a complete stop only within a short window (about 30 nt) and stop before the next gene. Emit T or TA only when fewer than three bases remain, or the next feature starts before a full codon, and the bases that are taken are not a complete sense codon. If the first base is `T` and a full codon is present and is not a stop, do not shorten it. Check `TAA`/`TAG` before `TA`.

### 2.4 high — reverse-strand nuclear rDNA is stranded from the subject, not the query

**Where:** `mito_nuc.py` 487–488, 524–537, 563, 624–628, 785, 815–818. The regression that states the intended result is `tests/test_mito_nuc.py` 136–149 (`test_rdna_complete_unit_reverse_strand`, decorated `@needs_tools`). It reverse-complements a complete rDNA unit and expects every gene on strand `"-"`.

The blastn query is the contig and the subject is the reference (`annotate_rdna`, outfmt `qstart qend sstart send`). The hit tuple is `(qstart, qend, sstart, send, bits, qseq, sseq)`.

```python
def _strand(h) -> int:
    return 1 if h[3] >= h[2] else -1
```

```python
qp, sp = qs - 1, ss - st
for x, y in zip(qa, sa):
    if x != "-":
        qp += 1
    if y != "-":
        sp += st
        if sp == b:
            return qp if x != "-" else qp + 0.5
```

```python
lo, hi, st = c[0][0], c[-1][1], _strand(c[0])
...
if hi_p < lo_p or (hi_p - lo_p < 20 and not (p_lo or p_hi)):
    continue
```

**Why:** For blastn, a hit on the minus strand of the query has **qstart > qend** and **sstart < send** (the reference stays plus). `_strand` looks only at subject coordinates, so `send >= sstart` and the hit is labelled plus. `map_pos` then always walks the query toward larger coordinates. The alignment string walks from `qstart` toward `qend`, which is downward. A single HSP therefore has `lo > hi` and is dropped at line 628, or, if several HSPs are chained with the plus-strand test (both query coordinates must increase, lines 503–505), the 18S and 28S ends are projected on the wrong side of the locus.

This is not the histone path. Histone tblastn uses the contig as the subject, so `sstart > send` really is the minus strand, and `test_histone_minus_strand_translation` is a different coordinate system. Do not change that path to fix this one.

**Input:** The contig built by `test_rdna_complete_unit_reverse_strand`: one rDNA unit, reverse-complemented. Capture one blastn line. It should show `qstart > qend` and `sstart < send`. `_strand` of that tuple is `+1`.

**Fix:** Query strand is −1 when `qstart > qend`. In `map_pos` and `project`, step the query coordinate from `qstart` toward `qend` (decrement when `qstart > qend`). Chain HSPs in that query direction. Keep subject order for "does this HSP cover this reference base", not for the query strand. After the change, `test_rdna_complete_unit_reverse_strand` should report strand `"-"` and the reverse-complement coordinates already computed in the test (`n - end + 1`, `n - start`).

### 2.5 low — same-strand tRNA overlap uses a length that is short by one

**Where:** `trna_ml.py` 332–333, 341, and 432–435.

```python
def _ov(a, b) -> int:
    return min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
```

```python
if any(x["strand"] == c["strand"] and _ov(x, c) > 0.5 * min(x["end"] - x["start"], c["end"] - c["start"])
       for x in final):
    continue
```

**Why:** `_ov` is a closed-interval length (`+ 1`). `_loci` at line 341 compares it with `end - start + 1`. `call` compares it with `end - start`, which is one less than the true length. The threshold is half a base too small, so a call is dropped when the overlap is just over half of `(length - 1)` and not over half of `length`.

**Input:** Two same-strand calls, shorter length 70 (`end - start` = 69). True half is 35. The code uses 34.5. An overlap `_ov` of 35 is rejected (`35 > 34.5`) even though it is not more than half of 70.

**Fix:** Use `0.5 * min(x["end"] - x["start"] + 1, c["end"] - c["start"] + 1)`, the same length `_loci` already uses.

### 2.6 medium — the atp8 / nad4L fallback ORF cannot cross the origin

**Where:** `mito_core.py` 92–93, 491–510, called at 480–484. Circular search copies a prefix onto the end (`758`, `ORIGIN_EXT = 3000` at line 110), which extends genes that run off the end. It does not place the end of the genome in front of position 0.

```python
N = len(genome)
seq, bs = (genome, nb.start) if nb.strand == 1 else (rc, N - nb.end)
...
for e in range(max(3, bs - 25), min(N, bs + 11)):
    ...
    while i >= 0 and e - i <= longest and seq[i:i + 3] not in stops:
```

**Why:** `orf_before` only runs when blastx did not already call atp8 (before atp6; arthropods and echinoderms) or nad4L (before nad4; arthropods, annelids, chordates, nemerteans). It looks for a stop within 25 nt before the neighbor and 10 nt after it, then walks left, and it stops at index 0. If the neighbor starts near base 1, the upstream ORF is at the end of the circular chromosome. The appended prefix is a copy of the start, so that walk never sees it.

**Input:** A circular arthropod genome where atp6 starts within the first 30 nt, atp8 lies entirely in the last 200 nt, and blastx missed atp8. `orf_before` returns `None`.

**Fix:** Take `longest` nt of circular sequence immediately upstream of the neighbor (wrapping at 0), search that linear string, and map the ORF coordinates back onto the original record.

### 2.7 medium — `--stop_rule homolog` cuts a sense codon down to T or TA

**Where:** `mito_intron.py` 880–885. The default is `"first"` (`1573` and the function default on 853). This fires only when the user passes `--stop_rule homolog` and the alignment already reached the homolog's C-terminus (`qmiss_c <= 1`).

```python
if rule == "homolog" and qmiss_c <= 1:
    if s[b:b + 3] in stops and b + 3 <= limit:
        return (exons[:-1] + [(a, b + 3)]), False
    if s[b:b + 2] == "TA" and b + 2 <= limit:
        return (exons[:-1] + [(a, b + 2)]), False
    return (exons[:-1] + [(a, b + 1)] if s[b:b + 1] == "T" and b + 1 <= limit else exons), False
```

**Why:** Same order as 2.3. `TAT` becomes `TA`. `TTG` becomes `T`. `T` followed by `ATG` (`TATG`) becomes `TA`. A complete stop is handled first, which is right; the incomplete-stop branch then runs on every other codon that starts with T.

**Input:** `mitom intron --stop_rule homolog` on a CDS whose next in-frame codon after the homolog end is `TAT` or `TATG`.

**Fix:** If three bases are available and they are a stop, take all three. If they are a sense codon, do not shorten them. Take `TA` or `T` only when the feature boundary or the sequence end leaves fewer than three bases, and do not prefer `TA` when the third base would complete a sense codon that is still inside the gene.

### 2.8 medium — ARWEN / ARAGORN hits are stored with an empty anticodon

**Where:** `trna_ml.py` 349–357 and 385–388.

```python
out.append({..., "anticodon": "", "ac_known": 1, "ac_match": 1})
```

```python
acs = Counter(h["anticodon"] for h in v if h.get("anticodon"))
...
"anticodon": acs.most_common(1)[0][0] if acs else "",
```

**Why:** `ext_hits` copies the family (for Leu/Ser, from `EXT_LS` keyed by the real anticodon) and then stores the anticodon as `""`. It also forces `ac_known` and `ac_match` to 1, so the model is told the anticodon was checked. The vote on line 385 ignores empty strings. If no covariance-model hit is in that locus, the candidate's `anticodon` stays `""` and that is what is written out. The Leu/Ser family name can still be right while `/anticodon` is blank.

**Input:** A tRNA locus where the only calls are ARWEN and/or ARAGORN (no CM hit). The `Pred` objects have a three-letter anticodon; the candidate does not.

**Fix:** Keep `p.anticodon` on the pseudo-hit. Set `ac_known` / `ac_match` from that anticodon and the genetic code, using the same codon check the CM hits use (`trna_ml.py` around 320–325).

---

## 3. Circular coordinates, orientation, and the pipeline

### 3.1 high — a gene split across the origin is rotated to the 3′ piece and skips round 2

**Where:** `mitom_orient.py` 59–72, 226–253.

```python
if strand == -1:
    return max(same, key=lambda f: int(f.location.end)), gene
return min(same, key=lambda f: int(f.location.start)), gene
```

```python
moved = str(ro.seq) != str(a.seq)
...
if moved and not linear:
    edge = ...
cross = ... len(f.location.parts) > 1 and any(int(p.end) == n ...) and any(int(p.start) == 0 ...)
```

**Why:** The comment on lines 60–61 says the function wants the 5′-most piece. For a plus-strand gene it uses the smallest `start`. That is the 5′ exon only when both exons sit on one linear pass. When the 5′ exon is at the end of a circular record and the 3′ exon starts near base 1, `min(start)` selects the 3′ exon. If that exon starts at 0, the rotation is 0.

`reorient_check` then does not repair it:

- The cut-end test runs only when the sequence actually changed (`moved`).
- The origin test requires **one** feature with two parts, one touching 0 and one touching `n`. Two separate features do not match.
- `COX1` is in `STANDARD_GENES` (127–128), so a second COX1 feature means "cox1 is present" and does not request round 2.

`starts_with_gene` uses the same picker, so `cox1_first` can still be reported as yes. A single `join` location is a different path: `feature_start_pos` (`mito_common.py` 360–364) uses `parts[0]`, which is the 5′ part. This bug is the two-feature case, not the join case.

**Input:** A circular GenBank record with two plus-strand CDS, both named COX1: one at `15000..16000` (5′ exon) and one at `1..400` (3′ exon), no `join`. After `mitom reorient` the sequence is unchanged, the 5′ exon is still at the end, and `reannotation_decision.tsv` says `shifted_round1`.

**Fix:** If two features of the orientation gene are on one strand and one piece touches the start of the record while another touches the end, rotate to the 5′ piece (the piece at the end on the plus strand; the piece at the start on the minus strand). Even when the sequence does not move, record that pair as a split gene and send the record to round 2. Do not let `min(location.start)` define the 5′ end of a circular gene.

### 3.2 high — once any intron feature is added, every later CDS loses its origin-crossing intron note

**Where:** `mitom_orient.py` 326–365, inside `add_exon_gap_introns`.

```python
added = []
for cds in [f for f in rec.features if f.type == "CDS" and len(f.location.parts) > 1]:
    ...
    for k, (a, b) in enumerate(zip(parts, parts[1:]), 1):
        ...
        if not 0 < e - s < len(rec.seq) // 2 or s >= e:
            continue
        ...
        added.append(...)
    if added:
        cds.qualifiers["note"] = [x for x in cds.qualifiers.get("note", [])
                                  if not x.startswith(INTRON_NOTE)
                                  and "crosses the origin: given by the join only" not in x]
```

**Why:** `added` is one list for the whole record. The `if added` test is inside the CDS loop. After the first CDS appends an intron, every following CDS enters that branch. A later CDS whose only intron crosses the origin hits `continue` at line 350 (`s >= e` or the gap test), so no intron feature is created for it, and the note that was the only record of that intron is deleted. The intron is gone from the file.

**Input:** A record with two multi-exon CDS. The first has an internal intron that this function restores. The second has only an origin-crossing intron, recorded as a note (`intron N (...) crosses the origin: given by the join only`) and not as a feature. After the call, the second CDS has neither the feature nor the note.

**Fix:** Strip notes only from the CDS whose intron numbers were actually written as features in this call. Track that set per CDS, not the shared `added` list.

### 3.3 high — a RefSeq COMMENT drops the INSDC record without comparing the sequences

**Where:** `mitom_inputs.py` 385–398. The sequence-identity pass is later, at 400–408, and only sees records still in `by_acc`.

```python
for meta in by_acc.values():
    if is_refseq_acc(meta.acc_base):
        for src in refseq_source_accessions(meta.record):
            superseded.setdefault(src, meta.acc_versioned)
for key in [k for k in by_acc if k in superseded]:
    ...
    del by_acc[key]
```

**Why:** Every accession named in the RefSeq COMMENT is deleted. The sequences are not compared, and the versions are not compared. An updated GenBank record that no longer matches the `NC_` is removed, and the older `NC_` is what gets annotated. A second accession mentioned in the same comment is removed as well.

**Input:** An `NC_` record whose COMMENT says it was derived from `KX000001`, plus a newer `KX000001.2` whose sequence differs (a correction, or an added stretch). `KX000001.2` is deleted before the hash comparison runs.

**Fix:** Drop the INSDC record only when its sequence hash equals the `NC_` sequence, or when the COMMENT states the derivation and the INSDC version is older. If the sequences differ, or the INSDC version is newer, keep the INSDC record and say so in `dedup_report.tsv`.

### 3.4 medium — the process exits 0 unless every annotation failed

**Where:** `mitom_pipeline.py` 443–446 (`workflow_annotate`) and 854–857 (`workflow_run`).

```python
if total and len(summary.annotation_failed) == total:
    logging.warning(...)
    return 1
return 0
```

**Why:** One success out of 100 is exit 0. `final_write_failed` is written to `reports/final_write_failed.txt` (`784–785`) and is not part of this test, so a run that annotated records and then failed to write every final GenBank file still returns 0. `curated.ok`, `submitted.ok`, and `maps.ok` are display flags only. A caller that uses `mitom run && next_step` will treat a partial batch as success. The failure lists do exist; the exit code does not reflect them.

**Input:** A two-record input where the second record raises in the annotation core, or where final GenBank writing appends the sample to `final_write_failed`. The process status is 0.

**Fix:** Return 1 when `annotation_failed`, `final_write_failed`, or `reorient_failed` is non-empty. Keep writing the per-record lists either way.

### 3.5 medium — round 2 replaces round 1 whenever the round-2 file exists

**Where:** `mitom_pipeline.py` 694–698.

```python
if do_phase2 and not is_no_gene:
    r2 = results2.get(iid, {})
    if r2.get("best_gb") and Path(r2["best_gb"]).exists():
        chosen_gb = Path(r2["best_gb"])
        src_tag   = "phase2"
```

**Why:** A zero exit from the second annotation and a file on disk are enough. The gene set is not compared with the shifted round-1 file. If round 2 drops a standard gene that round 1 had, the worse file is the one that continues into curation, and the round-1 GenBank is not retained as the result.

**Input:** A record sent to round 2 whose second annotation writes a GenBank that is missing cox1 (or another `STANDARD_GENES` member) which the reoriented round-1 file still contains.

**Fix:** Parse both files. If round 2 is missing a standard gene that the shifted round-1 file has, keep round 1 and record that choice in `reannotation_decision.tsv`.

### 3.6 medium — a failed coordinate shift keeps the old location on the new sequence

**Where:** `mitom_orient.py` 103–118.

```python
new_seq = record.seq[start_pos:] + record.seq[:start_pos]
...
for feat in record.features:
    try:
        new_loc = shift_location(feat.location, start_pos, seq_len)
        new_features.append(...)
    except Exception:
        new_features.append(feat)
```

**Why:** The sequence has already been rotated (or, on the reverse-complement path, replaced) before the per-feature try. On any exception the original feature, with the original coordinates, is appended to the new record. There is no log line and the record is still written (`294–295`). The gene is then on the wrong bases.

**Input:** A feature whose location uses a position class that `shift_location` does not construct with one integer (`mito_common.py` 348–350 does `type(old)(new)` for `BeforePosition` and `AfterPosition` only), or whose qualifier shift raises inside `move_position_qualifiers`. The rotated sequence is written with that feature's old coordinates.

**Fix:** If any feature fails to shift, do not write the rotated record. Leave that record unchanged and add it to the reorient-failure list (see 3.4). Do not copy the old feature onto the new sequence.

### 3.7 medium — FASTA input never triggers the linear-mtDNA taxon rule

**Where:** `mitom_stage.py` 50–55; `mito_taxa.py` 53–60; `mitom_topology.py` 97–107.

```python
def sample_lineage(si) -> str:
    if not si.orig_gb_path:
        return ""
```

```python
LINEAR_MT_TAXA = ("Hydrozoa", "Scyphozoa", "Cubozoa", "Staurozoa", "Calcarea")
```

```python
if unlabelled and label == "circular":
    ...
    if complete_set(genes):
        return "circular", f"no topology label, ..."
```

**Why:** `linear_mtdna` only matches those names in a semicolon-separated lineage. FASTA has no original GenBank path, so the lineage is `""` and the taxon rule never fires. An unlabelled FASTA sequence with a complete gene set and no end-to-end gene is then called circular (lines 100–105) and can be rotated. Jellyfish, hydrozoan, and calcareous-sponge mitochondrial chromosomes are linear; rotating one moves a real end into the middle.

The GenBank path is different: a lineage that contains those taxa is set to linear unless a gene actually crosses the ends (`mitom_topology.py` 88–99). That order matches the docstring of `linear_mtdna` ("linear unless a gene runs across its ends"). Do not reorder it. The hole is the empty lineage.

**Input:** A FASTA of a complete *Hydra* or *Aurelia* mitochondrial chromosome, no `topology=` in the header, no taxonomy. The sequence is long enough and complete enough to pass `complete_set`.

**Fix:** Without a lineage, do not call a sequence circular only because the gene set looks complete. Require `topology=circular` in the header, or a gene that actually crosses the ends. Document that Hydrozoa / Scyphozoa / Cubozoa / Staurozoa / Calcarea FASTA should carry `topology=linear`.

### 3.8 medium — a GenBank record with no topology is stored as circular

**Where:** `mitom_inputs.py` 333–336; the label is then returned unchanged by `mitom_topology.py` 107 when nothing at the ends contradicts it.

```python
topo = (topology_from_header(rec.description) if fmt == "fasta"
        else rec.annotations.get("topology", "circular").lower())
if topo not in ("circular", "linear"):
    topo = "circular"
```

**Why:** Biopython omits `annotations["topology"]` when the LOCUS line has neither word. The default filled in here is `circular`. For a GenBank record the value is a label, not the FASTA "unlabelled" branch, so `decide()` returns "as labelled" when no gene is joined or cut at the ends. A linear draft is then eligible for rotation (`rotate_record` only refuses when the topology string is `linear`, `mitom_orient.py` 89–102).

**Input:** A linear mitogenome draft saved as GenBank with a LOCUS line that does not say `linear` or `circular`, and no gene running off either end.

**Fix:** If the GenBank record has no topology, do not invent `circular`. Leave it unset and run the same end-evidence decision used for unlabelled FASTA (`mitom_topology.py` 100–106), which calls a complete untouched gene set circular and everything else linear.

### 3.9 medium — a failed final recount still publishes the GenBank file

**Where:** `mitom_output.py` 288–297.

```python
try:
    pcg, trna, rrna = counts_from_genbank(draft_gb)
    orientation_ok = starts_with_gene(draft_gb, which_gene_first, force_forward)
except Exception:
    pass
os.replace(str(draft_gb), str(dest_gb))
```

**Why:** `counts_from_genbank` and `starts_with_gene` share one `try`. On any exception the draft is still moved into the final path. The returned PCG / tRNA / rRNA counts stay at whatever an earlier step stored (they are initialized earlier in the function and may already have been set when IDs were rewritten). `orientation_ok` can stay at the value from the previous successful step if `starts_with_gene` is the call that raised. `summary.tsv` can then disagree with the file that was just published, and a draft that did not parse is still the delivered GenBank.

**Input:** A draft whose GenBank parse raises (truncated feature table, illegal location). `os.replace` still runs.

**Fix:** If this block raises, do not replace `dest_gb`. Leave the draft aside, append the sample to `final_write_failed`, and return the failure to the exit-code check in 3.4.

### 3.10 medium — on Windows every genome map fails, and the run can still exit 0

**Where:** `mitom_steps.py` 386–396.

```python
with mp.get_context("fork").Pool(...) as pool:
    for k, d, msgs in pool.imap(_plot_file, jobs):
        ...
except Exception as exc:
    logging.error(f"[PLOT] failed: {exc!r}")
    ok = False
```

**Why:** Windows has no `fork`. `get_context("fork")` raises `ValueError` before any map is drawn. The exception is caught and `ok` is false, so the step is recorded as failed, but `workflow_run` does not look at `maps.ok` when it chooses the exit code (3.4). A Windows run with `--plot` left on produces no plots and returns 0.

**Input:** `mitom run` on Windows with the default plot step and at least one GenBank to draw.

**Fix:** Use `spawn` when `fork` is not available (`multiprocessing.get_all_start_methods()`). Keep the pool, or draw in-process if there is a single job. Combine with the exit-code change in 3.4 so a failed plot step is visible to a script.

### 3.11 low — duplicate counters do not match the per-record reason

**Where:** `mitom_inputs.py` 377–383 and 470–473.

```python
else:
    kind = "dup_version" if meta.version_int < prev.version_int else "dup_exact"
    stats["dup_exact"] += 1
    stats["dup_exact_list"].append(...)
...
dropped.append((meta, "older version" if kind == "dup_version" else "same accession listed twice", ...))
```

```python
n_dup = n_total - n_kept - n_no_seq
```

**Why:** When the later record is the older version, the kept accession is correct and the row reason can say `older version`, but the counter and the list that the log prints are `dup_exact`. The summary line therefore says "same accession listed twice" for an older version. Separately, `n_dup` treats every record that is not kept and not sequence-less as a duplicate, which includes `parse_errors`. A truncated file is counted as a duplicate in the log line at 475–479.

**Input:** The same accession twice, the older version read second. And a file that fails to parse, in the same run. Compare `dedup_report.tsv` column `reason` with the `dup_exact=` / `dropped_dup=` figures in the log.

**Fix:** Increment `stats[kind]`, and append to `stats[f"{kind}_list"]`. Subtract `parse_errors` (and any other non-duplicate drop) inside `n_dup`.

### 3.12 low — the merged GenBank file is not written via a temporary file

**Where:** `mitom_output.py` 76.

```python
merged_out.write_text("".join(chunks))
```

**Why:** `Path.write_text` truncates the destination first. A kill or a full disk leaves a partial `all.final.gbf`, which a later step can pick up as `final_annotation_all.gb`. The same direct `SeqIO.write` into the reorient output is at `mitom_orient.py` 294–297. One record's final file uses `os.replace` (`mitom_output.py` 296); these two writers do not.

**Input:** Interrupt `mitom run` while it is writing `all.final.gbf`.

**Fix:** Write the merged text and the reorient GenBank to a temporary file in the same directory, then `os.replace` it into place.

### 3.13 low — an origin-split `/transl_except` is not moved when annotations are compared

**Where:** `mito_compare.py` 111–121.

```python
m = re.match(r"\(pos:(complement\()?(\d+)(?:\.\.(\d+))?\)?,(.*)$", text)
if not m:
    return text
```

**Why:** A one-base or simple range is rewritten. `pos:join(999..1000,1..2)` does not match, so the qualifier is left on the pre-rotation coordinates while the feature location is moved. The gene-location columns of the comparison table are computed separately and can still be right; the lifted qualifier is not. `mito_common.move_positions` / the join logic already used for feature locations is the implementation to reuse.

**Input:** Compare two circular records after a rotation, with a CDS `/transl_except="(pos:join(16500..16560,1..2),aa:TERM)"`.

**Fix:** Parse `join(...)` and run each number through the same move used for a simple position. If the parse fails, drop the qualifier or mark the row, instead of copying stale coordinates.

---

## Checklist

Tick these in order. High items first. After each fix, run the existing test named on that line if there is one, and add a regression for the input in the item.

**High**

- [ ] 1.1 Trematoda → code 21, not via Platyhelminthes → 9
- [ ] 1.2 Tunicate / Urochordata → 13, lancelet → 5; `--clade Chordata` must not override them
- [ ] 1.3 Curation writes the code it used into `/transl_table`; no silent fallback to 5
- [ ] 1.5 Feature table emits a valueless `pseudo`, not `pseudo	None`
- [ ] 2.1 `coding_frame` does not move a short CDS on two stops; NC_086642.1 nad5 still moves
- [ ] 2.2 Intron limits use each `location.parts` entry, including origin wraps and genes inside introns
- [ ] 2.3 `gene_end` stops at an early stop, does not walk 300 nt into the next gene, and does not cut `TAT` / `TATG`
- [ ] 2.4 rDNA query-minus (`qstart > qend`, `sstart < send`) is strand `"-"`; `test_rdna_complete_unit_reverse_strand`
- [ ] 3.1 Two-piece circular cox1 rotates to the 5′ piece and is sent to round 2
- [ ] 3.2 Origin-crossing intron notes are removed only when that CDS got an intron feature
- [ ] 3.3 RefSeq dedup keeps an INSDC record whose sequence differs

**Medium**

- [ ] 1.4 `validate_record` puts the same code in `[mgcode=]`
- [ ] 1.6 Code 14 is possible; existing `/transl_table=14` is never replaced by 9
- [ ] 1.7 Rhabdopleuridae → 24, Cephalodiscidae → 33
- [ ] 1.8 Submission names use NCBI trnL1/L2 and trnS1/S2; MITOS2 numbers stay internal
- [ ] 1.9 Bare `atp` is not renamed to ATP6
- [ ] 1.10 `/transl_except` keeps one-base and `join` forms
- [ ] 1.11 `pcg_check` uses `/codon_start`
- [ ] 1.12 Apicoplast is not written as chloroplast
- [ ] 1.13 `mitom fetch` matches `complete mitochondrial genome`
- [ ] 2.6 atp8 / nad4L fallback searches across the origin
- [ ] 2.7 `--stop_rule homolog` does not shorten a sense codon
- [ ] 2.8 ARWEN / ARAGORN anticodons are kept
- [ ] 3.4 Exit 1 if any record failed to annotate or to write
- [ ] 3.5 Round 2 does not replace a better round 1
- [ ] 3.6 A shift exception does not publish old coordinates on a new sequence
- [ ] 3.7 FASTA without a lineage is not called circular only for a full gene set
- [ ] 3.8 Missing GenBank topology is decided from the ends, not defaulted to circular
- [ ] 3.9 A failed recount does not `os.replace` the final GenBank
- [ ] 3.10 Windows plots use `spawn`

**Low**

- [ ] 1.14 *Radopholus* is not assumed to be code 5
- [ ] 1.15 `subunit 1` maps to COX1 in curation
- [ ] 2.5 tRNA half-overlap uses `end - start + 1`
- [ ] 3.11 `dup_version` is counted as `dup_version`; parse errors are not duplicates
- [ ] 3.12 Merged GenBank and reorient output use `os.replace`
- [ ] 3.13 Comparison lifts `pos:join(...)` in `/transl_except`
