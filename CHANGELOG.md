# Changelog

Notable changes of MitoM. Versions follow semantic versioning (major.minor.patch); `mitom --version` prints the
installed one.

## 1.0.0 - 2026-10-28

### Added
- `--engine mitoz|mitos2` (`mitom run`, `mitom annotate`): MitoZ 3.6 or MITOS2 2.1.10 find the genes when they are
  installed (they are not installed with MitoM); every later MitoM step runs on their annotation (reorientation,
  round 2, curation and NCBI names, NCBI submission files, genome maps). `mito_engines.py`.
- `mitom import`: existing MitoZ results (result folders or .gbf files), MITOS2 result folders (result.bed and the
  genome FASTA) or any annotated GenBank -> the same MitoM steps, without annotating again.
- Circular genome maps (`mito_plot.py`, `--plot`, default yes): genes by strand and class, GC content, GC skew; read
  depth with `--infq` (bwa; minimap2 when it is installed, about 8x faster).
- `mitom setup`: downloads the MITOS2 reference data refseq89m (md5 checked) and writes the config file.
- `mitom --version`, `CHANGELOG.md`, tests on GitHub Actions (`.github/workflows/tests.yml`), a bioconda recipe
  (`packaging/bioconda/`).
- Introns from the intron sites of annotated relatives (`mito_intron_ref.py`, run after `mito_intron.py`; report
  `intron_ref_report.tsv`): the site library `models/intron_sites.*` (79 KB) is built from a gold standard of all 193
  local NCBI records with introns, each NCBI intron checked (spliced product, junction against homologs, group I / II
  splice rules, related records): 309 introns, 97 placed differently from NCBI (`benchmark/results/
  intron_gold_2026-10-29/`). `mitom run` against it: 309 / 309 exact with the whole library, 299 / 309 when the
  record is left out (a new genome), no false intron on 76 intron-free records, no NCBI ERROR.
- `benchmark/eval_introns.py` (introns vs NCBI or the gold standard, nested introns, per gold class) and
  `benchmark/bench_intron_ref.py` (own / leave-one-out / leave-species / leave-genus).

### Changed
- One conda environment (`env/mitotools.yml`) holds everything, the Python packages included (Biopython,
  scikit-learn 1.7.2, matplotlib); `mitom` runs on that environment's Python and uses its programs first; nothing is
  installed into the system Python any more.
- The main program is `mitom_pipeline.py` (was `batch_mitoz.py`; the old name still works).
- Work folders `01.anno_round1/`, `03.anno_round2_reoriented/`; default output folder `mitom_out`.
- The Docker image is built from that one environment, so it includes scikit-learn (MitoM-ML tRNAs instead of the
  MiTFi fallback).
- Java (MiTFi; its classes are Java 5, so any Java >= 8 runs it): the computer's own Java is used; install.sh adds
  openjdk 8 to the environment only when there is none (openjdk 25 alone was 0.33 GB). Environment 1.6 GB, 157 packages.
- Step 6 of `run` prints what the intron step did (introns written, low-confidence candidates not written, or skipped
  because exonerate / blastp were not found - now a warning on the console, before only in run.log).

- Second copies of a protein gene are kept when they are >= 95 % identical to it (length within 3 nt): the core keeps
  a single such copy (`COPY_IDENTITY`; before only copies in a block of >= 2 genes), and curation rule D1 no longer
  removes them (before it needed both copies to end with a stop, which the core's genes lack at that point).
  NCBI annotates such copies (heron / cormorant / turaco nad6, frog cob, squid cox2): on the records concerned
  MitoM gains 9 genes NCBI has, 1 it does not (NC_049135.1 atp8, 98 %); imported MitoZ / MITOS2 results gain 9 / 7,
  all in NCBI; nothing else changed.
- `mito_intron.py`: a gene across the origin of a circular record is modelled on a rotated copy (the cox1-first
  rotation put the origin inside annelid cox1 introns, half of them were missed); exons chained in protein order only;
  a model whose intron holds a copy of its own exons is refused (duplicated genes, Hexamermis atp6); trans-spliced
  genes written as partial CDS pieces; partial core genes searched with the core's reference protein; intron ends
  refined independently by the median score over up to five homologs; introns next to an inserted mini-exon placed
  again. Against the gold standard 285 -> 296 of 309 exact on its own.
- After the intron step a circular record is rotated again so that the orientation gene's first exon starts it;
  rotation keeps partial ends (< / >; dropped before, which gave table2asn ERRORs on trans-spliced pieces).

### Removed
- MitoZ as a part of MitoM: `--core mitoz|hybrid`, `--trna`, `--circos`, `--conda_env_mitoz`, `--mitoz_path`,
  `--mitoz_extra`, `env/mitoz.yml`, `scripts/mitoz_fast.py` (MitoZ can still be the engine, see above).

## 0.9.0 - 2026-10-27
- NCBI-submittable output (no table2asn ERROR on 2,497 benchmark records), MiTFi run as MITOS2 runs it, tRNA second
  opinion (MitoM-ML with MiTFi), `--rrna_method hmm|cm`; MitoM's own annotation core `mito_core.py` as the default.
