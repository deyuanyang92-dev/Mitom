Checked 2026-10-23.

# Nuclear markers in contigs: rDNA 18S-ITS1-5.8S-ITS2-28S and histones (mito_nuc.py)

Scope: how `mitom nuc` annotates nuclear rDNA units and histone genes in assembled contigs (genome skims, SPAdes
assemblies), how it was tuned, what it was compared with, and its limits. Every bullet ends with a source:
`(path:line)`, `(https://...)`, `(doi:...)` or 未确证. [our run] = measured on 2026-10-23 on this machine (72 CPUs);
scratch files of the tuning runs are not kept, the final numbers are reproduced by `benchmark/nuc_bench.py`.

## 1. Truth sets (NCBI, metazoa, fetched 2026-10-23)

- rDNA query `18S[Title] AND 5.8S[Title] AND 28S[Title] AND "complete sequence"[Title] NOT partial[Title]` under txid33208 gave 667 records; 542 carry exactly one complete 18S, ITS1, 5.8S, ITS2 and 28S feature on one strand, 378 species after keeping the first record per species [our run] (mito_nuc.py:85-89) (benchmark/nuc_bench.py:35-62)
- The rest annotate the whole unit as one misc_RNA "contains 18S ribosomal RNA, internal transcribed spacer 1, ..." (118 records, mostly Cnidaria) or have partial ends [our run] (benchmark/nuc_bench.py:35-62)
- Taxon mix of the 378: Arthropoda and Nematoda dominate, then Chordata and Lophotrochozoa (annelids: Hesionidae, Polynoidae, Sigalionidae; no Nereididae), a few Porifera, Cnidaria, Gnathifera [our run] (models/nuc/README.md:1)
- One species can have ~100 records (Nacobbus aberrans 97); without a cap they fill BLAST's top 50 targets, and after leaving the genus out a fragment had no reference hit at all, so the library keeps <= 3 records per species [our run] (mito_nuc.py:71)
- Histone query `"histone H3/H4/H2A/H2B/H1"[Title] AND "complete cds"[Title]` gave 631 records; 677 CDS remain after dropping partial ones, > 400 aa ones and histone-modifying proteins ("histone H4 transcription factor" 7 CDS of 1.5 kb, DOT1 and other methyltransferases) [our run] (mito_nuc.py:83-84) (benchmark/nuc_bench.py:65-80)
- Submitter boundaries are not consistent: the 18S 3' end is GATCATTA in 313 of 378 records, the 28S 3' end TGATTCGT in 185, the 5.8S 5' end AAACCATT (96) or AAACCACT (78) or ACAACTCT (22); some records are shifted as a whole (ON014401.2: all three genes -302 nt; the predicted 18S start is the canonical ...TGGTTGATCCTGCCAG, the annotated one is not) [our run] 未确证 which convention is right

## 2. rDNA method and tuning

- barrnap 0.9 `--kingdom euk` (nhmmer with Rfam-derived HMMs) on the 378 records: boundaries exact 1.8 %, <= 5 nt 42.1 %, <= 20 nt 79.9 %; the 28S model ends early (28S 5' end within 20 nt in 5.6 %) [our run] (https://github.com/tseemann/barrnap) (benchmark/nuc_bench.py:113-125)
- The barrnap on bioconda today is 1.10.5, a rewrite without `--kingdom euk`; 0.9 must be pinned [our run] (https://github.com/tseemann/barrnap)
- Infernal cmsearch with Rfam RF01960 (SSU, 1831 columns), RF00002 (5.8S) and RF02543 (LSU, 3401 columns) did not finish 30 records in 10 min on 16 CPUs, so covariance models were not used [our run] (https://rfam.org/family/RF02543)
- Chosen method (MitoFinder-like reference projection): blastn of the contig against complete reference units, best collinear HSP chain per reference, each annotated reference boundary mapped column by column through the alignment, weighted vote of the TOPN = 5 best references [our run] (mito_nuc.py:271-330) (mito_nuc.py:333-420) (doi:10.1111/1755-0998.13160)
- Tuning on the 378 records, genus left out (exact / <= 5 nt): first prototype with dc-megablast and 10 votes 82.5 / 89.8; HSP chaining + 30 nt projection 82.2 / 90.3; 5 votes 83.2 / 90.8; vote weight by local identity at the boundary or dropping references that disagree with the others did not help (82.1-83.2) [our run] 未确证 beyond this set
- Projection without alignment for a boundary no reference reaches inside the contig: 45 of 2,252 boundaries, only 18 exact; kept, but the feature note and the tsv evidence column say "projected" [our run] (mito_nuc.py:333-420)
- A false 28S in a real assembly came from stacked HSPs of one short repeat (the same 48 nt of the contig hit eight places of the intergenic spacer); chains now require both sequences to advance and a gene needs >= 80 aligned nt (contig-end pieces next to a supported gene excepted) [our run] (mito_nuc.py:282) (mito_nuc.py:394-401)
- Stage 2 with megablast word 16 instead of dc-megablast: 85.2 / 92.6 vs 84.3 / 92.5, and one third of the CPU time [our run] (mito_nuc.py:441)
- Final (benchmark/nuc_bench.py, 378 records, genus left out), exact / <= 1 / <= 5 / <= 20 nt: 18S 91.3 / 94.0 / 96.0 / 97.8, 5.8S 81.7 / 87.6 / 92.6 / 95.6, 28S 82.5 / 86.4 / 89.2 / 90.3, all 85.2 / 89.3 / 92.6 / 94.6 [our run] (benchmark/nuc_bench.py:150-200)

## 3. Histones

- tblastn of the reference proteins, exonerate protein2genome per locus with the 3 best proteins: whole model exact 82.7 % of 677 CDS; most errors were C-terminal tails shorter in the reference and missed N-terminal methionines [our run] (mito_nuc.py:511-570)
- Extending in frame to the next stop and, when the model does not start at protein position 1 with ATG, to the farthest upstream in-frame ATG: 97.3 % (nearest ATG instead: 96.3 %) [our run] (mito_nuc.py:476-508)
- Also running protein2dna (no introns) and multiplying spliced scores by 0.9 per intron removes false introns in H1 tails: 97.6 % (prototype) / 97.5 % (module, protein cover filter 0.3) [our run] (mito_nuc.py:78)
- Remaining misses are divergent H1 variants (oocyte H1) whose only same-type references are in the left-out genus [our run] 未确证 whether a Pfam linker-histone HMM would find them
- No comparable tool was run for histones (barrnap does not annotate protein genes); the prototype numbers above are the baseline [our run] 未确证

## 4. Massive assemblies and consensus

- SRR16188826 SPAdes scaffolds (1,477,287 sequences, 721 Mb): prefilter (makeblastdb + megablast of one unit per genus + tblastn of one protein per gene and phylum) 80 s, 3,048 candidate contigs; whole run 2 min 18 s on 32 threads; barrnap 0.9 on the same file 39 s [our run] (mito_nuc.py:234-260)
- On that assembly barrnap reports 18S 5, 28S 3, 5.8S 2 (and 64 5S); mitom nuc 18S 9, 28S 24, histone genes 22 (the rDNA is split over many short high-coverage contigs) [our run] (mito_nuc.py:842-890)
- Single blastn processes over many short queries used ~4 of 32 CPUs; splitting the regions into one chunk per thread fixed it [our run] (mito_nuc.py:422-450)
- Consensus with MAFFT --addfragments on a divergent reference backbone spread short fragments over the alignment and dropped a 98 %-identity fragment (13.5 % identity to the consensus); replaced by overlap joining of the sample fragments plus a blastn pileup [our run] (mito_nuc.py:785-840)
- Simulation: 40 complete NCBI units cut into overlapping contigs of 300-2,500 nt (coverage 500-3,000), half with a low-coverage copy carrying 3 errors, plus a 1,500 nt fragment of another phylum at coverage 2-20; the consensus of every gene found equals the original sequence base for base (the only N are N of the original records); ITS1 and ITS2 were not found in 1 record each [our run] (mito_nuc.py:842-890)
- SPAdes neighbours share an exact end overlap of about k (25 nt seen) that is below the 40 nt blastn merge; they were joined with an N and the overlap repeated (PV200780 28S, 27 extra nt); exact end overlaps >= 15 nt are now joined directly [our run] (mito_nuc.py:839)
- GetOrganelle 1.7.7.1 comparison on simulated skims (8 species of 6 phyla; 150-copy rDNA array with 0.05 % copy variation in 3 Mb background, 5x 150 bp PE ART reads, plus a primate rDNA contaminant at 5x; seeds and references without the species' genus; GetOrganelle -F anonym with our rDNA library as seed and label database): SPAdes + mitom nuc consensus 0 aligned errors in 40 genes, GetOrganelle paths + mitom nuc 671; GetOrganelle assembled the primate contaminant instead of the sponge (AB511881.1) and the jellyfish (KY212123.1), on the other 6 species both are equal; run time per species SPAdes ~50 s + mitom 21 s, GetOrganelle ~55 s [our run] (https://github.com/Kinggerm/GetOrganelle)
- Using also NCBI records with only a complete 18S, 28S or 5.8S (and a per-boundary choice of the voting references) was started but stopped: the library keeps complete 18S-ITS1-5.8S-ITS2-28S units only; per-boundary choice on that library: exact 84.5 % vs 85.2 % [our run] 未确证 with fragment references
- Histone copies are not merged: the representative is the complete CDS of the best-covered contig, because the copies of a histone cluster are paralogs [our run] 未确证 for taxa with identical copies

## 5. Held-out test and second tuning round (2026-10-23)

- Held-out set: `txid33208 AND 18S[Title] AND 5.8S[Title] AND 28S[Title] AND 3000:60000[SLEN]` (1354 records) minus the 667 library-query records -> 898 records, 71 with all five parts complete (one per species, none of these species in the library; 21 from one oribatid-mite study LC8173xx); histone genomic records (title histone H1/H2A/H2B/H3/H4, core histone or histone gene cluster/repeat, `biomol_genomic`, NOT "partial cds", 2372) minus the library records -> 172 complete CDS in 128 records, 9 spliced [our run] (benchmark/nuc_bench.py:150-200)
- With the first library the held-out rDNA scored exact / <= 5 / <= 20 nt 41.5 / 76.8 / 82.2 vs 85.2 / 92.6 / 94.6 on the tuning set: the 28S 3' end came out 270-640 nt short in mites because related references (ticks KY4575xx, one study) end their 28S ~530 nt before the conserved ...GATTCGT end [our run] (models/nuc/README.md:38)
- Library QC by other studies (`qc_rdna`, studies = first publication title, one vote per study, tolerance 15 nt, >= 3 studies and > 2x those agreeing): 133 boundaries set not ok; tolerance 30 nt flagged 88 and left the tick 18S 5' end (19 nt late) in; choosing the consensus as a +-5 nt window instead of the exact mode let it slide toward the annotation and missed the same case [our run] (mito_nuc.py:196)
- A side whose top-5 references all have unusable boundaries now takes the next 5 references that have one; before, the gene was dropped (21 firefly 28S lost in the first QC try) [our run] (mito_nuc.py:444)
- 5.8S of divergent taxa (Myxozoa, Bursaphelenchus, Acanthocephala) had no megablast hit; the spacer between found 18S and 28S is searched again with blastn word 11 and chains >= 50 bits: held-out 5.8S <= 20 nt 82.4 -> 91.5 % (barrnap 90.8), tuning set 95.6 -> 98.0 %; ITS1/ITS2 are no longer drawn across a missing 5.8S [our run] (mito_nuc.py:503)
- exonerate protein2dna writes only a `similarity` GFF line, so the intronless histone model had never been read; with it a false ~50 nt intron in the low-complexity H1 tail (16 Limenitis H1) is gone: held-out histone exact 88.4 -> 97.7 %, tuning set unchanged 97.5 %; a simpler rule (merge any in-frame stop-free intron) broke real H3 introns (KR133446) and was dropped [our run] (mito_nuc.py:599)
- Final, exact / <= 1 / <= 5 / <= 20 nt: tuning set all 85.1 / 89.6 / 92.9 / 95.1 (18S 90.7/93.5/95.4/97.8, 5.8S 82.8/89.7/95.0/98.0, 28S 81.6/85.4/88.2/89.6); the 0.1-0.6 point loss in exact vs before is records whose own annotation the QC rejected; held-out all 42.3 / 55.4 / 82.6 / 90.6, barrnap 2.6 / 24.2 / 58.0 / 71.8 [our run] (benchmark/nuc_bench.py:150-200)
- Left: the oribatid study's own conventions (18S +2/-3, 5.8S +4, 28S end -18); Myxozoa 28S 3' end 270 nt short (no myxozoan reference; the cnidarian references end correctly at ...GATTTGT); arachnid 28S ends "CAGAATTC" in 5 other studies vs ~570 nt later in LC8173xx, unresolved [our run] 未确证

## 6. Gold standard: curated RefSeq rRNA (2026-10-24)

- Gold = RefSeq NR_ rRNA of the same species mapped onto GenBank records (blastn >= 98 %, gene end aligned); 2470 GenBank records of 11 species fetched, 558 gold boundaries in 336 records [our run] (benchmark/nuc_gold_bench.py:1)
- RefSeq NR_ records are not all consistent: mouse 18S NR_003278.3 starts ACCTGG and ends GATCATTAA (1 nt shifted at both ends vs human, rat, macaque, rabbit, cow), mouse 28S ends 13 nt before human/rat, hamster 18S has 4 extra 5' nt, pig "18S" is 2302 nt, rat 5.8S ends 1 nt before human/macaque/cow/hamster; these are left out of gold and of the curated references [our run] (models/nuc/README.md:32)
- Curated references as trusted voters (a curated reference aligned over a boundary decides it) and 500 BLAST targets also in normal mode (short 5.8S-only references fell out of the top 50): genus left out, all gold boundaries exact 76.0 -> 85.8 %, <= 1 nt 91.5 -> 94.3 %; mammal 5.8S was +5 / -8 nt before [our run] (mito_nuc.py:463)
- Head-to-head on the 53 gold boundaries NCBI submitters annotated: submitters exact 96.2 %, mitom genus left out 79.2 % (C. elegans and Drosophila without a curated relative), mitom full library 100 % (circular: the species' own curated record is in the library); submitters annotate only 9.3 % of the gold boundaries, mitom 95.5 % [our run] (benchmark/nuc_gold_bench.py:1)
- Real assembly SRR16188826: the 28S consensus was 4714 nt ending in spacer (...GACGTCC); with the curated library it is 3768 nt ending at the conserved ...GAGATTTGT; run time 2 min 23 s [our run] (mito_nuc.py:592)
- Submitter-truth sets with the curated library: tuning set exact 85.1 -> 83.2 % (curated ends differ from some submitter conventions), held-out unchanged 42.3 %; histones unchanged (97.5 / 97.7 %) [our run] (benchmark/nuc_bench.py:150-200)
- Histone disagreement X57138.1: NCBI product "histone H2A.2" on a CDS whose translation MPDPAKSAPAPKKGSKKAV... is H2B (title "H2B.2 and H2A.1 genes"); mitom calls it H2B; a curated-RefSeq histone gold was not run [our run] 未确证 beyond this record

## 7. Limits

- Boundaries follow the conventions of the references; a group whose submitters use another convention (e.g. the 5.8S of the PQ4078xx chordate study, -5 / +8 nt) is off by that amount [our run] 未确证
- Phyla with few references (Porifera, Cnidaria, Ctenophora, Platyhelminthes) depend on projection from distant units; add own annotated records with `mitom nuc build --gb ... my_refs.gb` [our run] (models/nuc/README.md:1)
- Contigs without SPAdes `cov_` in the name all weigh 1, so the consensus cannot separate a contaminant from the sample by coverage, only by identity where they overlap [our run] (mito_nuc.py:80)
