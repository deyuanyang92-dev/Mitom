# Experimental evidence for mature animal mt-rRNA 5′/3′ ends: literature and data survey for MitoM benchmarking

Compiled 2026-10-26. Companion table (machine-readable): `rrna_end_anchors.tsv` (this directory)

**How things were checked**
- Every DOI in this file was resolved with the Crossref REST API (`api.crossref.org/works/<doi>`) on 2026-10-26 and the returned title/journal matched. Items I could not verify are marked **not verified**.
- PMIDs/PMCIDs and abstracts come from the Europe PMC REST API. Open-access full texts were read through Europe PMC `fullTextXML` or PMC HTML.
- GenBank coordinates were fetched with NCBI E-utilities (efetch, `gbwithparts`) on 2026-10-26.
- Several key older papers are scanned or paywalled: Elliott & Jacobs 1989, HsuChen 1984, Clary & Wolstenholme 1985, Okimoto 1992/1994, Cantatore 1990, Van Etten 1980/1983, Dubin 1982, Bakke & Johansen 2002, Kyriakou 2014, Neira-Oviedo 2011 and Cameron 2014. For these I used only the abstract plus the authors' GenBank records. PMC's PDF download for the scanned papers sits behind a browser proof-of-work check, and I did not try to get around it. Any number taken from these papers should be checked against the PDF before it is used as truth.

**Provenance codes used for coordinates**
- **T**: number copied from the paper's text or table.
- **F**: read by me from a paper figure.
- **A**: a relative statement in the abstract (for example "extends 7–13 nt into tRNA-Glu"), turned into coordinates by me using that accession's GenBank features.
- **M**: the authors' deposited RNA, cDNA or clone sequence, which I mapped onto the RefSeq/GenBank genome with blastn plus local alignment.
- **G**: current GenBank value only.

---

## 0. Key take-aways

1. **Non-mammal experimental data at roughly 1–10 nt resolution exist for 13 species. These can serve as a truth set:**
   - sea urchins: *Strongylocentrotus purpuratus*, *Paracentrotus lividus*
   - Atlantic cod: *Gadus morhua*
   - ascidian: *Halocynthia roretzi*
   - ticks: *Ixodes ricinus*, plus EST-supported 3′ ends for about 10 other tick species
   - mosquito: *Aedes albopictus*
   - *Drosophila yakuba* and *D. melanogaster*
   - the hemipteran *Erthesina fullo*
   - the tick *Dermacentor silvarum* (small-RNA-seq)
   - nematodes: *C. elegans* and *Ascaris suum* (5′ ends only)
   - the mussel *Mytilus galloprovincialis*

   Mouse and human serve as calibration. I found nothing for annelids, birds or amphibians. For cnidarians there are only qualitative RNA-seq coverage profiles.
2. **5′ ends almost always follow the punctuation rule.** The rRNA starts at the nucleotide after the upstream tRNA or gene, sometimes with a 1-nt offset:
   - cod 12S major 5′ end at +1 nt
   - mouse 12S has a 1-nt spacer after tRNA-Phe
   - cod 16S 5′ end at either the last tRNA-Val nt or the next nt
3. **3′ ends often do not follow the rule.** They:
   - **extend into the downstream gene**:
     - *S. purpuratus*: 12S runs 7–13 nt into tRNA-Glu, and 16S runs 3–5 nt into the cox1 coding region.
     - *P. lividus*: overlaps of 5 and 10 nt.
     - *I. ricinus*: rrnL runs 14–19 nt into trnL(CUN).
     - Mouse 16S: 0–7 nt into tRNA-Leu(UUR).
     - Cod 16S minor form: 5 nt into tRNA-Leu.
   - **stop well short of it** where termination motifs are present:
     - *Halocynthia* 16S ends 66 nt before the annotated end, at a TGGCAKA "termination box".
     - *Mytilus* l-rRNA ends 47–48 nt before the previous annotation.
     - Tick rrnL ends at the "Tick-Box" motif.
   - **carry non-templated oligo(A)**. The 12S tail is usually short (0–5 A in mouse, mosquito, cod and *Drosophila*). The 16S tail is longer, with oligo- or poly(A) in mouse, cod, mosquito, ticks, ascidians, *E. fullo*, zebrafish and *Drosophila*.

   An exact 3′ match to a contiguous-with-tRNA model will therefore be wrong by several nt in several lineages, even when the GenBank record is "right".
4. **Many of these GenBank records already encode the experimental ends.** Examples:
   - *S. purpuratus* X12631: 12S overlaps trnE by 13 nt and 16S overlaps cox1 by 5 nt.
   - *I. ricinus* JN248424.2: rrnL carries `/experiment="3' RACE"`.
   - Cod X99772 and *Aedes albopictus* NC_006817.

   Others do not:
   - *H. roretzi* NC_002177 is still 66 nt too long at the 16S 3′ end.
   - *D. melanogaster* NC_024511.2 12S is 2 nt longer than the later sRNA-seq annotation.
5. **Ordinary random-primed Illumina RNA-seq gives an accurate 5′ end but not a 3′ end.** The 5′ end is the pile-up of sense-read 5′ termini. The 3′ end is fixed precisely only by reads that cross the RNA–poly(A) junction. Oligo-dT-based data under-represent 12S (Gissi & Pesole 2003; Torres et al. 2009).
6. **No cross-species database of experimentally supported mt transcript ends exists.** The nearest substitutes:
   - GenBank features with `/experiment` qualifiers or transcript features
   - the Gao lab "precise annotation" records (sRNA-seq, single-nt; MK374364, MN347015, NC_024511 corrections) and their SRA data
   - raw public 3′-end/poly(A) datasets (Nano3P-seq PRJEB53494 for zebrafish and mouse; PAL-seq, 3P-seq)
7. **Inconsistency literature.** None of the papers I found measures rRNA end disagreement among near-identical GenBank records as you did.
   - Related work quantifies tRNA overlaps in RefSeq (more than 3,700 cases), gene overlaps >10 bp (2,085 cases), and general error rates.
   - MITOS2 tolerates rRNA overlaps by design: the Galaxy wrapper default is `--maxrrnaovl 50` ("Allow that a tRNA/rRNA overlaps with another feature by this number of nucleotides").

---

## 1. Papers that experimentally mapped mature mt-rRNA ends

### 1.1 Master table: coordinates on public accessions

All coordinates are 1-based, on the accession given. "5′" means the first nt of the mature rRNA. "3′" means the last templated nt, with any non-templated tail listed separately. The strand column gives the rRNA strand. A machine-readable version with more columns is in the TSV.

| Taxon / species | Accession (strand) | Gene / end | Experimental position | Relation to neighbour | Tail | Method | Source (where) | Prov. |
|---|---|---|---|---|---|---|---|---|
| Echinoid *S. purpuratus* | X12631 = NC_001453 (+) | rrnS 5′ | **74** | butt-joined to trnF (1–73) | – | S1 mapping, primer extension | Elliott & Jacobs 1989 (abstract) | A |
| | | rrnS 3′ | **945–951** (GenBank 951) | 7–13 nt **into trnE** (939–1009) | n.s. | S1 | same | A |
| | | rrnL 5′ | **4260** | butt-joined to nad2 mRNA (CDS 3201–4259) | – | S1/PE | same | A |
| | | rrnL 3′ | **5787–5789** (GenBank 5789) | 3–5 nt **into cox1 CDS** (5785–); 8–10 nt overlap with COI mRNA 5′ end (GenBank `precursor_RNA` "COI mRNA" 5778–7348) | n.s. | S1 | same | A |
| Echinoid *P. lividus* | J04815 = NC_001572 (+) | rrnS 3′ | **956** (GenBank 956) | overlaps trnE (952–1020) by **5 nt** | n.s. | transcript mapping | Cantatore et al. 1990 (abstract) | A |
| | | rrnL 3′ | overlaps COI mRNA by **10 nt**; GenBank end 5840 (COI CDS 5838) | overlap | n.s. | same | same | A (exact nt not derivable) |
| Ascidian *H. roretzi* | AB024528 = NC_002177 (+) | rrnL 3′ | **13099** (118/151 poly(A) ESTs); minor 13100 (8), 13095–13098 (18; A-run ambiguity), 13093–13094 (4), 13091–13092 (1), 13089–13090 (2). **GenBank 13165** | ends just before the termination box TGGCATA (13100–), 66 nt before the annotated 16S end / nad1 start (13166) | poly(A) | poly(A)-start sites in oligo-dT ESTs | Gissi & Pesole 2003, Fig. 2 + Results text | F + T |
| | | rrnS | no ESTs: "The only mitochondrial gene absent from both the *Halocynthia* and *Ciona* EST collections was the 12S rRNA" (likely not polyadenylated) | | | | same, Discussion | T |
| Teleost *G. morhua* | X99772 = NC_002081 (+) | rrnS 5′ | **70** (454: 46/46 reads; also PE); 69 (PE only). GenBank 69 | trnF ends 68, so the major form starts 1 nt after the punctuation site | – | primer extension + 454 read 5′ ends | Coucheron et al. 2011 text + Fig. 2B (PE from Bakke & Johansen 2002) | F + M |
| | | rrnS 3′ | **1018** + 1 non-templated A | abuts trnV (1019) | mono-A | RT-PCR sequencing | Coucheron 2011 Fig. 2B; Bakke & Johansen 2002 abstract ("homogenous and mono-adenylated") | F + M |
| | | rrnL 5′ | **1090 or 1091** (GenBank 1091) | trnV annotated 1019–1090 | – | PE + 454 | Coucheron 2011: Fig. 2B puts 454 at 1090 and PE at 1091; the text says "119 of the 120 mtLSU rRNA reads mapped one position downstream that of the PE site" (inconsistent) | F + M |
| | | rrnL 3′ | **2759** (+AAAAA, abuts trnL(UUR) 2760) and **2764** (+AAAA, 5 nt into trnL(UUR)). GenBank 2759 | | oligo-A | RT-PCR sequencing | Coucheron 2011 Fig. 2B. The legend says only "L-rRNA-2" is tRNA-processed and "L-rRNA-1" is an mTERF-terminated H1 transcript; in the figure, L-rRNA-1 ends at 2759 and L-rRNA-2 at 2764, so the labels look swapped | F + M |
| Tick *I. ricinus* | JN248424.2 (−) | rrnL 3′ | **12032–12037** (most clones, "14–19 bp inside trnL(CUN)"); 12039–12040 (1 clone, "11–12 bp inside"). GenBank rrnL 12032–13281 with `/experiment="3' RACE"` | **overlaps trnL(CUN)** (11990–12050); ends just before the Tick-Box | poly(A) | 3′ RACE (6 clones) + 13 poly(A) ESTs | Montagna et al. 2012, Results text, Table 3, Fig. 4 | T + A |
| Tick *I. scapularis* | AB161439 | rrnL 3′ | "two sites corresponding exactly to those found in *I. ricinus*" (8 poly(A) ESTs) | as above | poly(A) | ESTs | Montagna 2012 Table 3, Fig. 4 | T (nt only in Fig. 4) |
| 9 other ticks (Amblyomma ×3, Hyalomma ×2, *Dermacentor andersoni*, *Boophilus microplus* AF110619, *Rhipicephalus sanguineus* NC_002074, *Ornithodoros*) | | rrnL 3′ | Most poly(A) ESTs end at the 5′ end of the Tick-Box. In Metastriata this is "a few bp upstream of the 3′-end of the currently annotated rrnL". A minority end at the 5′ end of the downstream trnL (punctuation) | | poly(A) | ESTs | Montagna 2012 Table 3 (counts), Fig. 4 (sequences) | T |
| Tick *D. silvarum* | MN347015 (UNVERIFIED record, no features) (−) | rrnL | **6965–8185** | contiguous with trnL (6904–6964) and trnV (8186–8247) | – | sRNA-seq 5′/3′-end small RNAs + PacBio | Chen et al. 2020, Table 2 | T |
| | | rrnS | **8248–8949** | contiguous with trnV and CR1 (8950–9258) | – | same | same | T |
| Mosquito *Aedes albopictus* (C6/36 cells) | clone X01078 (rrnL 61..1395) mapped to NC_006817 (−) | rrnL 5′ | **14017** (= X01078 pos. 61) | abuts trnV (14018). X01078 has a 1-nt spacer (pos. 60) that is an indel relative to NC_006817 | – | S1 protection + RNA sequencing | HsuChen, Kotin & Dubin 1984 (abstract + their GenBank X01078 features) | M |
| | | rrnL 3′ | **12685** (= X01078 pos. 1395) | abuts trnL(CUN) (12684) | poly(A); "remarkable degree of homogeneity" | S1, RNA sequencing, 3′-end labelling | HsuChen 1984; Dubin & HsuChen 1985 (abstracts) | M |
| | RNA sequence M32716 mapped to NC_006817 (−) | rrnS 3′ | **14090** | abuts trnV (14089, first nt C) | 3/4 of molecules end "…GA", 1/4 "…GAA" (+1 non-templated A) | direct RNA sequencing of the 3′-terminal 101 nt | Dubin & HsuChen 1983 (abstract); GenBank M32716 | M (2 mismatches in the 20-nt end segment) |
| *Drosophila yakuba* | X03240 = NC_001322 (−) | rrnS 5′ | GenBank **14942**. The paper says "The 5′ end of the small rRNA gene was located by S1 protection analysis"; the coordinate is not checked in the paper | next to the A+T region; "a tRNA gene was not found at the 5′ end" | – | S1 | Clary & Wolstenholme 1985 (abstract) | G |
| *D. melanogaster* | cDNA X53506 mapped to NC_024511.2 (−) | rrnL 5′ / 3′ | 5′ **14058**; 3′ **12735** (cDNA ends …ATAAAA, so 12734 is ambiguous) | the 5′ end is the last nt of RefSeq trnV (14058–14130), a 1-nt overlap; the 3′ end abuts trnL(CUN) | – | full-length cDNA clone | Kobayashi & Okada 1990 (GenBank X53506) | M (cDNA ends need not be mature ends) |
| | NC_024511.2 (−) | rrnS | **14131–14914** (RefSeq still 14131–14916) | 5′ end abuts the new control-region start 14915 | – | sRNA-seq + PacBio + ONT cDNA | Liang et al. 2022 (bioRxiv v1 text: "modified from 14131–14916 to 14131–14914") | T (read via WebFetch; Table 1 not accessed) |
| | | tails | 12S "retained on average four adenines on full-length transcripts"; 16S oligoadenylated. Torres et al.: 12S has "a very short poly(A) tail (approximately 6-nt long)" | | | cRT-PCR; 454 | Bratic et al. 2016 text; Torres et al. 2009 text | T |
| Hemipteran *Erthesina fullo* | MK374364 (−) | rrnL | **12466–13742** | contiguous with trnL (12401–12465) and trnV (13743–13810) | poly(A) | sRNA-seq end small RNAs (+PacBio) | Ji et al. 2019, Table 1 (identical to the MK374364 features) | T |
| | | rrnS | **13811–14624** | contiguous with trnV; 5′ end abuts the control region (14625–) | poly(A) | same | same. The text notes that a 15-bp 12S 5′ segment "was missed in JQ743673.1" | T |
| | own PacBio consensus genome (Gao 2016, Suppl. file 1; not GenBank) | rRNAs | 16S 12,464–13,742 (overlapping trnL by 2 nt in that table). Three 12S isoforms share the 5′ end 14,624, with 3′ ends 13,814 / 13,865 / 14,051. "full-length polyadenylated transcripts of two rRNAs" | | poly(A) | PacBio Iso-Seq | Gao et al. 2016, Table 1 + Results | T |
| Nematode *C. elegans* | X54252 = NC_001328 (+) | rrnS 5′ | **898** | adjacent to trnE (842–897) | – | primer extension | Okimoto et al. 1994 (abstract) | A |
| | | rrnL 5′ | **10403** | adjacent to trnH (10348–10402) | – | primer extension | same | A |
| Nematode *A. suum* | X54253 = NC_001327 (+) | rrnS 5′ / rrnL 5′ | **901** / **11346** | adjacent to trnE (842–900) / trnH (11291–11345) | – | primer extension | same | A |
| Mussel *M. galloprovincialis* F | AY497292.2 = NC_006886.2 (+) | rrnS 5′ / rrnL 5′ | "coincided with the first nucleotide of the annotated genes". Current GenBank: 14070 / 15501 | after trnF (14002–14069) / trnD (15436–15500) | – | RNA circularization + RT-PCR | Kyriakou et al. 2014 (abstract) | A (annotation version not verified) |
| | | rrnL 3′ | "47 or 48 nucleotides upstream from the one assigned by a previous annotation". If that annotation is the current one (16744), the 3′ end is ≈ **16696–16697** (**my inference, not verified**) | the control-region VD1 is correspondingly longer | 3′ ends "heterogeneous", both polyadenylated and non-adenylated; truncated forms present | same | same | A (inference) |
| Mouse (calibration) | NC_005089 (+); same numbering as J01420 here | rrnS 5′ / 3′ | **70** / **1024** (major), 1023 (minor) | 1-nt spacer (69) after trnF; the 3′ end abuts trnV | 12S: "1 to about 5" A | direct RNA sequencing of the termini, S1, fingerprinting | Van Etten et al. 1980, 1983 (abstracts) | A |
| | | rrnL 5′ / 3′ | **1094** / **2675–2682** | 3′ end "from the nucleotide immediately adjacent to the 5′-end of the tRNALeuUUR gene to any position up to 7 nucleotides downstream within the tRNALeuUUR gene" | oligo(A), "somewhat greater" than 12S | same | same | A |

*n.s.* = not stated in the abstract.

### 1.2 Notes per study (non-mammals first)

**Echinoderms**
- **Elliott DJ, Jacobs HT (1989).** Mutually exclusive synthetic pathways for sea urchin mitochondrial rRNA and mRNA. *Mol Cell Biol* 9:1069–1082. doi:10.1128/mcb.9.3.1069-1082.1989 (verified); PMC362697 (scanned only).
  - Species: *S. purpuratus* embryos.
  - Methods: RNA blot hybridization, S1 mapping, primer extension.
  - Abstract (verbatim): "although they appeared to be butt-joined at their 5′ ends to the upstream transcripts, tRNA-Phe 5′ to the small subunit (12S) rRNA and NADH dehydrogenase subunit 2 mRNA 5′ to the large subunit (16S) rRNA, respectively, their 3′ ends were found to overlap the 5′ ends of the downstream transcripts. 12S rRNA was found to extend 7 to 13 nucleotides into the sequence of tRNA-Glu; 16S rRNA was shown to terminate 3 to 5 nucleotides inside the coding region of cytochrome oxidase subunit 1 (COI) and 8 to 10 nucleotides from the mapped 5′ end of COI mRNA."
  - The GenBank record X12631 already uses the longest form (12S 74–951; 16S 4260–5789). It also carries transcript features ("COI mRNA" `precursor_RNA` 5778–7348, "ND1 mRNA start site", "ND2 primary transcript", promoters).
- **Cantatore P, Roberti M, Loguercio Polosa P, Mustich A, Gadaleta MN (1990).** Mapping and characterization of *Paracentrotus lividus* mitochondrial transcripts: multiple and overlapping transcription units. *Curr Genet* 17:235–245. doi:10.1007/bf00312615 (verified).
  - Abstract: "The 3′ ends of the two rRNAs (12S and 16S) overlap with the downstream transcripts (tRNAGlu and CoI mRNA) by 5 and 10 nt respectively", and "Several mtRNAs were found to have 5′ and 3′ termini which differ from those inferred through DNA sequencing".
  - J04815 12S (74–956) overlaps trnE (952–) by exactly 5 nt, consistent with the paper.

**Fish**
- **Bakke I, Johansen S (2002).** Characterization of mitochondrial ribosomal RNA genes in gadiformes… *Mol Phylogenet Evol* 25:87–100. doi:10.1016/s1055-7903(02)00220-8 (verified; paywalled).
  - Abstract: "The 5′ and 3′ ends of Atlantic cod mt-rRNAs were precisely mapped. The 3′ ends of mt-SSU rRNAs were found to be homogenous and mono-adenylated, whereas that of the mt-LSU rRNAs were heterogenous and oligo-adenylated. The 5′ ends of mt-SSU rRNAs appeared to be heterogenous, corresponding to the presumed first and second positions of the gene."
- **Coucheron DH, et al. (2011).** Characterization of mitochondrial mRNAs in codfish reveals unique features compared to mammals. *Curr Genet* 57:213–222. doi:10.1007/s00294-011-0338-2 (verified); PMC3097352 (open access).
  - Fig. 2B shows all four rRNA ends at sequence level. I mapped them onto X99772/NC_002081 (results in the table).
  - The 454 reads come from random-hexamer-primed cDNA of poly(A)-enriched RNA. Text: "Pyrosequencing identified 46 and 120 reads corresponding to the 5′ ends of mtSSU rRNA and mtLSU rRNA… a single 5′ end in mtSSU rRNA (46/46 reads) corresponding to one of the PE sites. Similarly, 119 of the 120 mtLSU rRNA reads mapped one position downstream that of the PE site."
  - This is direct evidence that ordinary random-primed reads recover rRNA 5′ ends to within 0–1 nt of primer extension.
  - Also: "The ribosomal RNAs are oligoadenylated at their 3′ ends according to downstream tRNA acceptor structures".
- **Zebrafish (and mouse): Begik O, et al. (2023).** Nano3P-seq. *Nat Methods* 20:75–85. doi:10.1038/s41592-022-01714-w (verified); PMC9834059.
  - "many 16S mitochondrial rRNA reads contained a polyA tail… present not only in zebrafish… but also in mouse".
  - No coordinates are given. The raw data (ENA PRJEB53494) could be used to call 16S 3′ ends in zebrafish.

**Tunicates**
- **Gissi C, Pesole G (2003).** Transcript mapping and genome annotation of ascidian mtDNA using EST data. *Genome Res* 13:2203–2212. doi:10.1101/gr.1227803 (verified); PMC403730.
  - Results (verbatim): "In the 16S rRNA contig, there are at least six different sites in which poly(A) stretches start. All are concentrated in an AT-rich region, 11 bp long… with the most frequent poly(A) start site found in 118 ESTs of a total of 151", and "a shorter 16S rRNA gene than has been proposed, with a 3′ end located about 50 bp upstream of the annotated end (Yokobori et al. 1999) (position 13100 against the annotated position 13165)".
  - My reading of Fig. 2 puts the 118-EST site immediately before T13100, the first nt of the boxed termination motif. The last templated nt is therefore 13099, which is 66 nt short of the RefSeq end. RefSeq NC_002177 is still unchanged.
  - The paper also notes that the "rRNA transcription termination box" sits "exactly where the 16S rRNA transcript terminates" in sea urchin.

**Ticks and other chelicerates**
- **Montagna M, Sassera D, Griggio F, Epis S, Bandi C, Gissi C (2012).** Tick-Box for 3′-end formation of mitochondrial transcripts in Ixodida, basal chelicerates and *Drosophila*. *PLoS One* 7:e47538. doi:10.1371/journal.pone.0047538 (verified); PMC3471875.
  - Methods: 3′ RACE (oligo(dT) adaptor) for *I. ricinus*, plus EST poly(A)-start mapping. The rule used: "'A' or 'T' stretches >10 bp located at the end of EST sequences were considered equivalent to the polyA tail". The authors note that genomic A-runs allow mapping only "in a range of 2–5 nucleotides".
  - Results: "In *I. ricinus*, the rrnL polyadenylated transcript ends at two alternative sites, separated by 1 bp and located inside trnL(CUN), immediately before the 5′-end of the Tick-Box motif… most rrnL 3′ RACE clones stop 14–19 bp inside trnL(CUN), while only one clone stops 11–12 bp inside trnL(CUN)… these ESTs do not provide support for the existence of rrnL transcripts terminating at the 5′-end of trnL(CUN), as predicted by the tRNA punctuation model."
  - Methods also state that rrnS and other rRNA boundaries "were settled as adjacent to those of the flanking genes", so only the rrnL 3′ end is experimental.
  - Table 3 gives per-species counts of poly(A) ESTs at the Tick-Box versus "other" sites. Fig. 4 gives the sequences.
- **Chen Z, Xuan Y, Liang G, et al. (2020).** Precise annotation of tick mitochondrial genomes reveals multiple copy number variation of short tandem repeats and one transposon-like element. *BMC Genomics* 21:488. doi:10.1186/s12864-020-06906-2 (verified); PMC7367389.
  - Species: *D. silvarum*. Method: sRNA-seq end small RNAs (SRA SRP178347) plus PacBio.
  - Table 2 gives the rRNA coordinates in the master table. They are contiguous with the flanking tRNAs.
  - Caveat: the Gao-lab pipeline aims for annotations "without gaps or overlaps", and its contiguity with the tRNA-Leu boundary differs from Montagna's EST result for metastriate ticks (end at the Tick-Box). Treat this as moderately independent evidence.
- **Warren JM, Sloan DB (2021).** *NAR Genom Bioinform* 3:lqaa111. doi:10.1093/nargab/lqaa111 (verified). This is tRNA-seq of mite mt-tRNAs: tRNA ends and processing, but not rRNA.

**Insects**
- **HsuChen CC, Kotin RM, Dubin DT (1984).** Sequences of the coding and flanking regions of the large ribosomal subunit RNA gene of mosquito mitochondria. *Nucleic Acids Res* 12:7771–7785. doi:10.1093/nar/12.20.7771 (verified); PMC320199 (scanned).
  - *Aedes albopictus* cells: "located the ends of the gene by S1 protection analysis and by comparison with RNA sequences… approximately 1335 residues long."
  - GenBank X01078 (deposited by these authors) has LSU 61..1395, trnV <1..59 and trnL(CUN) 1396..1463. Mapped onto NC_006817 these become 14017..12685, exactly the current RefSeq l-rRNA.
- **Dubin DT, HsuChen CC (1983).** The 3′-terminal region of mosquito mitochondrial small ribosomal subunit RNA. *Plasmid* 9:307–320. doi:10.1016/0147-619x(83)90008-2 (verified).
  - "The extreme 3′ end of the mosquito sequence is heterogeneous, three-fourths of the molecules ending in … GA and one-fourth ending in … GAA… some, at least, of the 3′-terminal A residues may be added post-transcriptionally".
  - The RNA sequence was deposited as GenBank M32716. Mapped to NC_006817, "…GA" ends at 14090, exactly the RefSeq s-rRNA 3′ end and abutting trnV. The next genomic nt is C, so the extra A is non-templated.
- **Dubin DT, HsuChen CC (1985).** The 3′ end of large ribosomal subunit RNA from mosquito mitochondria: homogeneity of transcribed moieties. *Plasmid* 13:139–144. doi:10.1016/0147-619x(85)90066-6 (verified).
  - 16S is "post-transcriptionally polyadenylated", with "a remarkable degree of homogeneity compared to homologous mammalian mitochondrial systems".
- **Clary DO, Wolstenholme DR (1985).** The ribosomal RNA genes of *Drosophila* mitochondrial DNA. *Nucleic Acids Res* 13:4029–4045. doi:10.1093/nar/13.11.4029 (verified); PMC341294 (scanned).
  - *D. yakuba*: "The 5′ end of the small rRNA gene was located by S1 protection analysis… a tRNA gene was not found at the 5′ end".
  - This is an example of a 12S 5′ end that no tRNA defines.
- **Kobayashi S, Okada M (1990).** Complete cDNA sequence encoding mitochondrial large ribosomal RNA of *Drosophila melanogaster*. *Nucleic Acids Res* 18:4592. doi:10.1093/nar/18.15.4592 (verified). The cDNA is GenBank X53506 (mapping in the table).
- **Benkel BF, et al. (1988).** *Nucleic Acids Res* 16:9880. doi:10.1093/nar/16.20.9880 (verified). Drosophila lrRNA is an abundant polyadenylated RNA.
- **Torres TT, Dolezal M, Schlötterer C, Ottenwälder B (2009).** Expression profiling of *Drosophila* mitochondrial genes via deep mRNA sequencing. *Nucleic Acids Res* 37:7509–7518. doi:10.1093/nar/gkp856 (verified); PMC2794191.
  - Data: 454 sequencing of oligo-dT-primed cDNA. No coordinates.
  - Key statement: "the short subunit has a very short poly(A) tail (approximately 6-nt long)… prevented an efficient cDNA synthesis using poly(T) oligos". This is a warning for oligo-dT-based data.
- **Bratic A, et al. (2016).** Mitochondrial polyadenylation is a one-step process required for mRNA integrity and tRNA maturation. *PLoS Genet* 12:e1006028. doi:10.1371/journal.pgen.1006028 (verified); PMC4866704.
  - cRT-PCR in *D. melanogaster* larvae. 12S "retained on average four adenines on full-length transcripts". Tails were measured relative to annotated termini (Fig. 1F, Fig. 2D).
- **Gao S, et al. (2016).** PacBio full-length transcriptome profiling of insect mitochondrial gene expression. *RNA Biol* 13:820–825. doi:10.1080/15476286.2016.1197481 (verified); PMC5014001.
  - *E. fullo*. Both rRNAs are polyadenylated, and three 12S isoforms were found (coordinates in the table).
  - The authors also claim "possible 5′ m7G caps of rRNAs". This is unconventional; treat with caution.
- **Ji H, et al. (2019).** Using high-resolution annotation of insect mitochondrial DNA to decipher tandem repeats in the control region. *RNA Biol* 16:830–837. doi:10.1080/15476286.2019.1591035 (verified); PMC6546396.
  - sRNA-seq "5′ and 3′ end small RNAs" at "1 bp resolution". GenBank MK374364; SRA SRP174926.
  - The method is described in **Xu X, et al. (2019)**, *Front Genet* 10:105, doi:10.3389/fgene.2019.00105 (verified).
- **Liang G, et al. (2022).** Precise annotation of *Drosophila* mitochondrial genomes leads to insights into AT-rich regions. *Mitochondrion* 65:145–149. doi:10.1016/j.mito.2022.06.006 (verified). Preprint doi:10.1101/2021.06.05.447193.
  - Covers *D. melanogaster* NC_024511.2, *D. simulans* (new record MN611461) and *D. grimshawi* BK006341.1.
  - sRNA-seq sets: SRP017617, SRP199841, ERP011832, ERP012885, ERP106913. PacBio: SRP302633. ONT: ERP117816.
  - The full Table 1 (all genes) was not accessed.
- **Neira-Oviedo M, et al. (2011).** The RNA-Seq approach to studying the expression of mosquito mitochondrial genes. *Insect Mol Biol* 20:141–152. doi:10.1111/j.1365-2583.2010.01053.x (verified).
  - *Ae. aegypti* and *Anopheles*, 454 and Solexa. The abstract covers abundance only. Whether transcript ends are reported is **not verified** (paywalled).
- **Dunn E, Shao R (2026).** Transcription of the extensively fragmented mitochondrial genomes of human lice. *Biology* 15:296. doi:10.3390/biology15040296 (verified). Expression levels from public RNA-seq only; no rRNA end coordinates.

**Nematodes**
- **Okimoto R, Macfarlane JL, Wolstenholme DR (1994).** *J Mol Evol* 39:598–613. doi:10.1007/bf00160405 (verified).
  - Abstract: "Using primer extension analysis, the 5′ end nucleotides of the mt-s-rRNA and mt-l-rRNA genes were determined to be adjacent to the 3′ end nucleotides of the tRNA(Glu) and tRNA(His) genes, respectively".
  - The genomes are in **Okimoto et al. 1992**, *Genetics* 130:471–498, doi:10.1093/genetics/130.3.471 (verified).
  - I found no experimental 3′ ends for nematode rRNAs.

**Molluscs**
- **Kyriakou E, Chatzoglou E, Zouros E, Rodakis GC (2014).** *Gene* 540:78–85. doi:10.1016/j.gene.2014.01.080 (verified; paywalled).
  - F and M genomes; RT-PCR plus RNA circularization.
  - "In all these genes the 5′ ends of the transcripts coincided with the first nucleotide of the annotated genes, but the 3′ ends were heterogeneous. The l-rRNA 3′ end is 47 or 48 nucleotides upstream from the one assigned by a previous annotation… We have observed s-rRNA and l-rRNA transcripts with truncated 3′ end".
  - The F-type GenBank record AY497292.2 was last updated in 2012, before this paper, and still ends l-rRNA at 16744.
- **Ghiselli F, et al. (2013).** *Genome Biol Evol* 5:1535–1554, doi:10.1093/gbe/evt112 (verified). *Ruditapes* RNA-seq; transcription levels only, no rRNA ends.

**Cnidaria and sponges** (qualitative only)
- **Kayal E, Bentlage B, Collins AG (2016).** *RNA Biol* 13:799–809. doi:10.1080/15476286.2016.1194161 (verified).
  - *Alatina alata*, TruSeq poly(A) RNA-seq coverage: "The dramatic drop in read coverage is observed at the 3′end of rns, which can fold into a stem-loop structure". No coordinate.
- **Shimpi GG, et al. (2017).** *BMC Mol Biol* 18:16. doi:10.1186/s12867-017-0093-0 (verified). Octocorals; RACE and cRT-PCR for mRNAs only. rRNAs were the most abundant transcripts.
- **Kayal E, et al. (2015).** *PeerJ* 3:e1403. doi:10.7717/peerj.1403 (verified). Hydrozoan RNA-seq transcription model; rRNA ends "delimited… by alignment".
- **Lavrov DV, et al. (2013).** *Mol Biol Evol* 30:865–880. doi:10.1093/molbev/mss274 (verified). Calcarean sponge with fragmented rRNAs. Whether the fragment ends were RNA-mapped is **not verified**.

**Annelids**
- I found no study with mapped mt-rRNA ends.
- The closest is **Selifanova M, et al. (2023)**, *Genome Biol Evol* 15:evad219, doi:10.1093/gbe/evad219 (verified): RNA-seq coverage on the *Polydora* mitogenome, with no end analysis.

**Mammals** (calibration; besides Ojala 1981 and Mercer 2011, which you already have)
- **Crews S, Attardi G (1980).** *Cell* 19:775–784, doi:10.1016/s0092-8674(80)80053-5 (verified). Human 12S and tRNA-Phe are "joined end to end".
- **Dubin DT, Montoya J, Timko KD, Attardi G (1982).** *J Mol Biol* 157:1–19, doi:10.1016/0022-2836(82)90510-1 (verified). Precise human rRNA 3′ ends; the content is **not verified** because no abstract is in Europe PMC.
- **Van Etten RA, Walberg MW, Clayton DA (1980).** *Cell* 22:157–170, doi:10.1016/0092-8674(80)90164-6 (verified), and **Van Etten RA, Bird JW, Clayton DA (1983).** *J Biol Chem* 258:10104–10110, doi:10.1016/s0021-9258(17)44611-4 (verified). Mouse; details in the table.
- **Kobayashi M, et al. (1981).** *Gene* 16:297–307, doi:10.1016/0378-1119(81)90085-8 (verified). Rat 12S is contiguous with tRNA-Phe and tRNA-Val. Whether this rests on RNA mapping is not stated in the abstract.
- **Marková S, Filipi K, Searle JB, Kotlík P (2015).** *BMC Genomics* 16:870, doi:10.1186/s12864-015-2103-2 (verified). In the bank vole, "deep sequencing strategies specifically targeting transcript ends will be required to describe the ends of… rrnS and rrnL transcripts". In other words, ordinary poly(A) RNA-seq did not settle the rRNA ends there.
- **Pearce SF, et al. (2017).** *eLife* 6:e27596, doi:10.7554/elife.27596 (verified). Human cells, cRT-PCR plus MPAT (Fig. 5; per-clone data in Supplementary file 1): "the extensions to the 3′ end of 16S [are] usually of 4 or 5 nt in length in control PDE12+/+ cells", while 12S 3′ extensions "were generally heterogenous".
- **Kuznetsova I, et al. (2017).** *Nucleic Acids Res* 45:5487–5500, doi:10.1093/nar/gkx104 (verified). Mouse cRNA-seq. 12S rRNA is susceptible to 5′-phosphate-dependent exonuclease, i.e. it has a 5′-monophosphate, so 5′P-ligation methods (PARE, cRNA-seq) capture it.
- **Jin X, et al. (2020).** *RNA Biol* 17:395–402, doi:10.1080/15476286.2019.1709746 (verified). sRNA-seq annotation of human, chimpanzee, macaque and mouse (tables in the supplement, not accessed).

### 1.3 What this means for MitoM's punctuation rule
- **5′ ends:** the tRNA punctuation prediction is supported in sea urchin, nematodes, mosquito, mussel, *E. fullo*, ticks (sRNA-seq) and mouse/human. The main deviations are ±1 nt (cod 12S major 5′ end at +1; mouse 1-nt spacer; cod and *D. melanogaster* 16S 5′ end possibly including the last nt of tRNA-Val). With no upstream tRNA (insect and ascidian 12S next to the control region; *D. yakuba*), only RNA data can place the 5′ end. For *E. fullo* 12S, the earlier annotation (JQ743673.1) was missing 15 nt.
- **3′ ends:** the punctuation rule holds in cod, mosquito, *Drosophila* (Stewart & Beckenbach 2009, already known to you) and the sRNA-seq-based annotations. It fails by +5 to +19 nt (into the downstream gene) in sea urchins, ticks (RACE/EST), and minor mouse and cod forms. It fails by −47 to −66 nt where a termination motif sits upstream of the next gene (*Halocynthia* 16S, *Mytilus* l-rRNA). The 3′ ends are frequently heterogeneous and oligo/poly-adenylated.
- **Practical suggestions:**
  - Score 5′ and 3′ ends separately.
  - Report agreement at exact, ±1, ±3 and ±10 nt.
  - Treat experimental 3′ ends as intervals (minimum–maximum of the observed sites).
  - Do not count a 1–2-nt overlap of a 3′ end with the downstream tRNA as automatically wrong: some lineages genuinely do this.
  - The MITOS-style overlaps you see in submitter files are, however, not evidence-based.

---

## 2. Calling transcript ends from ordinary Illumina RNA-seq

### 2.1 What each library type preserves (mechanism plus evidence)

| Library / chemistry | 5′ end of the RNA | 3′ end of the RNA | Notes and evidence |
|---|---|---|---|
| **Random-hexamer-primed, fragmented (TruSeq mRNA/total, dUTP stranded)** | **Precise, in principle.** RT runs off the 5′ terminus of each RNA fragment, so every fragment containing the native 5′ end starts exactly there. In dUTP libraries the sense read (R2) starts at the RNA fragment 5′ end, giving a read-start pile-up ("cliff"). Possible 1–3 non-templated nt from MMLV-type RT (Wulf et al. 2019, doi:10.1074/jbc.RA119.010676) can show up as soft clips. | **Imprecise.** The fragment's 3′ boundary in the library is the hexamer priming site (R1 start in dUTP), not the RNA 3′ end. The native 3′ end is seen only if a hexamer primed in the last ~6 nt or in a poly(A) tail. Coverage therefore tapers over the last few nt. The exception is reads crossing the templated → non-templated poly(A) junction, which are exact. | Hexamer sequence bias in the first ~13 nt of reads: Hansen, Brenner & Dudoit 2010, doi:10.1093/nar/gkq224. Positional and fragment biases: Li, Jiang & Wong 2010, doi:10.1186/gb-2010-11-5-r50; Roberts et al. 2011, doi:10.1186/gb-2011-12-3-r22; Love et al. 2016, doi:10.1038/nbt.3682; Lahens et al. 2014 (IVT-seq), doi:10.1186/gb-2014-15-6-r86. End completeness by protocol (dUTP: 62% / 73% of yeast genes with complete 5′ / 3′ coverage): Levin et al. 2010, doi:10.1038/nmeth.1491; dUTP method: Parkhomchuk et al. 2009, doi:10.1093/nar/gkp596. Empirical mt evidence: cod 454 random-primed reads matched primer-extension rRNA 5′ ends to within 0–1 nt (Coucheron 2011). The asymmetry argument (5′ precise, 3′ smeared) is my mechanistic reading of the chemistry, not a cited measurement. Validate it on species with known ends (cod, mosquito, mouse, *Drosophila*) before trusting it. |
| **Oligo(dT)-primed cDNA / ESTs / PacBio Iso-Seq / 3′ RACE** | Depends on full-length capture. Template switching adds non-templated G's. | **Exact** at the poly(A) junction (Gissi & Pesole 2003; Montagna 2012). Genomic A-runs blur this by 2–5 nt. Internal priming at A-rich internal sites produces false ends (Gissi & Pesole: 16% of *Halocynthia* 16S ESTs). | Under-represents short-tailed or untailed RNAs, especially 12S (Torres 2009; Gissi & Pesole 2003). |
| **Poly(A)-selected RNA-seq** | as for random-primed | poly(A)-junction reads (soft-clipped A's) give exact 3′ ends | Marková et al. 2015 called mt 3′ ends from "stretches of three or more adenine residue that do not map back to the genome". They observed at least 40 tail-reads at each H-strand mRNA end, but did not resolve rRNA ends ("…will be required to describe the ends… of rrnS and rrnL transcripts"). Generic tools: KLEAT (Birol et al. 2015, doi:10.1142/9789814644730_0034); ContextMap2 poly(A)-read mapping (Bonfert & Friedel 2017, doi:10.1371/journal.pone.0170914). |
| **rRNA-depleted total RNA** | as for random-primed | as above | Probe-based depletion kits are species-specific. Commercial human/mouse/rat kits are reported to remove mt rRNA as well (**not verified here**; check the kit), whereas non-model species usually keep mt rRNA reads. |
| **Small-RNA-seq (5′P/3′OH adapter ligation)** | exact | exact | Gao lab "5′ and 3′ end small RNAs" approach (Xu et al. 2019; Ji et al. 2019; Chen et al. 2020; Liang et al. 2022). Ligation bias: Hafner et al. 2011, doi:10.1261/rna.2799511; Fuchs et al. 2015, doi:10.1371/journal.pone.0126049. The method also assumes gap-free tiling, which is a modelling constraint. |
| **PARE / degradome / 5′P-seq** | exact for 5′-monophosphate ends (processed mt RNAs, including 12S) | – | German et al. 2008, doi:10.1038/nbt1417. Applied to human mt by Mercer 2011 and protocolized in Rackham & Filipovska 2014, *Methods Mol Biol* 1125:263–275, doi:10.1007/978-1-62703-971-0_21 (verified). |
| **Circularized RNA-seq / cRT-PCR** | exact (needs 5′-P) | exact, including the tail sequence | Couttet et al. 1997, doi:10.1073/pnas.94.11.5628; Slomovic et al. 2005, doi:10.1128/MCB.25.15.6427-6435.2005; Kuznetsova et al. 2017. |
| **CAGE / cap-trapping** | capped TSSs only | – | Shiraki et al. 2003, doi:10.1073/pnas.2136655100. Processed rRNA 5′ ends are not capped, so this is unsuitable. |
| **Nanopore direct RNA** | **Truncated.** Workman et al. 2019 (doi:10.1038/s41592-019-0617-2): the motor enzyme "is 10–15 nt from the nanopore sensor… prevents reading the terminal 10–15 nt"; "The sharpest coverage drop-off is typically at 10 nt from the 5′ transcript end." | exact (adapter ligated to the poly(A) end); needs a tail long enough for the splint | Includes human mt RNAs, RNR1/RNR2 among them. |
| **Nano3P-seq (template switching at the 3′ end)** | not the focus | exact, regardless of poly(A) status | Begik et al. 2023; zebrafish and mouse 16S tails. |
| **TIF-seq (both ends of the same molecule)** | exact | exact | Pelechano et al. 2013, doi:10.1038/nature12121. Not yet applied to animal mt, as far as I found. |

### 2.2 Recommended pipeline for existing public paired-end RNA-seq (mt-rRNA ends)
1. **Map** to a rotated or duplicated mitogenome with a local aligner that keeps soft clips (bowtie2 `--local` or BWA-MEM).
   - Keep only properly paired reads that are uniquely mt.
   - Mask or flag NUMT-like hits. The bank vole and Kuznetsova studies filtered multireads; the Gao lab argues NUMTs are not transcribed.
   - Determine library strandedness.
2. **5′ ends:** for each position, count fragments whose RNA-sense 5′ terminus starts there (R2 start for dUTP; R1 start for "forward" or ligation kits). Compute an end ratio: starts at *p* divided by fragments covering *p*. Call the mature 5′ end at a sharp local maximum.
   - Inspect 1–3-nt soft clips at the start (non-templated additions).
   - Expect 0–1-nt agreement with primer extension (cod).
3. **3′ ends:** collect reads with ≥3 (better ≥4) non-templated 3′ A's (Marková 2015). In dUTP data these are R1 reads starting with T's. Record the last templated nt.
   - Where the genome itself has A's after the site, report the ambiguity interval, as Gissi & Pesole and Montagna did.
   - Reject internal-priming candidates, i.e. sites followed by genomic A-rich stretches in oligo-dT data.
   - If tails are too short (12S), fall back to the coverage cliff, the last nt with ≥X% of the plateau coverage. Report it as an interval with an expected underestimate of a few nt.
4. **Heterogeneity:** report the dominant site, its read support and the range. Mature 3′ ends are often ragged (*Halocynthia*: six sites in 11 bp; mouse 16S: 0–7 nt).
5. **Validate first** on genomes with experimental ends (human/mouse; cod X99772; *Ae. albopictus* NC_006817; *D. melanogaster* with Stewart & Beckenbach; *I. ricinus* JN248424; *H. roretzi* AB024528 ESTs) before applying to annelids.
6. **Coverage-based APA tools** (DaPars, Xia et al. 2014, doi:10.1038/ncomms6274; APAtrap, Ye et al. 2018, doi:10.1093/bioinformatics/bty029; TAPAS, Arefeen et al. 2018, doi:10.1093/bioinformatics/bty110; QAPA, Ha et al. 2018, doi:10.1186/s13059-018-1414-4; IsoSCM, Shenker et al. 2015, doi:10.1261/rna.046037.114) detect 3′ ends by change-points in coverage. They resolve tens of nt at best and are designed for nuclear 3′ UTRs, so they are not suitable for single-nt rRNA truth.
7. **Transcript-assembly-based end inference** (e.g. GRIT, Boley et al. 2014, doi:10.1038/nbt.2850, which needs CAGE/RAMPAGE plus poly(A)-seq) is not appropriate for compact, polycistronic mt transcripts.

**Mt-specific confounders**
- Unprocessed polycistronic precursors and antisense transcripts produce reads that span mature ends. Use start/end ratios rather than raw coverage.
- Short oligo(A) on 12S (0–5 A) is often below the soft-clip threshold.
- rRNA abundance (10–100× the mRNAs) saturates pile-ups in a useful way, but also amplifies RT and PCR artefacts.
- In oligo-dT data, internal priming inside A-rich 16S regions creates false 3′ ends.

### 2.3 Best options if new data can be generated (one sample per species is enough)
- **cRNA-seq or cRT-PCR:** both ends plus tail composition on the same molecule.
- **PARE/5′P-seq** for 5′ ends, combined with a **poly(A)-independent 3′ method** (Nano3P-seq, or 3′ adapter ligation after in vitro polyadenylation) for 3′ ends.
- **Nanopore DRS** for 3′ ends and tails, remembering that 5′ ends are truncated by about 10–15 nt.

---

## 3. Public compiled data on experimentally supported mt transcript ends
- **No dedicated cross-species database found.** MITOMAP and similar resources are human-only. MitoZoa (Lupi et al. 2010, doi:10.1016/j.mito.2010.01.004; MitoZoa 2.0, D'Onorio de Meo et al. 2012, doi:10.1093/nar/gkr1144) re-curates annotations, with an emphasis on CDS. It is not transcript-based. The original web address is caspur.it, and whether it is still online is **not verified**.
- **GenBank itself, mineable locally:**
  - Features carrying `/experiment="…"` qualifiers, e.g. *I. ricinus* JN248424.2 rrnL `/experiment="3' RACE"`.
  - Transcript features (`precursor_RNA`, "mRNA start site", promoters), e.g. *S. purpuratus* X12631.
  - Deposited RNA/cDNA end sequences, e.g. M32716 (mosquito SSU 3′-terminal RNA sequence), X53506 (*D. melanogaster* lrRNA cDNA), X01078 (mosquito LSU clone with S1-mapped ends).

  Suggestion: grep your local GenBank mirror (mtDB) for `/experiment` within rRNA/tRNA features to find more.
- **Gao lab "precise annotation" records** (sRNA-seq end small RNAs, 1-nt resolution; caveat: a gap-free tiling assumption):
  - *E. fullo* MK374364 (SRA SRP174926)
  - *D. silvarum* MN347015, with coordinates in Table 2 of Chen 2020 because the GenBank record lacks features (SRP178347)
  - *Drosophila* spp. (Liang 2022: MN611461 plus the SRA sets listed in §1.2)
  - human, chimp, macaque and mouse (Jin 2020 supplement)
- **Per-paper compilations:**
  - Montagna 2012 Table 3 and Fig. 4: rrnL 3′ ends for 11 ticks plus *Limulus* nad1.
  - Gissi & Pesole 2003 Fig. 2: *Halocynthia* poly(A) sites.
  - Coucheron 2011 Figs. 1–3: cod and saithe mRNA and rRNA ends; saithe mitogenomes FR751399–FR751401.
- **Raw 3′-end and tail datasets in non-mammal model animals** (they could be re-mined for mt reads; whether chrM reads were retained is **not verified**):
  - Nano3P-seq, zebrafish and mouse: ENA PRJEB53494 (16S tails shown).
  - PAL-seq, Subtelny et al. 2014, doi:10.1038/nature13007 (zebrafish, *Xenopus*, mammalian cells).
  - *Drosophila* PAL-seq, Eichhorn et al. 2016, doi:10.7554/eLife.16955.
  - *C. elegans* 3P-seq, Jan et al. 2011, doi:10.1038/nature09616.
  - Zebrafish 3P-seq, Ulitsky et al. 2012, doi:10.1101/gr.139733.112.
  - PolyASite 2.0, Herrmann et al. 2020, doi:10.1093/nar/gkz918 (human, mouse, worm atlas; chrM inclusion **not verified**).
- **SRA in general:** RNA-seq of almost any animal contains abundant mt-rRNA reads (12S and 16S were the most abundant transcripts in octocoral, cubozoan and louse datasets), so the §2.2 pipeline can generate ends at scale for annelids.

---

## 4. Literature on inconsistency of rRNA/tRNA boundary annotations
- **Direct quantification of rRNA end disagreement among near-identical GenBank/RefSeq records** (your 70% figure): **none found**.
- **Bernt M, et al. (2013).** MITOS. *Mol Phylogenet Evol* 69:313–319, doi:10.1016/j.ympev.2012.08.023 (verified). RefSeq annotation "is still plagued by misassignments of reading directions, erroneous gene names, and missing as well as false positive annotations in particular for the RNA genes".
- **MITOS2 tolerates rRNA overlaps by design.** In the Galaxy tools-iuc MITOS2 wrapper (MITOS 2.1.10; https://github.com/galaxyproject/tools-iuc/blob/main/tools/mitos/mitos2.xml, accessed 2026-10-26; not a peer-reviewed source):
  - `--maxrrnaovl` default **50**, help text "Allow that a tRNA/rRNA overlaps with another feature by this number of nucleotides"
  - `--maxtrnaovl` 50
  - `--finovl` 50, "Maximum number of nucleotides by which genes of different types may overlap"

  The 1–2-nt rRNA/tRNA overlaps in raw MITOS output are therefore not resolved by the software.
- **Jühling F, et al. (2012).** *Nucleic Acids Res* 40:2833–2845, doi:10.1093/nar/gkr1131 (verified).
  - MiTFi "by default regards predictions that overlap not more than 10 nt as distinct loci".
  - "more than 3700 cases of overlaps between tRNA genes in all 1876 metazoan mitogenomes… Single nucleotide overlaps are most common."
  - RefSeq "shows a high level of misannotations" of Ser and Leu tRNA identity.
- **Donath A, et al. (2019).** *Nucleic Acids Res* 47:10543–10552, doi:10.1093/nar/gkz833 (verified).
  - "RefSeq 63 contains 2085 annotations (in 1342 species…) of genes overlapping more than 10 bp."
  - Translational exceptions "should be handled with caution, unless accompanied by direct experimental evidence".
- **Al Arab M, et al. (2017).** *Mol Phylogenet Evol* 106:209–216, doi:10.1016/j.ympev.2016.09.024 (verified). Protein-coding boundary accuracy.
- **Prada CF, Boore JL (2019).** *BMC Genomics* 20:73, doi:10.1186/s12864-019-5447-1 (verified). In 304 mammals, all 68 reported gene rearrangements are annotation errors, plus "impossibly short, incorrect gene annotations".
- **Montaña-Lozano P, et al. (2022).** *Sci Rep* 12:5479, doi:10.1038/s41598-022-09512-2 (verified). 389 annotation errors among 2,831 vertebrate mitogenomes: 94.6% tRNA, 4.3% rRNA (mostly rrnS). These are gene-order and identity errors, not end precision.
- **Morrison DA (2010/2014).** "How and where to look for tRNAs in Metazoan mitochondrial genomes…" arXiv:1001.3813; bioRxiv doi:10.1101/001875 (verified).
  - "very large number of erroneous annotations… wrongly annotated locations, false complements, and inconsistent criteria for assigning the 5′ and 3′ boundaries"
  - "The extent of overlap between genes is often greatly exaggerated due to inconsistent annotations"
- **Gao lab papers** (Ji 2019; Jin 2020) state that MITOS-type annotation methods "result in gaps and overlaps". This is an assertion, not a measured rate.
- **Single-genome demonstrations that RefSeq rRNA ends are wrong:**
  - *H. roretzi* 16S, 66 nt too long (Gissi & Pesole 2003)
  - *Mytilus* l-rRNA, 47–48 nt (Kyriakou 2014)
  - *E. fullo* 12S missing 15 nt (Ji 2019)
  - *D. melanogaster* 12S, 2 nt (Liang 2022)
  - tick rrnL (Montagna 2012)
- **Cameron SL (2014).**
  - *Syst Entomol* 39:400–411, doi:10.1111/syen.12071 (verified): how to annotate insect mitogenomes, including rRNA boundaries. Full text not accessed (HTTP 403), so its exact recommendations are **not verified**.
  - *Annu Rev Entomol* 59:95–117, doi:10.1146/annurev-ento-011613-162007 (verified).
- **Boore JL (2006).** *OMICS* 10:119–126, doi:10.1089/omi.2006.10.119 (verified). Standards for organelle genome databases (content not checked).

---

## 5. The five most useful references for building an experimental truth set
1. **Montagna et al. 2012** *PLoS One* 7:e47538, doi:10.1371/journal.pone.0047538. 3′ RACE plus EST-based rrnL 3′ ends in 11 ticks; demonstrates overlap into trnL(CUN); GenBank JN248424 carries an `/experiment` tag. It is also a template for EST/RNA-seq poly(A) mapping.
2. **Elliott & Jacobs 1989** *Mol Cell Biol* 9:1069–1082, doi:10.1128/mcb.9.3.1069-1082.1989, together with **Cantatore et al. 1990** *Curr Genet* 17:235–245, doi:10.1007/bf00312615. Both ends of both rRNAs in two sea urchins: 5′ ends butt-joined, 3′ ends overlapping trnE and cox1.
3. **Coucheron et al. 2011** *Curr Genet* 57:213–222, doi:10.1007/s00294-011-0338-2 (with Bakke & Johansen 2002, doi:10.1016/s1055-7903(02)00220-8). All four cod rRNA ends at nt resolution, including oligo(A). It also shows that ordinary 454 reads reproduce primer-extension 5′ ends to within 0–1 nt.
4. **Gissi & Pesole 2003** *Genome Res* 13:2203–2212, doi:10.1101/gr.1227803. The poly(A)-start method; *Halocynthia* 16S 3′ end 66 nt upstream of the RefSeq end; 12S missing from oligo-dT data.
5. **Marková et al. 2015** *BMC Genomics* 16:870, doi:10.1186/s12864-015-2103-2. A ready-made method for calling mt 3′ ends from standard poly(A) RNA-seq with soft-clipped A's, including its stated failure for rRNAs.

   For 5′ ends and both-end methods, add Kuznetsova et al. 2017 (doi:10.1093/nar/gkx104). For insect and tick sRNA-seq-based whole-genome end annotations, add Ji et al. 2019 (doi:10.1080/15476286.2019.1591035) and Chen et al. 2020 (doi:10.1186/s12864-020-06906-2).

---

## 6. Verification log
- **Resolved via Crossref** (title and journal match):
  - all DOIs cited above, **including** 10.1128/mcb.9.3.1069-1082.1989, 10.1007/bf00312615, 10.1093/nar/12.20.7771, 10.1016/0147-619x(83)90008-2, 10.1016/0147-619x(85)90066-6, 10.1093/nar/13.11.4029, 10.1007/bf00160405, 10.1093/genetics/130.3.471, 10.1016/s1055-7903(02)00220-8, 10.1007/s00294-011-0338-2, 10.1016/j.gene.2014.01.080, 10.1371/journal.pone.0047538, 10.1101/gr.1227803, 10.1080/15476286.2016.1197481, 10.1080/15476286.2019.1591035, 10.1186/s12864-020-06906-2, 10.1016/j.mito.2022.06.006, 10.1080/15476286.2019.1709746, 10.3389/fgene.2019.00105, 10.1093/nar/gkp856, 10.1371/journal.pgen.1006028, 10.1111/j.1365-2583.2010.01053.x, 10.1186/s12867-017-0093-0, 10.1080/15476286.2016.1194161, 10.7717/peerj.1403, 10.1016/0092-8674(80)90164-6, 10.1016/s0021-9258(17)44611-4, 10.1016/0022-2836(82)90510-1, 10.1016/s0092-8674(80)80053-5, 10.1093/nar/gkx104, 10.1186/s12864-015-2103-2, 10.1038/s41592-022-01714-w, 10.7554/elife.27596, 10.1073/pnas.91.12.5368, 10.1093/nar/18.15.4592, 10.1093/nar/16.20.9880, 10.1016/0378-1119(81)90085-8, 10.1093/nargab/lqaa111, 10.3390/biology15040296, 10.1093/molbev/mss274, 10.1093/gbe/evad219, 10.1093/gbe/evt112
  - all §2 method DOIs
  - all §4 DOIs
  - the known references 10.1038/290470a0, 10.1016/j.cell.2011.06.051, 10.1016/j.gene.2009.06.006, 10.1093/molbev/mst141, 10.3389/fgene.2023.1250907
- **Additional relevant paper:** Valverde JR, Marco R, Garesse R (1994), *PNAS* 91:5368–5371, doi:10.1073/pnas.91.12.5368 (verified). It defines the conserved rRNA "termination box" heptamer TGGCAGA near the 16S 3′ end "independent of the adjacent gene", which is useful as a 3′-end feature.
- **Not verified:**
  - full-text coordinates of all paywalled or scanned papers listed at the top
  - the Mytilus 3′ inference
  - whether PolyASite or the PAL-seq/3P-seq data retain chrM reads
  - Ribo-Zero mt-rRNA depletion claims
  - MitoZoa availability
  - Cameron 2014 content
  - Liang 2022 Table 1
