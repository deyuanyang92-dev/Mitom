# Nuclear marker reference library (mito_nuc.py)

Built 2026-10-23 with `python3 mito_nuc.py build --gb rdna.gb histone.gb -o models/nuc --threads 32` from the NCBI records
that `mitom nuc fetch` downloads (queries in `mito_nuc.QUERIES`, fetched 2026-10-23):

- rDNA: `txid33208[Organism:exp] AND 18S[Title] AND 5.8S[Title] AND 28S[Title] AND "complete sequence"[Title]
  NOT partial[Title]` -> 667 records; 548 have at most one rRNA feature per gene (18S / 5.8S / 28S, named by
  product, note or gene) on one strand; the first 3 per species are kept (441 records; one species has 97
  records in NCBI). Records that annotate the unit only as one misc_RNA "contains 18S ...", mostly Cnidaria,
  have no rRNA features and are skipped.
- histones: `txid33208[Organism:exp] AND ("histone H3"[Title] OR ... H4, H2A, H2B, H1) AND "complete cds"[Title]
  AND 300:20000[SLEN]` -> 631 records; 677 CDS named histone H1/H2A/H2B/H3/H4, not partial, codon_start 1, with a
  translation of <= 400 aa, and not a histone-modifying or -binding protein (product without factor,
  transferase, acetyl, chaperone, binding, kinase, demethylase, ligase, DOT1, lysine, regulator, interact,
  domain, methyl).

Files:

| file | content |
|---|---|
| rdna_ref.fa.gz | each rDNA record, reverse-complemented when its rRNAs are on the minus strand |
| rdna_ref.tsv | ref, organism, genus, taxonomy, study (first publication title), gene, start, end (on the oriented sequence), start_ok / end_ok (0 = that end is marked partial or failed the QC below and is not used as a boundary) |
| rdna_qc.tsv | reference boundaries set not ok by the QC: annotated position, position most other studies map to, number of studies agreeing with each |
| hist_ref.faa | histone proteins, id = accession_CDS-start |
| hist_ref.tsv | ref, accession, organism, genus, taxonomy, gene, length_aa |

Taxon mix (rDNA records): Arthropoda most, then Nematoda, Chordata, Lophotrochozoa (annelids: Hesionidae,
Polynoidae, Sigalionidae and others; no Nereididae), a few Porifera and Cnidaria. Boundaries are as the
submitters annotated them; conventions differ (e.g. 5.8S 5' end AAACCATT... vs ACAACTCT...), which caps the exact
agreement with any one record.

Curated (`--curated curated.gb`, query `mito_nuc.QUERIES["curated"]`): 28 RefSeq NR_ rRNA whose ends agree between
species (human N1-N5, rat 18S/28S, macaque, rabbit, cow, hamster 5.8S, Drosophila 18S, C. elegans 18S/5.8S). Left out
because their ends disagree with the others: mouse NR_003278-80 (18S shifted 1 nt, 28S/5.8S 3' ends), rat 5.8S
NR_046238 and human N4 5.8S (1 nt short), hamster 18S (+4 nt), pig 18S (2302 nt), horse (truncated), macaque
NR_146159, Drosophila 28S (hidden break), Saccoglossus. A curated reference aligned over a boundary decides it.

QC (`mito_nuc.qc_rdna`): each reference is aligned (megablast) with the references of other studies, their
boundaries are mapped onto it, and a boundary is not used when the position most studies map to is > 15 nt away and
supported by >= 3 studies and > 2x the studies agreeing with the annotation. Records of one study vote once, so a
study that annotated all its records the same wrong way (e.g. a 28S 3' end ~530 nt before the conserved ...GATTCGT
end, ticks KY4575xx; a firefly study whose 28S ends scatter by up to 2 kb) cannot outvote the others.

Own references (e.g. Geneious exports of Perinereis rDNA): annotate 18S / 5.8S / 28S as rRNA features with
those names and histones as CDS with product "histone H3" etc., then
`mitom nuc build --gb rdna.gb histone.gb my_refs.gb -o my_lib` and `mitom nuc annotate ... --lib my_lib`.

Rfam models (`rfam/`, downloaded 2026-10-25 from https://rfam.org/family/<acc>/cm, Rfam data are CC0): RF01960
SSU_rRNA_eukarya (CLEN 1831), RF00002 5_8S_rRNA (154), RF02543 LSU_rRNA_eukarya (3401), Infernal 1.1.4. `mitom nuc
annotate` places every rRNA end at the model's first / last column plus the nt the curated RefSeq rRNAs have outside
it (18S 0 / 0, 5.8S 1 / 3, 28S 2 / 1; `mito_nuc.CM_OFFSET`); `--no_cm` turns this off.
