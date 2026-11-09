# RNAcode as a coding test for mitochondrial ORFs (2026-11-10)

RNAcode 0.3.1 (bioconda; Washietl et al. 2011) scores a nucleotide alignment for protein-coding selection. It has the
standard genetic code only (hard-coded, src/code.c), so codons whose mitochondrial meaning differs were recoded in
frame to a standard synonym (TGA -> TGG, AGR -> Ser codon, ATA -> ATG ...). `compare.py`: each candidate ORF with
homologs (the true ORFs of ../make_truth.py, record left out; all candidates of the 362 intron-free records) gets a
protein-guided codon alignment (MAFFT) with up to 12 homologs and RNAcode's p for frame +1, next to mito_orf's
p_open and codon-pattern p (`*_candidates.tsv`).

| 规则 / rule | 真 ORF / true (258) | 假候选 / intron-free candidates (683, minus fragments of the 13) |
|---|---|---|
| RNAcode p <= 1e-3 | 239 | 45 |
| RNAcode <= 1e-4, p_open <= 0.05, codon pattern <= 0.05 | 224 | 6 |
| p_open <= 0.05 and codon pattern <= 0.05 (no RNAcode; adopted) | 231 | 8 |
| p_open <= 0.01 and codon pattern <= 0.01 | 204 | 4 |

RNAcode alone calls many control-region stretches; with the two mito_orf tests it adds almost nothing, so it is not
used (one dependency less). Of the 8 intron-free "false" calls, 4 are probably real (Aurelia polB and its neighbour,
a Pseudanodonta F-ORF, a glass-sponge ORF) and 3 are in bivalves / an octocoral.
