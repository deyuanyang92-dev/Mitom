"""rRNA, tRNA and name rules of mito_curate.Curator: R1 rRNA ends, T1 weak tRNAs, C1 tRNAs inside protein genes, N1 NCBI names.

Code moved verbatim from mito_curate.py; the methods use the Curator helpers (self.log, self.tx, ...)."""
from __future__ import annotations

import re
from typing import Dict, List

from Bio.SeqFeature import AfterPosition, BeforePosition, FeatureLocation, SeqFeature

from mito_curate_base import (
    INTERNAL_NOTE,
    MAX_INSIDE_TRIM,
    NO_START_NOTE,
    NO_STOP_NOTE,
    WEAK_TRNA_E,
    gene_feature,
    simple,
    span,
    spans,
)
from mito_names import AA1_TO_3, AA3, PCG_NAMES, PCG_SYN, RRNA_PRODUCT, trna_symbol
from mito_names import PRODUCT as RRNA_PRODUCT_NAMES
from mito_names import RRNA_ALIAS as RRNA_NAMES
from mito_names import name_key as key


class RnaRules:
    """Mixin of mito_curate.Curator."""

    def rrna(self, f: SeqFeature) -> None:
        # an rRNA across the origin is one interval s..e with e > L (MitoM's frame starts at cox1, so the echinoderm
        # rrnL before cox1 crosses it: NC_013881.1 ran 19 nt into cox1 and left 46 nt to nad2 unfilled)
        origin = self.origin_span(f)
        if not simple(f) and origin is None:
            return
        g = gene_feature(self.rec, f)
        s, e = old = origin or span(f)
        nf = [(sp, x) for x in self.others(f, types=("tRNA", "CDS", "rRNA")) for sp in spans(x)]
        if self.rec.annotations.get("topology") == "circular":        # neighbours across the origin
            nf += [((a + k, b + k), x) for (a, b), x in list(nf) for k in (-self.L, self.L)]
        neigh = [sp for sp, _ in nf]
        who = {sp: x.qualifiers.get("gene", ["?"])[0] for sp, x in nf}
        noted = self.noted            # R3 'inside' rows once per rRNA pass pair (rrna() runs twice in run())
        changed = True
        while changed:                                # R1: no overlap with any other gene (repeat until stable)
            changed = False
            for ns, ne in sorted(neigh):
                if ne < s or ns > e:
                    continue
                if ns <= s and ne >= e:               # rRNA inside another gene: left to R3 / manual check
                    continue
                if ns <= s:                           # overlaps the 5' (left) end
                    s = ne + 1
                elif ne >= e:                         # overlaps the 3' (right) end
                    e = ns - 1
                else:                                 # gene inside the rRNA: keep the longer side, but only if what is
                    keep_left = ns - s >= e - ne      # lost is <= 10 % or <= MAX_INSIDE_TRIM nt (a tRNA-sized overhang
                    lost = e - ns + 1 if keep_left else ne - s + 1   # of the core; a spurious tRNA call deep inside
                    if lost * 10 > e - s + 1 and lost > MAX_INSIDE_TRIM:   # must not cut a core rRNA)
                        if (ns, ne) not in noted:
                            noted.add((ns, ne))
                            self.log(f, "R3", old, old, f"gene {who[(ns, ne)]} inside rRNA, not trimmed (would lose "
                                     f"{lost} nt) - check manually (spurious tRNA call?)")
                        continue
                    if keep_left:
                        e = ns - 1
                    else:
                        s = ne + 1
                changed = True
        if self.a.rrna_bounds == "fill":
            lefts = [ne for ns, ne in neigh if ne < s]          # fill only up to a neighbouring gene, never to the
            rights = [ns for ns, ne in neigh if ns > e]         # end of the sequence (2026-11-20 review)
            if lefts and 0 < s - max(lefts) - 1 <= self.a.max_fill:
                s = max(lefts) + 1
            if rights and 0 < min(rights) - e - 1 <= self.a.max_fill:
                e = min(rights) - 1
        cds_ov = max([min(e, ne) - max(s, ns) + 1 for x in self.others(f, types=("CDS",)) for ns, ne in spans(x)] or [0])
        if cds_ov > 50:
            self.log(f, "R3", old, (s, e), f"rRNA overlaps a CDS by {cds_ov} nt - boundary or assembly problem?")
        if (s, e) != old and e - s + 1 >= 200:
            # an end at the end of a linear sequence that R1 leaves in place keeps its < / >: the rRNA runs on beyond
            # the sequenced part (DQ517436.1, a partial genome starting inside rrnL)
            linear = self.rec.annotations.get("topology") != "circular"
            keep = linear and (isinstance(f.location.start, BeforePosition) and s == old[0] == 1,
                               isinstance(f.location.end, AfterPosition) and e == old[1] == self.L)
            self.set_circular_span(f, s, e)
            if keep and any(keep):
                f.location = FeatureLocation(BeforePosition(s - 1) if keep[0] else f.location.start,
                                             AfterPosition(e) if keep[1] else f.location.end, f.location.strand)
            if g is not None and (simple(g) or origin):
                g.location = f.location
            self.log(f, "R1", old, (s, e), "trimmed to flanking genes" + (" and filled to them" if self.a.rrna_bounds == "fill" else ""))

    def weak_trnas(self) -> None:
        """T1: a weak tRNA call (MiTFi E > WEAK_TRNA_E, or a second copy: the core's note) lying wholly inside a
        same-strand protein gene is dropped - tRNAs punctuate the transcripts, a whole one inside a gene on its strand is
        a false call or a wrong gene end (table2asn SEQ_FEAT.BadCDScomponentOverlapTRNA; 2026-10 runs: 7 of 9 such calls
        absent from NCBI). Only a second copy of a family, or a call in the middle third of the gene: one in its first or
        last third is more often a gene end placed across it (NC_060482.1 / NC_060708.1: nad6 extended to the stop
        inside trnD, NCBI ends nad6 right before trnD) and P1 / P3 move that end instead. Strong calls are never dropped.
        Run before the protein gene rules, else P3 ends a gene with the T / TA it finds right before such a call."""
        for t in [x for x in self.rec.features if x.type == "tRNA" and id(x) in self.weak]:
            host = next((c for c in self.rec.features if c.type == "CDS" and c.location.strand == t.location.strand
                         and all(any(cs <= ts and te <= ce for cs, ce in spans(c)) for ts, te in spans(t))), None)
            if host is None:
                continue
            a, b = self.tx(host) if simple(host) else (0, 0)
            ts, te = self.tx(t) if simple(t) else (0, 0)
            n = b - a + 1
            if simple(host) and simple(t) and id(t) not in self.copies and (ts - a < n / 3 or b - te < n / 3):
                continue                               # at a gene end: P1 / P3 move the end
            g = gene_feature(self.rec, t)
            self.rec.features = [x for x in self.rec.features if x is not t and x is not g]
            self.log(t, "T1", span(t), span(t), f"weak tRNA call inside {host.qualifiers.get('gene', ['?'])[0]} "
                     f"{span(host)[0]}..{span(host)[1]} (same strand) removed")

    def trna_conflicts(self) -> None:
        """C1: a tRNA still lying wholly inside a same-strand protein gene after P1 / P3 (table2asn
        BadCDScomponentOverlapTRNA): a gene without start or stop codon (P5 / P6) or not translatable (P10) is a false
        call and is removed (NC_061193.1 flatworm 'atp8' around tRNA-Glu, NC_013814.1 chaetognath 'atp6' around
        tRNA-Phe: neither phylum has the gene); a complete gene is kept and the conflict reported (P4)."""
        for c in [x for x in self.rec.features if x.type == "CDS"]:
            # inside one PART of the CDS: the hull of a join across the origin is 1..L and 'contains' every tRNA
            # (2026-10-31 review: partial cox1 / nad5 across the origin were removed by C1)
            inside = [t for t in self.rec.features if t.type == "tRNA" and t.location.strand == c.location.strand
                      and simple(t) and any(a <= span(t)[0] and span(t)[1] <= b for a, b in spans(c))]
            if not inside:
                continue
            notes = c.qualifiers.get("note", [])
            if "pseudo" in c.qualifiers or NO_START_NOTE in notes or NO_STOP_NOTE in notes:
                g = gene_feature(self.rec, c)
                self.rec.features = [x for x in self.rec.features if x is not c and x is not g]
                self.log(c, "C1", span(c), span(c), f"incomplete or untranslatable gene around "
                         f"{inside[0].qualifiers.get('gene', ['?'])[0]} {span(inside[0])[0]}..{span(inside[0])[1]} removed")
            else:
                self.log(c, "P4", span(c), span(c), f"contains {inside[0].qualifiers.get('gene', ['?'])[0]} "
                         f"{span(inside[0])[0]}..{span(inside[0])[1]} on its strand - check manually")

    # ── N1: NCBI-style names (every feature, also joins / origin-spanning ones) ──
    def names(self, f: SeqFeature) -> None:
        g = gene_feature(self.rec, f)
        if f.type == "CDS" and simple(f):                  # P11 reads the core's alignment (its note is dropped below)
            m = re.search(r"residues (\d+)-(\d+)/(\d+)", " ".join(f.qualifiers.get("note", [])))
            if m:
                self.aln[id(f)] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
            self.core_start[id(f)] = self.tx(f)[0]
        for n in f.qualifiers.get("note", []) if f.type == "tRNA" else []:
            m = re.match(r"MiTFi \S+ bits [\d.]+ E ([0-9.eE+-]+)( copy)?", n)
            if m and (float(m.group(1)) > WEAK_TRNA_E or m.group(2)):
                self.weak.add(id(f))                   # T1 (the note itself is dropped below)
                if m.group(2):
                    self.copies.add(id(f))
        notes = [n for n in f.qualifiers.get("note", []) if not INTERNAL_NOTE.match(n)]
        if notes:
            f.qualifiers["note"] = notes
        else:
            f.qualifiers.pop("note", None)
        cur = f.qualifiers.get("gene", [""])[0]
        q: Dict[str, List[str]] = {}
        if f.type == "CDS":
            name = PCG_SYN.get(key(cur)) or PCG_SYN.get(key(f.qualifiers.get("product", [""])[0]))
            if name:
                q = {"gene": [name], "product": [PCG_NAMES[name][0]]}
        elif f.type == "rRNA":
            prod = f.qualifiers.get("product", [""])[0]
            name = RRNA_NAMES.get(key(cur)) or RRNA_PRODUCT_NAMES.get(prod.strip().lower()) or RRNA_NAMES.get(key(prod)) \
                or {"12sribosomalrna": "rrnS", "16sribosomalrna": "rrnL"}.get(key(prod))   # (2026-10-31 review: l-rRNA,
            #                                                                 small subunit ribosomal RNA were not read)
            if name:
                q = {"gene": [name], "product": [RRNA_PRODUCT[name]]}
        elif f.type == "tRNA":
            # one of three readers of a tRNA's amino acid / anticodon, kept apart on purpose (see mito_names.trna_name)
            txt = " ".join(f.qualifiers.get("product", []) + [cur])
            m = re.search(r"tRNA-([A-Z][a-z]{2})", txt)
            aa = AA3.get(m.group(1), "") if m else ""
            if not aa and re.match(r"trn([A-Z])", cur):
                aa = cur[3]
            ac = ""
            m = re.search(r"\(([acgtuACGTU]{3})\)", cur) or \
                re.search(r"anticodon:\s*([acgtuACGTU]{3})", " ".join(f.qualifiers.get("note", []))) or \
                re.search(r"seq:([acgtuACGTU]{3})", " ".join(f.qualifiers.get("anticodon", [])))
            if m:
                ac = m.group(1).lower().replace("u", "t")
            if aa in AA1_TO_3:
                name = trna_symbol(aa, ac)
                if not ac and re.fullmatch(r"trn[LS][12]", cur):
                    name = cur                         # anticodon unknown: the core numbered the gene from the MiTFi model
                notes = [n for n in f.qualifiers.get("note", []) if not n.startswith("anticodon:")]
                q = {"gene": [name], "product": [f"tRNA-{AA1_TO_3[aa]}"],
                     "note": ([f"anticodon: {ac.upper().replace('T', 'U')}"] if ac else []) + notes}
        if not q:
            if f.type in ("CDS", "rRNA", "tRNA"):
                self.log(f, "N1", span(f), span(f), f"gene name '{cur}' not recognised - kept")
            return
        if g is None and not cur and f.type != "CDS":
            q.pop("gene", None)          # no gene feature and no /gene: none added (table2asn GeneXrefWithoutGene x 22
            #                              on the tRNAs / rRNAs of NCBI NC_001636.1, 2026-10-31 review)
        for k, v in q.items():
            if v:
                f.qualifiers[k] = v
            else:
                f.qualifiers.pop(k, None)
        if g is not None:
            g.qualifiers["gene"] = q["gene"]
