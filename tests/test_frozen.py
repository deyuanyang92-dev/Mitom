"""benchmark/frozen: leak audit / test-set selection (frozen_set.py) and the pre-registered metrics (frozen_stats.py).

Run: python3 -m pytest -q tests/test_frozen.py
"""
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "benchmark" / "frozen"))
import frozen_set as fs  # noqa: E402
import frozen_stats as st  # noqa: E402


def rec(acc, genus, family, phylum="Annelida", src="", length=15000, cox1=True, n_cds=13):
    return {"accession": acc, "source_accession": src, "organism": f"{genus} sp.", "genus": genus, "family": family,
            "phylum": phylum, "lineage": f"Metazoa;{phylum};X;{family};{genus}", "length": length,
            "n_cds": n_cds, "cox1": cox1}


def test_lineage_taxa():
    assert fs.lineage_taxa("Metazoa;Annelida;Polychaeta;Nereididae;Platynereis".split(";")) == ("Platynereis", "Nereididae")
    assert fs.lineage_taxa([]) == ("", "")


def test_repository_accessions_count_as_seen():
    accs = fs.seen_accessions()
    assert "NC_000931" in accs                    # tests/data, intron gold, tRNA splits ...
    assert not any(v.startswith("benchmark/frozen") for v in accs.values())


def test_audit_and_select():
    pool = [rec("NC_1.1", "Seen", "Fa"), rec("NC_2.1", "Newa", "Fb", src="MK1"), rec("NC_3.1", "Newb", "Fb"),
            rec("NC_4.1", "Newc", "Fb"), rec("NC_5.1", "Newd", "Fc", cox1=False), rec("NC_6.1", "Newe", "Seenfam")]
    accs = {"MK1": "benchmark/x.tsv"}                                # its source (GenBank) record was evaluated
    taxa = ({"Seen sp."}, {"Seen"}, {"Fa", "Seenfam"})
    aud = fs.audit(pool, accs, taxa, {"Fb"})
    by = {r["accession"]: r for r in aud}
    assert by["NC_2.1"]["acc_seen"] == "benchmark/x.tsv" and by["NC_1.1"]["genus_seen"] == "yes"
    assert by["NC_3.1"]["clade_cm"] == "Fb" and by["NC_3.1"]["clade_cm_rank"] == "2"
    got = [r["accession"] for r in fs.select(aud, {r["accession"]: r for r in pool}, "genus", 40, 2, 1, 1)]
    assert sorted(got) == ["NC_3.1", "NC_4.1", "NC_6.1"]             # NC_5 no cox1; family Fb capped at 2
    got = [r["accession"] for r in fs.select(aud, {r["accession"]: r for r in pool}, "family", 40, 1, 1, 1)]
    assert len(got) == 1 and got[0] in ("NC_3.1", "NC_4.1")           # family level: Seenfam out, Fb capped at 1


def test_stats_paired_difference(tmp_path):
    feats = [{"tool": t, "accession": a, "kind": "PCG", "class": c}
             for t, a, c in [("m", "A", "exact"), ("m", "A", "exact"), ("m", "B", "exact"), ("m", "B", "le3"),
                             ("r", "A", "exact"), ("r", "A", "gt3"), ("r", "B", "missed"), ("r", "B", "le3")]]
    c = st.counts(feats, [{"tool": "r", "accession": "B"}], ["m", "r"])
    assert st.value([c["m"]["A"], c["m"]["B"]], "pcg_exact") == 0.75
    assert st.value([c["r"]["A"], c["r"]["B"]], "pcg_found") == 0.75
    assert st.value([c["r"]["A"], c["r"]["B"]], "err_free") == 0.5
    for name, rows in (("features", feats), ("errors", [{"tool": "r", "accession": "B", "error": "x"}])):
        with (tmp_path / f"{name}.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    (tmp_path / "set.tsv").write_text("accession\tfamily\nA\tF1\nB\tF2\n")
    out = tmp_path / "stats.tsv"
    assert st.main(["--eval", str(tmp_path), "--set", str(tmp_path / "set.tsv"), "--tools", "m,r", "--ref", "r",
                    "--boot", "200", "-o", str(out)]) == 0
    rows = {(r["metric"], r["tool"]): r for r in csv.DictReader(out.open(), delimiter="\t")}
    assert rows[("pcg_exact", "m - r")]["value_pct"] == "+50.0" and rows[("pcg_exact", "m")]["families"] == "2"
