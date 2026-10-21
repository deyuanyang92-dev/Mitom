"""mito_compare.py: original annotation vs MitoM, gene by gene."""
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_compare as mc  # noqa: E402


def _rec(seq, feats):
    r = SeqRecord(seq, id="X1.1", name="X1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular")
    r.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source")] + feats
    return r


def test_summary_counts_each_gene_once():
    n = Counter(records=1, same=2, ends_differ=3)
    assert mc.summary_text(n).startswith("1 records, 5 genes: same 2, ends_differ 3")


def test_record_without_original_annotation_adds_no_rows(tmp_path):
    seq = Seq("ATG" + "GCC" * 50 + "TAA" + "GCC" * 30)
    cds = SeqFeature(FeatureLocation(0, 156, 1), type="CDS",
                     qualifiers={"gene": ["COX1"], "product": ["cytochrome c oxidase subunit 1"], "transl_table": ["5"]})
    (tmp_path / "o").mkdir()
    (tmp_path / "m").mkdir()
    SeqIO.write(_rec(seq, []), str(tmp_path / "o" / "X1.1.gb"), "genbank")
    SeqIO.write(_rec(seq, [cds]), str(tmp_path / "m" / "X1.1.gb"), "genbank")
    n = mc.compare_dirs(tmp_path / "o", tmp_path / "m", tmp_path / "out.tsv")
    assert n["records"] == 1 and n["records_without_original_annotation"] == 1 and n["only_mitom"] == 0
    assert (tmp_path / "out.tsv").read_text().count("\n") == 1          # header only
    assert "1 records had no original annotation" in mc.summary_text(n)


def test_same_gene_is_same_after_rotation(tmp_path):
    seq = Seq("ATG" + "GCC" * 50 + "TAA" + "GCC" * 30)
    cds = SeqFeature(FeatureLocation(0, 156, 1), type="CDS",
                     qualifiers={"gene": ["COX1"], "product": ["cytochrome c oxidase subunit 1"], "transl_table": ["5"]})
    rot = Seq(str(seq[30:] + seq[:30]))
    cds_rot = SeqFeature(FeatureLocation(216, 246, 1), type="CDS", qualifiers=dict(cds.qualifiers))
    (tmp_path / "o").mkdir()
    (tmp_path / "m").mkdir()
    SeqIO.write(_rec(seq, [cds]), str(tmp_path / "o" / "X1.1.gb"), "genbank")
    rec_m = _rec(rot, [])
    rec_m.features.append(SeqFeature(FeatureLocation(216, 246, 1) + FeatureLocation(0, 126, 1), type="CDS",
                                     qualifiers=dict(cds.qualifiers)))
    SeqIO.write(rec_m, str(tmp_path / "m" / "X1.1.gb"), "genbank")
    n = mc.compare_dirs(tmp_path / "o", tmp_path / "m", tmp_path / "out.tsv")
    assert n["same"] == 1 and n["records"] == 1, dict(n)
