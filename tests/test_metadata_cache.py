"""The metadata source of a multi-record input is parsed once per run, not once per sample."""
import mitom_metadata as mm
from mitom_output import finalize_genbank_output


def _record(i: int, with_meta: bool = True) -> str:
    acc = f"MN{i:06d}"
    seq = "ACGT" * 10
    return (
        f"LOCUS       {acc}                  40 bp    DNA     circular INV 01-JAN-2020\n"
        f"DEFINITION  Testus species {i} mitochondrion, complete genome.\n"
        f"ACCESSION   {acc}\nVERSION     {acc}.1\n"
        "SOURCE      mitochondrion Testus species\n  ORGANISM  Testus species\n"
        "            Eukaryota; Metazoa.\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..40\n"
        '                     /organism="Testus species"\n'
        '                     /mol_type="genomic DNA"\n'
        "ORIGIN\n"
        f"        1 {seq[:10]} {seq[10:20]} {seq[20:30]} {seq[30:40]}\n//\n")


def test_source_parsed_once(tmp_path, monkeypatch):
    n = 12
    src = tmp_path / "all.gb"
    src.write_text("".join(_record(i) for i in range(n)))
    mm.clear_source_caches()
    calls = {"index": 0, "org": 0, "text": 0}
    orig_index, orig_org, orig_read = mm.build_metadata_index, mm._build_organism_index, mm._source_text_cached.__wrapped__
    monkeypatch.setattr(mm, "build_metadata_index",
                        lambda t: (calls.__setitem__("index", calls["index"] + 1), orig_index(t))[1])
    monkeypatch.setattr(mm, "_build_organism_index",
                        lambda f: (calls.__setitem__("org", calls["org"] + 1), orig_org(f))[1])
    # cache wrappers were built with the original functions; rebuild them so the counters are seen
    import functools
    monkeypatch.setattr(mm, "_source_text_cached", functools.lru_cache(maxsize=2)(
        lambda k: (calls.__setitem__("text", calls["text"] + 1), orig_read(k))[1]))
    monkeypatch.setattr(mm, "_metadata_index_cached", functools.lru_cache(maxsize=2)(
        lambda k: mm.build_metadata_index(mm._source_text_cached(k))))
    monkeypatch.setattr(mm, "_meaningful_cached", functools.lru_cache(maxsize=2)(
        lambda k: any(mm._meta_is_meaningful(mm._split_meta_body(r)[0])
                      for r in mm._split_gb_records(mm._source_text_cached(k)))))
    monkeypatch.setattr(mm, "_organism_index_cached", functools.lru_cache(maxsize=2)(
        lambda ks: mm._build_organism_index([__import__("pathlib").Path(k[0]) for k in ks])))
    import mitom_output
    monkeypatch.setattr(mitom_output, "gb_file_has_meaningful_metadata", mm.gb_file_has_meaningful_metadata)
    monkeypatch.setattr(mitom_output, "build_organism_index", mm.build_organism_index)

    for i in range(n):
        body = tmp_path / f"body{i}.gb"
        body.write_text(_record(i))
        finalize_genbank_output(
            sample_name=f"S{i}", body_src_gb=body, metadata_src_gb=src, dest_gb=tmp_path / f"out{i}.gb",
            which_gene_first="COX1", force_forward=False, preserve_metadata=True)
    assert calls == {"index": 1, "org": 1, "text": 1}, calls


def test_cache_follows_file_change(tmp_path):
    src = tmp_path / "a.gb"
    src.write_text(_record(1))
    mm.clear_source_caches()
    assert mm.gb_file_has_meaningful_metadata(src)
    src.write_text(_record(1) + _record(2))
    assert len(mm.build_organism_index([src])) >= 2
    idx = mm.build_organism_index([src])
    idx.clear()                                   # a caller's copy never damages the cache
    assert mm.build_organism_index([src])
