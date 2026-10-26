"""Benchmark records -> taxid (from the GenBank file) -> ENA paired-end Illumina RNA-seq runs of that taxon."""
import csv, glob, os, re, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
sets = {"large": ("/media/deyuan/18T_B/yzz/MitoM/benchmark/sets/mito_bench_large_2026-10.tsv", "/media/deyuan/18T_B/yzz/mitom_work/bench_large"),
        "annelida": ("/media/deyuan/18T_B/yzz/MitoM/benchmark/sets/mito_bench_annelida_2026-10.tsv", "/media/deyuan/18T_B/yzz/mitom_work/bench_annelida")}
recs = []
for sname, (tsv, root) in sets.items():
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        gb = os.path.join(root, "gb", r["accession"] + ".gb")
        tx = ""
        with open(gb) as fh:
            for line in fh:
                m = re.search(r'/db_xref="taxon:(\d+)"', line)
                if m:
                    tx = m.group(1); break
        recs.append((sname, r["accession"], r["phylum"], r["species"], tx))
def q(rec):
    sname, acc, ph, sp, tx = rec
    if not tx: return rec + ("", 0)
    query = f'tax_eq({tx}) AND library_source="TRANSCRIPTOMIC" AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"'
    url = ("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&limit=50&fields="
           "run_accession,read_count,base_count,library_strategy,library_selection,nominal_length,fastq_ftp,first_public"
           "&query=" + urllib.parse.quote(query))
    for k in range(4):
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()
            lines = [l.split("\t") for l in txt.strip().split("\n")[1:] if l.strip()]
            return rec + (lines, len(lines))
        except Exception as e:
            time.sleep(2 + 3 * k)
    return rec + ("ERR", -1)
with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(q, recs))
with open("ena_runs.tsv", "w") as out:
    out.write("set\taccession\tphylum\tspecies\ttaxid\trun\tread_count\tbase_count\tstrategy\tselection\tnominal\tfastq_ftp\tfirst_public\n")
    n_sp = 0
    for sname, acc, ph, sp, tx, lines, n in res:
        if n > 0:
            n_sp += 1
            for l in lines:
                out.write("\t".join([sname, acc, ph, sp, tx] + l) + "\n")
print("records", len(res), "with RNA-seq", n_sp, "errors", sum(1 for x in res if x[-1] == -1))
