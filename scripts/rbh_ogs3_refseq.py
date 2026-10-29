#!/usr/bin/env python3
"""Mejor hit reciproco a nivel de GEN entre OGS3 (TC) y RefSeq Tcas5.2 (GeneID).
DIAMOND corrido con --id 90 --query-cover 70 --subject-cover 70 en ambos sentidos (job s4).
Isoformas colapsadas a gen tomando el maximo bitscore por par de genes. Empates en el mejor hit -> sin RBH."""
import csv, re, sys
from collections import defaultdict
p2g = dict(csv.reader(open(sys.argv[1]), delimiter="\t"))   # protein_id -> GeneID
tc = lambda s: re.search(r"(TC\d+)", s).group(1)

def mejores(path, qmap, smap):
    best = defaultdict(dict)
    for r in csv.reader(open(path), delimiter="\t"):
        q, s, bs = qmap(r[0]), smap(r[1]), float(r[7])
        if q and s: best[q][s] = max(best[q].get(s, 0), bs)
    out = {}
    for q, d in best.items():
        m = max(d.values()); tops = [s for s, v in d.items() if v == m]
        if len(tops) == 1: out[q] = tops[0]
    return out

a = mejores("score_insumos/ogs3_vs_refseq.tsv", tc, p2g.get)        # TC -> GeneID
b = mejores("score_insumos/refseq_vs_ogs3.tsv", p2g.get, tc)        # GeneID -> TC
rbh = sorted((t, g) for t, g in a.items() if b.get(g) == t)
with open("tabla_maestra/Tcas_TC_GeneID_rbh.tsv", "w") as f:
    for t, g in rbh: f.write(f"{t}\t{g}\n")
print(f"TC con hit directo: {len(a)}  GeneID con hit directo: {len(b)}  RBH: {len(rbh)}")
