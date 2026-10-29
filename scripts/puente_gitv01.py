#!/usr/bin/env python3
"""Tarea G: puente Helixer <-> GITV01 por mejor hit reciproco (DIAMOND blastp, --id 80 --query-cover 60 --subject-cover 60).
Proteinas GITV01 = TransDecoder (--single_best_only) sobre GITV01; se reportan como ID de transcripto (sin '.pN').
Empate en el mejor bitscore -> sin mejor hit. Columna extra: otros transcriptos GITV01 cuyo mejor hit es ese gen
(isoformas/alelos redundantes del ensamblado), util para buscar en la literatura."""
import csv
from collections import defaultdict
tr = lambda s: s.rsplit(".p", 1)[0]
def mejores(path, fq, fs):
    best = defaultdict(dict)
    for r in csv.reader(open(path), delimiter="\t"):
        q, s, b = fq(r[0]), fs(r[1]), float(r[7])
        if b > best[q].get(s, (0,))[0]: best[q][s] = (b, float(r[2]))
    out = {}
    for q, d in best.items():
        m = max(v[0] for v in d.values()); tops = [s for s, v in d.items() if v[0] == m]
        if len(tops) == 1: out[q] = (tops[0], d[tops[0]][1])
    return out
h2g = mejores("score_insumos/h2g.tsv", str, tr)
g2h = mejores("score_insumos/g2h.tsv", tr, str)
otros = defaultdict(list)
for g, (h, _) in g2h.items(): otros[h].append(g)
with open("mapeo_helixer_GITV01.tsv", "w") as f:
    f.write("helixer\tid_GITV01\tpident\totros_GITV01\n")
    n = 0
    for h, (g, pid) in sorted(h2g.items()):
        if g2h.get(g, (None,))[0] == h:
            n += 1
            f.write(f"{h}\t{g}\t{pid}\t{','.join(sorted(x for x in otros[h] if x != g))}\n")
print(f"Helixer con hit directo: {len(h2g)} | GITV01 con hit directo: {len(g2h)} | RBH: {n} de 22199 genes")
