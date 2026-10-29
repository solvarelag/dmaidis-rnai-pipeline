#!/usr/bin/env python3
"""Conteo del espacio de busqueda: letales iBeetle -> GeneID -> ortologo en D. maidis -> n de copias."""
import csv, re, sys
from collections import Counter, defaultdict

UMBRAL = 80  # dali11 >= UMBRAL
T = "tabla_maestra/"
ORTO = "orthofinder_pasada1/Orthologues_Dalbulus_maidis/Dalbulus_maidis__v__Tribolium_castaneum.tsv"

letales = set()
with open(T + "iBeetle_lethality.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter=";"):
        if r["dali11"] not in ("", "NA") and float(r["dali11"]) >= UMBRAL:
            letales.add(r["iB"])

ib_tc = defaultdict(set)
with open(T + "iB_TC.csv") as f:
    for line in f:
        if line.startswith("#"): continue
        ib, tc = line.strip().split(",")[:2]
        ib_tc[ib].add(tc)

tc_gid = defaultdict(set)
with open(T + "Tcas_GeneID_TC.tsv") as f:
    for gid, tc, _ in csv.reader(f, delimiter="\t"):
        if tc: tc_gid[tc].add(gid)
if "--solo-locus-tag" not in sys.argv:
    # complemento: TC de OGS3 con proteina identica en RefSeq
    # + mejor hit reciproco DIAMOND (id 90, cobertura 70/70)
    for fn in ("Tcas_TC_GeneID_secuencia.tsv", "Tcas_TC_GeneID_rbh.tsv"):
        with open(T + fn) as f:
            for tc, gid in csv.reader(f, delimiter="\t"):
                tc_gid[tc].add(gid)

# proteina Tcas (ID usado por OrthoFinder) -> GeneID
prot_gid = {}
for line in open("proteomas/Tribolium_castaneum.faa"):
    if line.startswith(">"):
        m = re.match(r">(\S+).*GeneID:(\d+)", line)
        prot_gid[m.group(1)] = m.group(2)

gid_dmai = defaultdict(set)
with open(ORTO) as f:
    next(f)
    for og, dm, tc in csv.reader(f, delimiter="\t"):
        dms = {x.strip() for x in dm.split(",")}
        for p in tc.split(","):
            gid_dmai[prot_gid[p.strip()]] |= dms

con_gid = {ib for ib in letales if any(tc_gid.get(tc) for tc in ib_tc.get(ib, ()))}
copias = {}
for ib in con_gid:
    dms = set()
    for tc in ib_tc[ib]:
        for g in tc_gid.get(tc, ()):
            dms |= gid_dmai.get(g, set())
    if dms: copias[ib] = len(dms)

dist = Counter(min(c, 3) for c in copias.values())
n = len(copias)
print(f"letales (dali11>={UMBRAL}): {len(letales)}")
print(f"  con TC: {sum(1 for ib in letales if ib_tc.get(ib))}")
print(f"  con GeneID: {len(con_gid)}")
print(f"  con ortologo en D. maidis: {n}")
for k, lab in ((1, "1 copia"), (2, "2 copias"), (3, "3+ copias")):
    print(f"    {lab}: {dist[k]} ({100*dist[k]/n:.1f}%)")
print(f"  max copias: {max(copias.values())}; iB con >1 GeneID: {sum(1 for ib in con_gid if len({g for tc in ib_tc[ib] for g in tc_gid.get(tc,())})>1)}")
