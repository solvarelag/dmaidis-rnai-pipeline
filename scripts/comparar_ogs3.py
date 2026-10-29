#!/usr/bin/env python3
"""
comparar_ogs3.py — compara el ranking principal (Tribolium RefSeq + mapeo TC -> GeneID) con el analisis separado
con Tribolium OGS3 (ranking_ogs3/, IDs TC nativos). Correr despues de score_candidatos.py y score_candidatos.py ogs3.

Comparacion por gen de D. maidis (dmai_representativo), que es comun a los dos analisis.
Letales = iB con dali11 >= 80 (mismo corte que conteo_letales.py).
"""
import pandas as pd

R = pd.read_csv("ranking/ranking_completo.tsv", sep="\t", index_col=0)
O = pd.read_csv("ranking_ogs3/ranking_completo.tsv", sep="\t", index_col=0)
let = pd.read_csv("tabla_maestra/iBeetle_lethality.csv", sep=";", encoding="utf-8-sig")
ib_let = set(let.iB[let.dali11 >= 80])

def letales(U):
    return {ib for s in U.iB for ib in s.split(",")} & ib_let
lr, lo = letales(R), letales(O)
print(f"universo filtrado: RefSeq {len(R)} | OGS3 {len(O)}")
print(f"iB letales (dali11 >= 80, {len(ib_let)}) en el ranking: RefSeq {len(lr)} | OGS3 {len(lo)} | "
      f"solo OGS3 {len(lo - lr)} | solo RefSeq {len(lr - lo)}")

# copias Dmai de los letales en cada analisis (antes del filtro no esta en la tabla; se usa copias_dmai del ranking)
for n, U in (("RefSeq", R), ("OGS3", O)):
    m = U[[bool(set(s.split(",")) & ib_let) for s in U.iB]]
    print(f"  {n}: genes letales en ranking {len(m)}, copias {m.copias_dmai.value_counts().sort_index().to_dict()}")

kr = R.drop_duplicates("dmai_representativo").set_index("dmai_representativo")
ko = O.drop_duplicates("dmai_representativo").set_index("dmai_representativo")
t_r, t_o = set(kr.index[kr["rank"] <= 100]), set(ko.index[ko["rank"] <= 100])
com = kr.index.intersection(ko.index)
print(f"\ntop100 compartido (por gen Dmai): {len(t_r & t_o)}/100; genes Dmai comunes {len(com)}; "
      f"Spearman de rank en comunes {kr.loc[com, 'rank'].corr(ko.loc[com, 'rank'], method='spearman'):.3f}")
print("controles (rank RefSeq -> OGS3):")
for e in ("PSMB5", "Vha55", "Vha36-1"):
    a = R[R.entrenamiento.fillna("").str.contains(e, regex=False)]["rank"]
    b = O[O.entrenamiento.fillna("").str.contains(e, regex=False)]["rank"]
    print(f"  {e}: {a.tolist()} -> {b.tolist()}")

nuevos = ko[(ko["rank"] <= 100) & ~ko.index.isin(t_r)].copy()
nuevos["rank_refseq"] = kr["rank"].reindex(nuevos.index)
nuevos = nuevos[["rank", "rank_refseq", "dmel", "iB", "dali11", "dali22", "copias_dmai"]].sort_values("rank")
nuevos.to_csv("ranking_ogs3/top100_solo_ogs3.tsv", sep="\t", float_format="%.4g")
print(f"\nen el top100 OGS3 y no en el RefSeq: {len(nuevos)} (rank_refseq vacio = no estaba en el ranking RefSeq)")
print(nuevos.to_string(max_colwidth=20))
