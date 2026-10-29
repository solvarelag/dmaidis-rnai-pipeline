#!/usr/bin/env python3
"""
rnaseq_integracion.py — brief_rnaseq: agrega al ranking evidencia de RNA-seq y de GITV01.

Columnas nuevas (por gen Tcas del ranking; se evalua el ortologo Dmai representativo, y para paralogo_resuelto ambas copias):
  rnaseq_cobertura   fraccion de bases exonicas (Helixer, incl. UTR) con >= 1 lectura (bedtools coverage -split, BAM de
                     19 corridas; -split para que el intron de una lectura empalmada no cuente como cobertura)
  rnaseq_lecturas    fragmentos con algun bloque en un exon del gen (read1 primario, MAPQ >= 10), suma de las 19 corridas
                     (rnaseq_insumos/conteos/, sbatch r4)
  rnaseq_soporte     cobertura >= 0.5 y lecturas >= 10 (ver UMBRALES)
  paralogo_resuelto  solo genes con 2 copias: ambas_con_lecturas / solo_una / ninguna
  estabilidad_estadios    CV entre ninfa, hembra y macho (PRJNA318849; = CV_estadios, que ya esta en el subscore de expresion)
  estabilidad_condiciones CV entre las 4 condiciones de PRJNA631706 (replicas tecnicas colapsadas); senal debil, informativa,
                          NO entra al score. CV menor = mas estable. Ambas de expresion/ (scripts/expresion_matrices.py)
  log2FC_MRFV_4h/_7d      descriptivos (sin test ni umbral), PRJNA631706 colapsado
  id_GITV01          mejor hit reciproco proteina Helixer <-> proteina GITV01 (tarea G)
  soporte_GITV01     RBH o, si no, un transcripto GITV01 alineado (minimap2) sobre exones del gen
  confianza_gen      alta / media / media_revisar / baja (tabla del brief)
  respuesta_MRFV     up_4h / down_4h / up_7d / down_7d (Tablas 3 y 4 de Xu et al. 2021), via GITV01 -> gen
  loci_genoma        Tarea F: loci con alineamiento del mRNA (identidad >= 0.90, cobertura de query >= 0.5)

UMBRALES de rnaseq_soporte: con ~390 M de pares, un gen expresado aunque sea bajo en adultos acumula decenas de lecturas;
>= 10 lecturas evita contar lecturas sueltas (ADN genomico contaminante, multimapeo en TEs), y >= 50% de las bases
exonicas cubiertas exige que la transcripcion recorra el modelo y no un fragmento.
"""
import csv, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

RS = Path("rnaseq_insumos")
MIN_COB, MIN_LECT = 0.5, 10
MIN_ID_LOCUS, MIN_QCOV_LOCUS = 0.90, 0.5

U = pd.read_csv("ranking/ranking_completo.tsv", sep="\t", index_col=0)
U = U.drop(columns=[c for c in ("rnaseq_cobertura", "rnaseq_lecturas", "rnaseq_soporte", "paralogo_resuelto",
                                "estabilidad_expresion", "estabilidad_estadios", "estabilidad_condiciones",
                                "log2FC_MRFV_4h", "log2FC_MRFV_7d", "id_GITV01", "soporte_GITV01", "confianza_gen",
                                "respuesta_MRFV", "loci_genoma", "paralogos_prot", "max_ident_paralogo") if c in U.columns])
gid = lambda s: s.replace("Dmai|", "")

# ---------- cobertura (Tarea B) ----------
ex = pd.read_csv(RS / "cobertura_exones_split.txt", sep="\t", header=None,
                 names=["chr", "ini", "fin", "gen", "n", "cubiertas", "largo", "frac"])
g = ex.groupby("gen")[["cubiertas", "largo"]].sum()
cob = (g.cubiertas / g.largo).rename("cob")
lect = sum(pd.read_csv(f, sep="\t", header=None, index_col=0)[1] for f in sorted((RS / "conteos").glob("SRR*.tsv")))
soporte = (cob.reindex(lect.index).fillna(0) >= MIN_COB) & (lect >= MIN_LECT)

# ---------- puente GITV01 (Tareas G, H) ----------
rbh = pd.read_csv("mapeo_helixer_GITV01.tsv", sep="\t")
rbh["gen"] = rbh.helixer.map(gid)
id_gitv = rbh.set_index("gen").id_GITV01
# locus: transcripto GITV01 alineado al genoma que cae sobre exones del gen
exb = defaultdict(list)
for r in ex.itertuples(index=False): exb[r.chr].append((r.ini, r.fin, r.gen))
for k in exb: exb[k].sort()
tx_gen, gen_tx = {}, defaultdict(set)
for row in csv.reader(open(RS / "gitv01_genoma_c.paf"), delimiter="\t"):
    q, ql, qs, qe, ch, ts, te, m, al = row[0], int(row[1]), int(row[2]), int(row[3]), row[5], int(row[7]), int(row[8]), int(row[9]), int(row[10])
    if qe <= qs or m / (qe - qs) < MIN_ID_LOCUS or (qe - qs) / ql < MIN_QCOV_LOCUS: continue  # identidad sobre la query (el bloque del PAF incluye intrones)
    solape = defaultdict(int)
    for a, b, gn in exb.get(ch, ()):
        if a > te: break
        if b > ts: solape[gn] += min(b, te) - max(a, ts)
    if solape:
        best = max(solape, key=solape.get)
        if solape[best] >= 100:
            tx_gen[q] = best; gen_tx[best].add(q)
soporte_gitv = pd.Series({gn: True for gn in set(id_gitv.index) | set(gen_tx)})

# ---------- MRFV (Tarea I) ----------
mrfv = pd.read_csv("mrfv/mrfv_DE.tsv", sep="\t")
gitv_rbh_inv = {v: k for k, v in id_gitv.items()}
# 1) locus genomico, 2) RBH de proteina, 3) mejor hit directo de proteina GITV01 -> Helixer (>= 80% id, cob 60/60)
g2h_best = {}
for r in csv.reader(open("score_insumos/g2h.tsv"), delimiter="\t"):
    t, b = r[0].rsplit(".p", 1)[0], float(r[7])
    if b > g2h_best.get(t, ("", 0))[1]: g2h_best[t] = (gid(r[1]), b)
mrfv["gen"] = [tx_gen.get(t) or gitv_rbh_inv.get(t) or g2h_best.get(t, (None,))[0] for t in mrfv.id_GITV01]
mrfv["via"] = ["locus" if t in tx_gen else "rbh" if t in gitv_rbh_inv else "hit_directo" if t in g2h_best else ""
               for t in mrfv.id_GITV01]
resp = mrfv.dropna(subset=["gen"]).assign(et=lambda d: d.sentido + "_" + d.tiempo) \
           .groupby("gen").et.apply(lambda s: ",".join(sorted(set(s))))

# ---------- estabilidad (Tarea E) ----------
EST = pd.read_csv("expresion/estadios_318849.tsv", sep="\t", index_col=0)
CON = pd.read_csv("expresion/condiciones_631706.tsv", sep="\t", index_col=0)
EST.index, CON.index = EST.index.map(gid), CON.index.map(gid)

# ---------- loci (Tarea F) ----------
loci = defaultdict(list)
for row in csv.reader(open(RS / "todos_loci_c.paf"), delimiter="\t"):
    q, ql, qs, qe, ch, ts, te, m, al = row[0], int(row[1]), int(row[2]), int(row[3]), row[5], int(row[7]), int(row[8]), int(row[9]), int(row[10])
    if qe > qs and m / (qe - qs) >= MIN_ID_LOCUS and (qe - qs) / ql >= MIN_QCOV_LOCUS:
        g0 = re.sub(r"\.\d+$", "", q)
        if not any(c == ch and ts < b and te > a for c, a, b in loci[g0]):
            loci[g0].append((ch, ts, te))

# ---------- paralogos a nivel proteina (complementa Tarea F) ----------
# otro gen Dmai con identidad >= 40% y cobertura >= 50% de ambas proteinas (DIAMOND --more-sensitive, e < 1e-10)
MIN_ID_PAR, MIN_COV_PAR = 40, 0.5
par = defaultdict(dict)
for r in csv.reader(open(RS / "dmai_vs_dmai.tsv"), delimiter="\t"):
    q, t, pid, ln, ql, sl = gid(r[0]), gid(r[1]), float(r[2]), int(r[3]), int(r[4]), int(r[5])
    if q != t and pid >= MIN_ID_PAR and ln / ql >= MIN_COV_PAR and ln / sl >= MIN_COV_PAR:
        par[q][t] = max(par[q].get(t, 0), pid)

# ---------- armar columnas ----------
rep = U.dmai_representativo.map(gid)
U["rnaseq_cobertura"] = rep.map(cob).round(3)
U["rnaseq_lecturas"] = rep.map(lect)
U["rnaseq_soporte"] = rep.map(soporte).fillna(False)
def resolver(row):
    if row.copias_dmai != 2: return ""
    s = [bool(soporte.get(gid(x), False)) for x in row.dmai_todos.split(",")]
    return {2: "ambas_con_lecturas", 1: "solo_una", 0: "ninguna"}[sum(s)]
U["paralogo_resuelto"] = U.apply(resolver, axis=1)
U["estabilidad_estadios"] = rep.map(EST.CV_estadios).round(3)
U["estabilidad_condiciones"] = rep.map(CON.CV_condiciones).round(3)
for c in ("log2FC_MRFV_4h", "log2FC_MRFV_7d"):
    U[c] = rep.map(CON[c]).round(3)
U["id_GITV01"] = rep.map(id_gitv).fillna("")
U["soporte_GITV01"] = rep.isin(soporte_gitv.index)
U["confianza_gen"] = [{(True, True): "alta", (True, False): "media", (False, True): "media_revisar",
                       (False, False): "baja"}[(bool(a), bool(b))] for a, b in zip(U.rnaseq_soporte, U.soporte_GITV01)]
U["respuesta_MRFV"] = rep.map(resp).fillna("")
U["loci_genoma"] = [len(loci[r]) if r in loci else np.nan for r in rep]
U["paralogos_prot"] = [len(par[r]) for r in rep]
U["max_ident_paralogo"] = [max(par[r].values()) if par[r] else np.nan for r in rep]
U.to_csv("ranking/ranking_completo.tsv", sep="\t", float_format="%.4g")

# ---------- reporte ----------
t100 = U.head(100)
allg = pd.Index(sorted(set(ex.gen)))
print(f"genes Helixer con soporte RNA-seq (cob>={MIN_COB}, lect>={MIN_LECT}): {int(soporte.reindex(allg).fillna(False).sum())}/{len(allg)}")
print(f"  top100 con soporte: {int(t100.rnaseq_soporte.sum())}/100")
print(f"puente GITV01: RBH {id_gitv.size}/{len(allg)}; por locus {len(gen_tx)}; cualquiera {soporte_gitv.size}")
print(f"  top100 con id_GITV01 (RBH): {(t100.id_GITV01 != '').sum()}; con soporte_GITV01: {int(t100.soporte_GITV01.sum())}")
print("confianza top100:", t100.confianza_gen.value_counts().to_dict())
print("confianza todos:", U.confianza_gen.value_counts().to_dict())
print("paralogo_resuelto (todos con 2 copias):", U[U.copias_dmai == 2].paralogo_resuelto.value_counts().to_dict())
print("paralogo_resuelto top100:", t100[t100.copias_dmai == 2].paralogo_resuelto.value_counts().to_dict())
print("paralogos_prot top100 (genes Dmai >=40% id, cob>=50%):", t100.paralogos_prot.clip(upper=3).value_counts().sort_index().to_dict(), "(3 = >=3)")
print("MRFV via:", mrfv.via.value_counts().to_dict())
print("loci_genoma top100:", t100.loci_genoma.value_counts(dropna=False).sort_index().to_dict())
print(f"MRFV: {len(mrfv)} transcriptos DE, mapeados a gen Helixer: {mrfv.gen.notna().sum()}; en el ranking: {(U.respuesta_MRFV != '').sum()}")
pc = lambda col, v: (U[col] <= v).mean()   # fraccion del ranking con CV <= (percentil; bajo = mas estable)
ctl = U[U.entrenamiento.fillna("").str.contains("PSMB5|Vha")]
for g, r in ctl.iterrows():
    print(f"estabilidad #{int(r['rank'])} {r.entrenamiento}: CV_estadios {r.estabilidad_estadios} (percentil {pc('estabilidad_estadios', r.estabilidad_estadios):.2f})"
          f" | CV_condiciones {r.estabilidad_condiciones} (percentil {pc('estabilidad_condiciones', r.estabilidad_condiciones):.2f})")
for c in ("estabilidad_estadios", "estabilidad_condiciones"):
    print(f"{c} mediana todos {U[c].median():.3f} / top100 {t100[c].median():.3f}")
print("Spearman CV_condiciones vs score:", round(U.estabilidad_condiciones.corr(U.score, method="spearman"), 3),
      "| CV_estadios vs CV_condiciones:", round(U.estabilidad_estadios.corr(U.estabilidad_condiciones, method="spearman"), 3))
print("\ntop100 con confianza baja:")
print(t100[t100.confianza_gen == "baja"][["rank", "dmel", "dmai_representativo", "rnaseq_cobertura", "rnaseq_lecturas", "TPM_minimo"]].to_string())
mrfv.to_csv("mrfv/mrfv_DE_mapeado.tsv", sep="\t", index=False)

# ---------- lista aparte de competencia vectorial (Tarea I) ----------
# genes del ranking con respuesta a MRFV (Tablas 3/4 de Xu et al. 2021: solo transcriptos ANOTADOS -> sesgo a genes
# conservados). Marcas: letalidad_norm >= 0.75; expresion estable = CV_estadios en el 25% mas bajo del ranking.
# log2FC de PRJNA631706 colapsado: descriptivos, sin umbral.
cv_p25 = U.estabilidad_estadios.quantile(0.25)
CVV = U[U.respuesta_MRFV != ""].assign(buena_letalidad=lambda d: d.letalidad_norm >= 0.75,
                                         expresion_estable=lambda d: d.estabilidad_estadios <= cv_p25)
CVV[["rank", "score", "respuesta_MRFV", "dmel", "letalidad_raw", "letalidad_norm", "estabilidad_estadios",
     "estabilidad_condiciones", "log2FC_MRFV_4h", "log2FC_MRFV_7d", "TPM_minimo", "confianza_gen", "copias_dmai",
     "paralogos_prot", "dmai_representativo", "id_GITV01", "buena_letalidad", "expresion_estable"]] \
    .to_csv("ranking/candidatos_competencia_vectorial.tsv", sep="\t", float_format="%.4g")
print(f"\ncompetencia vectorial: {len(CVV)} genes; con letalidad >= p75 o CV_estadios <= p25 ({cv_p25:.3f}): "
      f"{int((CVV.buena_letalidad | CVV.expresion_estable).sum())}")
print(CVV[CVV.buena_letalidad | CVV.expresion_estable][["rank", "dmel", "respuesta_MRFV", "letalidad_norm",
      "estabilidad_estadios", "log2FC_MRFV_4h", "log2FC_MRFV_7d"]].to_string())
