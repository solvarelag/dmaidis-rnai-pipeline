#!/usr/bin/env python3
"""
score_candidatos.py — ranking de blancos de RNAi en Dalbulus maidis a partir de iBeetle (Tribolium).

Universo: genes de Tribolium (GeneID RefSeq, GCF_000002335.3) con dato de letalidad en iBeetle (sin corte)
y con >= 1 ortologo en D. maidis (OrthoFinder pasada 1).

Copias en D. maidis (n de genes Dmai ortologos del gen Tcas):
    1 -> pasa;  2 -> pasa con flag riesgo_paralogo;  >= 3 -> se descarta.

Score = 0.35 letalidad + 0.30 firma GO + 0.20 expresion constitutiva + 0.15 conservacion

NORMALIZACION: cada columna se lleva a 0-1 por RANGO PERCENTIL dentro del universo filtrado
(rank promedio para empates; (rank-1)/(n-1)). Se eligio percentil y no min-max porque las cuatro variables
tienen colas largas (TPM, -log10 p del GO) y con min-max un solo outlier comprime al resto. Consecuencia:
el score refleja posiciones relativas, no magnitudes absolutas.

Letalidad = 0.5 dali22 + 0.3 dali11 + 0.2 dapi11 (valores por gen = promedio de sus iB).
Faltantes: dali11 y dapi11 se imputan con el promedio de la columna sobre los genes del universo que la tienen (no cero).
dali22 se imputa con max(dali11, promedio de dali22): la mortalidad es acumulada (dali22 >= dali11 en el 99.9% de las iB
con ambas medidas) y el promedio solo daria valores imposibles a los genes que matan rapido.
Sensibilidad: el top 100 se recalcula con letalidad sin dali22 (0.6 dali11 + 0.4 dapi11) -> ranking/sensibilidad_sin_dali22.tsv.

Firma GO: GO de Drosophila transferidos via ortologo Tcas-Dmel (OrthoFinder), propagados a ancestros
(is_a, part_of; go-basic.obo). Enriquecimiento: test hipergeometrico (cola superior) del conjunto de
entrenamiento contra el FONDO = genes Tcas con ortologo 1:1 en D. maidis. Terminos de la firma: p < 0.01 y
>= 2 genes de entrenamiento. Puntaje GO de un gen = suma de -log10(p) de los terminos de la firma que porta.
Los genes de entrenamiento se puntuan con leave-one-out (firma recalculada sin ellos) para no ser circulares.

Expresion constitutiva: SOLO PRJNA318849 (ninfa, hembra, macho; expresion/estadios_318849.tsv, de
scripts/expresion_matrices.py). Subscore = promedio de percentil(log10 TPM_minimo) y percentil inverso de CV_estadios;
despues se re-normaliza por percentil. sesgo_sexo = log2(hembra/macho) se informa y no entra al score.
PRJNA631706 (MRFV, replicas tecnicas) no se usa en el score.

Conservacion: % identidad (DIAMOND, mejor HSP) entre el gen Tcas y su ortologo Dmai representativo
(el de mayor identidad si hay 2 copias). Expresion se toma del mismo representativo.

Uso: python3 scripts/score_candidatos.py [ogs3]  (desde la raiz del proyecto). Salidas en ranking/.
Modo ogs3 (analisis separado): Tribolium = proteoma OGS3 (IDs TC, isoforma mas larga) en lugar de RefSeq, en una corrida
de OrthoFinder aparte con las mismas otras 9 especies (orthofinder_ogs3/). Los genes Tcas son TC: no hace falta el mapeo
TC -> GeneID. Salidas en ranking_ogs3/.
"""
import csv, glob, gzip, math, re, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata

T = Path("tabla_maestra")
OGS3 = "ogs3" in sys.argv[1:]
OF = Path("orthofinder_ogs3" if OGS3 else "orthofinder_pasada1")
SC = Path("score_insumos")          # salidas del job s4 copiadas del cluster
OUT = Path("ranking_ogs3" if OGS3 else "ranking"); OUT.mkdir(exist_ok=True)
PESOS = {"letalidad": 0.35, "go": 0.30, "expresion": 0.20, "conservacion": 0.15}
P_FIRMA, MIN_ENTR = 0.01, 2


def pct(x):
    """Rango percentil 0-1 (empates = rango promedio). NaN queda NaN."""
    x = pd.Series(x, dtype=float)
    ok = x.notna()
    r = pd.Series(np.nan, index=x.index)
    if ok.sum() > 1:
        r[ok] = (rankdata(x[ok]) - 1) / (ok.sum() - 1)
    return r


def leer_fasta_ids(path, patron):
    d = {}
    for line in open(path):
        if line.startswith(">"):
            m = re.match(patron, line)
            d[m.group(1)] = m.group(2)
    return d


# ---------------- mapeo iB -> TC -> GeneID ----------------
ib_tc = defaultdict(set)
for line in open(T / "iB_TC.csv"):
    if not line.startswith("#"):
        ib, tc = line.strip().split(",")[:2]
        ib_tc[ib].add(tc)

tc_gid, fuente_map = defaultdict(set), Counter()
for gid, tc, _ in csv.reader(open(T / "Tcas_GeneID_TC.tsv"), delimiter="\t"):
    if tc: tc_gid[tc].add(gid)
for f in ("Tcas_TC_GeneID_secuencia.tsv", "Tcas_TC_GeneID_rbh.tsv"):
    if (T / f).exists():
        for tc, gid in csv.reader(open(T / f), delimiter="\t"):
            tc_gid[tc].add(gid)
if OGS3:   # el gen Tcas ya es el TC
    tc_gid = {tc: {tc} for tcs in ib_tc.values() for tc in tcs}

# ---------------- letalidad por gen ----------------
let = pd.read_csv(T / "iBeetle_lethality.csv", sep=";", encoding="utf-8-sig")
# iBeetle codifica "no medido" como -1 (ademas de celdas vacias)
for c in ("dali11", "dali22", "dapi11"):
    let.loc[let[c] < 0, c] = np.nan
filas = []
for r in let.itertuples(index=False):
    for tc in ib_tc.get(r.iB, ()):
        for g in tc_gid.get(tc, ()):
            filas.append((g, r.iB, r.dali11, r.dali22, r.dapi11))
let_g = pd.DataFrame(filas, columns=["gid", "iB", "dali11", "dali22", "dapi11"])
ibs_por_gen = let_g.groupby("gid")["iB"].apply(lambda s: ",".join(sorted(set(s))))
let_g = let_g.groupby("gid")[["dali11", "dali22", "dapi11"]].mean()
let_g["iB"] = ibs_por_gen

# ---------------- ortologia Dmai <-> Tcas ----------------
tcas_prot_gid = leer_fasta_ids(OF / "Tribolium_castaneum.faa", r">((TC\d+))") if OGS3 else \
                leer_fasta_ids("proteomas/Tribolium_castaneum.faa", r">(\S+).*GeneID:(\d+)")
gid_tcasprot = {g: p for p, g in tcas_prot_gid.items()}
gid_dmai, uno_a_uno = defaultdict(set), set()
for og, dm, tc in csv.reader(open(OF / "Orthologues_Dalbulus_maidis/Dalbulus_maidis__v__Tribolium_castaneum.tsv"),
                             delimiter="\t"):
    if og == "Orthogroup": continue
    dms = [x.strip() for x in dm.split(",")]
    tcs = [x.strip() for x in tc.split(",")]
    for p in tcs:
        gid_dmai[tcas_prot_gid[p]] |= set(dms)
    if len(dms) == 1 and len(tcs) == 1:
        uno_a_uno.add(tcas_prot_gid[tcs[0]])

# ---------------- universo y filtro de copias ----------------
U = let_g[let_g.index.isin(gid_dmai.keys())].copy()
U["copias_dmai"] = [len(gid_dmai[g]) for g in U.index]
n_universo = len(U)
dist_copias = U["copias_dmai"].clip(upper=3).value_counts().sort_index()
U = U[U["copias_dmai"] <= 2].copy()
U["riesgo_paralogo"] = U["copias_dmai"] == 2

# imputacion por media (sobre el universo filtrado)
imput = {}
for c in ("dali11", "dapi11", "dali22"):
    falta = U[c].isna()
    imput[c] = int(falta.sum())
    U[c + "_imputado"] = falta
    if c == "dali22":   # monotonia: dali22 >= dali11
        U[c] = U[c].fillna(np.maximum(U.dali11, U[c].mean()))
    else:
        U[c] = U[c].fillna(U[c].mean())
U["n_imputadas"] = U[["dali11_imputado", "dali22_imputado", "dapi11_imputado"]].sum(axis=1)
U["letalidad_raw"] = 0.5 * U.dali22 + 0.3 * U.dali11 + 0.2 * U.dapi11
U["letalidad_sin_dali22"] = 0.6 * U.dali11 + 0.4 * U.dapi11

# ---------------- GO: Dmel -> Tcas, propagacion ----------------
padres = defaultdict(set)
cur, obsoleto = None, False
for line in open(T / "go-basic.obo"):
    line = line.strip()
    if line == "[Term]": cur = None
    elif line.startswith("id: GO:"): cur = line[4:]
    elif cur and line.startswith("is_a: "): padres[cur].add(line[6:16])
    elif cur and line.startswith("relationship: part_of "): padres[cur].add(line[22:32])
anc_cache = {}
def ancestros(t):
    if t not in anc_cache:
        s = {t}
        for p in padres.get(t, ()): s |= ancestros(p)
        anc_cache[t] = s
    return anc_cache[t]
sys.setrecursionlimit(10000)

fbgn_go = defaultdict(set)
for fb, go, asp in csv.reader(open(T / "Dmel_FBgn_GO.tsv"), delimiter="\t"):
    fbgn_go[fb] |= ancestros(go)
dmel_gid_fbgn = dict(csv.reader(open(T / "Dmel_GeneID_FBgn.tsv"), delimiter="\t"))
dmel_prot_gid = leer_fasta_ids("proteomas/Drosophila_melanogaster.faa", r">(\S+).*GeneID:(\d+)")
tcas_dmel = defaultdict(set)
for og, tc, dm in csv.reader(open(OF / "Tribolium_castaneum__v__Drosophila_melanogaster.tsv"), delimiter="\t"):
    if og == "Orthogroup": continue
    fbs = {dmel_gid_fbgn.get(dmel_prot_gid[p.strip()]) for p in dm.split(",")} - {None}
    for p in tc.split(","):
        tcas_dmel[tcas_prot_gid[p.strip()]] |= fbs
gid_go = {g: set().union(*[fbgn_go[f] for f in fbs]) for g, fbs in tcas_dmel.items() if fbs}
fb_sym = dict(csv.reader(open(T / "Dmel_FBgn_simbolo.tsv"), delimiter="\t"))

# ---------------- conjunto de entrenamiento ----------------
ENTR_TC = {  # Ulrich et al. 2015, Fig. 1 + Tabla S1
    "TC002003": "cact", "TC002574": "Srp54k", "TC011120": "Rop", "TC013571": "alphaSnap",
    "TC011058": "shi", "TC015321": "Pp1alpha-96A", "TC008263": "Inr-a", "TC004425": "Hsc70-3",
    "TC006375": "Rpn7", "TC006679": "gw", "TC007999": "Rpt3"}
ENTR_FB = {  # Knorr et al. 2018 (Diabrotica), PSMB5 (ledprona), V-ATPasa B y D (Yao et al. 2013, P. maidis)
    "FBgn0002183": "dre4", "FBgn0086707": "ncm", "FBgn0262955": "Polr2B/RpII140",
    "FBgn0029134": "Prosbeta5 (PSMB5)", "FBgn0005671": "Vha55 (V-ATPasa B)", "FBgn0022097": "Vha36-1 (V-ATPasa D)"}
entr = {}
for tc, n in ENTR_TC.items():
    for g in tc_gid.get(tc, ()): entr[g] = n
for g, fbs in tcas_dmel.items():
    for f in fbs & ENTR_FB.keys():
        entr.setdefault(g, ENTR_FB[f])
faltan_entr = [n for n in list(ENTR_TC.values()) + list(ENTR_FB.values()) if n not in entr.values()]

fondo = {g for g in uno_a_uno if g in gid_go}
def firma(excluir=frozenset()):
    fg = [g for g in entr if g in fondo and g not in excluir]
    N, n = len(fondo), len(fg)
    cnt_fondo = Counter(t for g in fondo for t in gid_go[g])
    cnt_fg = Counter(t for g in fg for t in gid_go[g])
    res = {}
    for t, k in cnt_fg.items():
        if k >= MIN_ENTR:
            p = hypergeom.sf(k - 1, N, cnt_fondo[t], n)
            if p < P_FIRMA: res[t] = (p, k, cnt_fondo[t])
    return res, n
FIRMA, n_entr_fondo = firma()
def go_score(g, f):
    return sum(-math.log10(f[t][0]) for t in gid_go.get(g, ()) if t in f)
U["go_raw"] = [go_score(g, firma({g})[0] if g in entr else FIRMA) for g in U.index]
U["entrenamiento"] = [entr.get(g, "") for g in U.index]

# ---------------- conservacion y expresion ----------------
ident = defaultdict(float)
for row in csv.reader(open(SC / ("dmai_vs_ogs3.tsv" if OGS3 else "dmai_vs_tcas.tsv")), delimiter="\t"):
    q, s, pid = row[0], row[1], float(row[2])
    ident[(q, s)] = max(ident[(q, s)], pid)
rep, cons = {}, {}
for g in U.index:
    tp = gid_tcasprot.get(g)
    mejor = max(gid_dmai[g], key=lambda d: ident.get((d, tp), 0))
    rep[g], cons[g] = mejor, ident.get((mejor, tp), np.nan)
U["dmai_representativo"] = pd.Series(rep)
U["dmai_todos"] = [",".join(sorted(gid_dmai[g])) for g in U.index]
U["conservacion_raw"] = pd.Series(cons)

EST = pd.read_csv("expresion/estadios_318849.tsv", sep="\t", index_col=0).reindex(U.dmai_representativo)
for c in ("tpm_ninfa", "tpm_hembra", "tpm_macho", "TPM_minimo", "CV_estadios", "sesgo_sexo"):
    U[c] = EST[c].values
U["expresion_raw"] = 0.5 * pct(np.log10(U.TPM_minimo + 0.01)).values + 0.5 * (1 - pct(U.CV_estadios)).values

# ---------------- normalizacion y score ----------------
for c in PESOS:
    U[c + "_norm"] = pct(U[c + "_raw"]).fillna(0).values   # sin dato -> 0 en esa columna (se reporta)
for c in PESOS:
    U[c + "_aporte"] = PESOS[c] * U[c + "_norm"]
U["score"] = sum(U[c + "_aporte"] for c in PESOS)
# sensibilidad: mismo score con letalidad sin dali22
U["score_sin_dali22"] = U.score - U.letalidad_aporte + PESOS["letalidad"] * pct(U.letalidad_sin_dali22).fillna(0).values
U["rank_sin_dali22"] = rankdata(-U.score_sin_dali22, method="ordinal").astype(int)
U = U.sort_values("score", ascending=False)
U["rank"] = np.arange(1, len(U) + 1)
U["tc"] = [",".join(sorted({tc for ib in U.at[g, "iB"].split(",") for tc in ib_tc[ib]})) for g in U.index]
U["dmel"] = [",".join(sorted(fb_sym.get(f, f) for f in tcas_dmel.get(g, ()))) for g in U.index]
U.index.name = "tcas_tc" if OGS3 else "tcas_geneid"

cols = (["rank", "score"] + [c + "_aporte" for c in PESOS] + [c + "_norm" for c in PESOS] +
        ["letalidad_raw", "dali11", "dali22", "dapi11", "dali11_imputado", "dali22_imputado", "dapi11_imputado",
         "n_imputadas", "go_raw", "conservacion_raw",
         "TPM_minimo", "CV_estadios", "sesgo_sexo", "tpm_ninfa", "tpm_hembra", "tpm_macho",
         "rank_sin_dali22", "copias_dmai", "riesgo_paralogo", "entrenamiento",
         "tc", "iB", "dmel", "dmai_representativo", "dmai_todos"])
U[cols].to_csv(OUT / "ranking_completo.tsv", sep="\t", float_format="%.4g")
# insumos para scripts/metamorfosis.py: GO propagados por gen (universo + entrenamiento) y set de entrenamiento
with open(OUT / "go_por_gen.tsv", "w") as f:
    f.write("gen\tgo_propagados\n")
    for g in sorted(set(U.index) | set(entr)):
        f.write(f"{g}\t{','.join(sorted(gid_go.get(g, ())))}\n")
pd.DataFrame([(g, n, ",".join(sorted(tcas_dmel.get(g, ()))), g in fondo, g in U.index) for g, n in entr.items()],
             columns=["gen", "nombre", "fbgn_dmel", "en_fondo_1a1", "en_ranking"]).to_csv(OUT / "entrenamiento.tsv",
                                                                                         sep="\t", index=False)
pd.DataFrame([(t, fb, p, k, K) for t, (p, k, K) in sorted(FIRMA.items(), key=lambda x: x[1][0])
              for fb in [""]], columns=["go", "_", "p", "k_entrenamiento", "k_fondo"]).drop(columns="_") \
  .to_csv(OUT / "firma_go.tsv", sep="\t", index=False)

# ---------------- reporte ----------------
print(f"universo (letalidad + ortologo Dmai): {n_universo}")
print("copias Dmai:", {int(k): int(v) for k, v in dist_copias.items()}, "(3 = >=3)")
print(f"tras filtro de copias: {len(U)} (riesgo_paralogo: {int(U.riesgo_paralogo.sum())})")
print("imputados por columna:", imput, "| genes con >=1 imputacion:", int((U.n_imputadas > 0).sum()))
print(f"entrenamiento mapeado: {len(entr)} genes Tcas; en fondo 1:1 con GO: {n_entr_fondo}; no mapeados: {faltan_entr}")
print(f"fondo GO (1:1 con GO): {len(fondo)}; terminos en la firma: {len(FIRMA)}")
print("sin dato -> 0:", {c: int(U[c + '_raw'].isna().sum()) for c in PESOS})
print("expresion: PRJNA318849 (3 corridas); sesgo_sexo |log2| > 2 en el universo:", int((U.sesgo_sexo.abs() > 2).sum()))
t100, s100 = set(U.index[:100]), set(U.sort_values("rank_sin_dali22").index[:100])
print(f"\nsensibilidad sin dali22: top100 compartido {len(t100 & s100)}/100; "
      f"Spearman rank (todos) {U['rank'].corr(U.rank_sin_dali22, method='spearman'):.3f}")
S = U[["rank", "rank_sin_dali22", "dmel", "dali11", "dali22", "dali22_imputado", "dapi11"]]
S = S[(S["rank"] <= 100) | (S.rank_sin_dali22 <= 100)].sort_values("rank")
S.to_csv(OUT / "sensibilidad_sin_dali22.tsv", sep="\t", float_format="%.4g")
print("salen del top100 sin dali22:", S[(S['rank'] <= 100) & (S.rank_sin_dali22 > 100)][['rank', 'rank_sin_dali22', 'dmel', 'dali11', 'dali22']].to_string())
print("entran:", S[(S['rank'] > 100) & (S.rank_sin_dali22 <= 100)][['rank', 'rank_sin_dali22', 'dmel', 'dali11', 'dali22']].to_string())
print("\nentrenamiento / controles en el ranking:")
for g in U.index[U.entrenamiento != ""]:
    r = U.loc[g]
    print(f"  #{int(r['rank']):5d}/{len(U)}  {r.entrenamiento:24s} score={r.score:.3f} "
          f"L={r.letalidad_norm:.2f} GO={r.go_norm:.2f} E={r.expresion_norm:.2f} C={r.conservacion_norm:.2f}"
          f"{'  riesgo_paralogo' if r.riesgo_paralogo else ''}")
for g, n in entr.items():
    if g in U.index: continue
    if g not in let_g.index: motivo = "sin dato de letalidad en iBeetle"
    elif g not in gid_dmai: motivo = "sin ortologo en D. maidis"
    else: motivo = f"descartado por copias ({len(gid_dmai[g])} en Dmai)"
    print(f"  fuera     {n:24s} GeneID {g}: {motivo}")
print("\nreparto del score (aporte medio en el top 50 vs. todos):")
for c in PESOS:
    print(f"  {c:13s} top50={U[c + '_aporte'].head(50).mean():.3f}  todos={U[c + '_aporte'].mean():.3f}")
