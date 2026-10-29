#!/usr/bin/env python3
"""
metamorfosis.py — brief_metamorfosis_y_cierre.md, tareas 0, 1, 2, 3 y 6 (la 4 es la regeneracion y la 5 sale de aca).

Correr DESPUES de la regeneracion completa (expresion_matrices -> score_candidatos [ogs3] -> comparar_ogs3 ->
rnaseq_integracion -> verificar_candidatos). Usa:
  ranking/ranking_completo.tsv, ranking_ogs3/ranking_completo.tsv  (tabla maestra; aca se le agregan columnas)
  ranking/go_por_gen.tsv, ranking_ogs3/go_por_gen.tsv               (GO de Drosophila propagados, por gen Tcas)
  ranking/entrenamiento.tsv, ranking/firma_go.tsv, ranking/candidatos_20_verificacion.tsv
  tabla_maestra/go-basic.obo (release fija 2026-07-26, la misma con la que se construyo la firma GO)
Salidas: resultados/metamorfosis/ y set_entrenamiento_GO.txt (junto al brief). El reporte va a stdout
(redirigir a resultados/metamorfosis/log_metamorfosis.txt).

Criterios (todos explicitos en el log):
  - Termino metamorfico ESTRICTO = descendiente (is_a, part_of) de una semilla o nombre con palabra clave, que no cae en
    la lista de exclusion. Exclusion: (a) muda/ecdisis/quitina/cuticula/embrion (transfieren a hemimetabolos);
    (b) nombre que aplica igual a un hemimetabolo: ecdysone, ecdysteroid, juvenile hormone, instar, larva(l), salvo que
    el nombre ademas diga metamorph/pupa/puparium/imaginal disc/wing disc/eclosion.
    Conjunto AMPLIO = solo con la exclusion (a); se reporta como sensibilidad.
  - Proteasomal (tarea 0) = el ortologo de Drosophila esta anotado (con propagacion) a GO:0000502 proteasome complex: si;
    anotado solo a catabolismo proteasomal (GO:0010498) sin ser parte del complejo: parcial; si no: no.
  - Generico (tarea 2) = profundidad <= 3 (camino minimo a la raiz por is_a/part_of) o anotado en > 25% del universo.
  - Alto (tarea 3) = >= percentil q de los valores OBSERVADOS (no imputados); q principal 0.75, alternativos 0.50 y 0.90.
    dali = dali22 observado; si dali22 es imputado, dali11 observado; si ambos imputados -> sin dato (no alto).
  - Ranking B (tarea 6) = sin GO; pesos redistribuidos proporcionalmente: 0.35/0.70, 0.20/0.70, 0.15/0.70 (= 50 / 28.6 / 21.4%).
"""
import re
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import sys
sys.setrecursionlimit(10000)

OUT = Path("resultados/metamorfosis"); OUT.mkdir(parents=True, exist_ok=True)
OBO = Path("tabla_maestra/go-basic.obo")

# ---------------- ontologia ----------------
nombre, ns, padres, obsoleto = {}, {}, defaultdict(set), set()
version = ""
cur = None
for line in open(OBO):
    line = line.rstrip("\n")
    if line.startswith("data-version:"): version = line.split(": ", 1)[1]
    if line.startswith("["): cur = None
    elif line.startswith("id: GO:"): cur = line[4:]
    elif cur is None: continue
    elif line.startswith("name: "): nombre[cur] = line[6:]
    elif line.startswith("namespace: "): ns[cur] = line[11:]
    elif line.startswith("is_a: "): padres[cur].add(line[6:16])
    elif line.startswith("relationship: part_of "): padres[cur].add(line[22:32])
    elif line == "is_obsolete: true": obsoleto.add(cur)
hijos = defaultdict(set)
for t, ps in padres.items():
    for p in ps: hijos[p].add(t)
RAICES = {"GO:0008150", "GO:0003674", "GO:0005575"}
prof = {r: 0 for r in RAICES}
def profundidad(t, pila=()):
    if t not in prof:
        ps = [p for p in padres.get(t, ()) if p not in pila]
        prof[t] = min(profundidad(p, pila + (t,)) for p in ps) + 1 if ps else 0
    return prof[t]
def descendientes(t):
    s, pend = set(), [t]
    while pend:
        x = pend.pop()
        for h in hijos.get(x, ()):
            if h not in s: s.add(h); pend.append(h)
    return s
print(f"go-basic.obo: {version} (archivo local, el mismo usado para la firma GO; no se volvio a descargar para no "
      f"desacoplar la firma del flag)")

# ---------------- TAREA 0: set de entrenamiento ----------------
REF = {
    "ulrich": "Ulrich J, Dao VA, Majumdar U, et al. (2015) Large scale RNAi screen in Tribolium reveals novel target genes "
              "for pest control and the proteasome as prime target. BMC Genomics 16:674 (Fig. 1 y Tabla S1; PMC4559001)",
    "knorr": "Knorr E, Fishilevich E, Tenbusch L, et al. (2018) Gene silencing in Tribolium castaneum as a tool for the "
             "targeted identification of candidate RNAi targets in crop pests. Sci Rep 8:2061 (PMC5794766); probado en "
             "T. castaneum, Diabrotica virgifera virgifera y Meligethes aeneus",
    "ledprona": "Rodrigues TB, Mishra SK, Sridharan K, et al. (2021) First sprayable double-stranded RNA-based biopesticide "
                "product targets proteasome subunit beta type-5 in Colorado potato beetle (Leptinotarsa decemlineata). "
                "Front Plant Sci 12:728652 (ledprona, blanco PSMB5)",
    "yao": "Yao J, Rotenberg D, Afsharifar A, Barandoc-Alviar K, Whitfield AE (2013) Development of RNAi methods for "
           "Peregrinus maidis, the corn planthopper. PLoS ONE 8(8) (PMC3737204); V-ATPasa B y D",
}
ENTR_INFO = {  # nombre en el codigo -> (simbolo, especie del experimento, id estable, ref)
    "cact": ("cact", "Tribolium castaneum", "TC002003", "ulrich"),
    "Srp54k": ("Srp54k", "Tribolium castaneum", "TC002574", "ulrich"),
    "Rop": ("Rop", "Tribolium castaneum", "TC011120", "ulrich"),
    "alphaSnap": ("alphaSnap", "Tribolium castaneum", "TC013571", "ulrich"),
    "shi": ("shi", "Tribolium castaneum", "TC011058", "ulrich"),
    "Pp1alpha-96A": ("Pp1alpha-96A", "Tribolium castaneum", "TC015321", "ulrich"),
    "Inr-a": ("Inr-a", "Tribolium castaneum", "TC008263", "ulrich"),
    "Hsc70-3": ("Hsc70-3", "Tribolium castaneum", "TC004425", "ulrich"),
    "Rpn7": ("Rpn7", "Tribolium castaneum", "TC006375", "ulrich"),
    "gw": ("gw", "Tribolium castaneum", "TC006679", "ulrich"),
    "Rpt3": ("Rpt3", "Tribolium castaneum", "TC007999", "ulrich"),
    "dre4": ("dre4", "T. castaneum / D. v. virgifera / M. aeneus", "FBgn0002183", "knorr"),
    "ncm": ("ncm", "T. castaneum / D. v. virgifera / M. aeneus", "FBgn0086707", "knorr"),
    "Polr2B/RpII140": ("RpII140 (Polr2B)", "T. castaneum / D. v. virgifera / M. aeneus", "FBgn0262955", "knorr"),
    "Prosbeta5 (PSMB5)": ("Prosbeta5 (PSMB5)", "Leptinotarsa decemlineata", "FBgn0029134", "ledprona"),
    "Vha55 (V-ATPasa B)": ("Vha55 (V-ATPasa B)", "Peregrinus maidis", "FBgn0005671", "yao"),
    "Vha36-1 (V-ATPasa D)": ("Vha36-1 (V-ATPasa D)", "Peregrinus maidis", "FBgn0022097", "yao"),
}
PROT_CX, PROT_CAT = "GO:0000502", "GO:0010498"
print(f"criterio proteasomal: {PROT_CX} = {nombre[PROT_CX]} -> si; solo {PROT_CAT} = {nombre[PROT_CAT]} -> parcial")

def leer_go(path):
    d = {}
    for line in open(path):
        if line.startswith("gen\t"): continue
        g, t = line.rstrip("\n").split("\t")
        d[g] = set(t.split(",")) if t else set()
    return d
GO = leer_go("ranking/go_por_gen.tsv")
def proteasomal(g):
    s = GO.get(str(g), set())
    return "si" if PROT_CX in s else "parcial" if PROT_CAT in s else "no"

E = pd.read_csv("ranking/entrenamiento.tsv", sep="\t", dtype={"gen": str})
E["proteasomal"] = E.gen.map(proteasomal)
firma = pd.read_csv("ranking/firma_go.tsv", sep="\t")
fondo_entr = E[E.en_fondo_1a1]
aporta = {t: [r.nombre for r in fondo_entr.itertuples() if t in GO[r.gen]] for t in firma.go}
prot_si = set(E.nombre[E.proteasomal == "si"])
solo_prot = {t for t, gs in aporta.items() if gs and set(gs) <= prot_si}

with open("set_entrenamiento_GO.txt", "w") as f:
    f.write("# Set de entrenamiento de la firma GO (17 genes) — generado por scripts/metamorfosis.py desde el codigo de\n"
            "# scripts/score_candidatos.py (ENTR_TC / ENTR_FB). Columnas separadas por tab.\n"
            f"# proteasomal: si = ortologo Dmel anotado a {PROT_CX} ({nombre[PROT_CX]}); parcial = solo {PROT_CAT} "
            f"({nombre[PROT_CAT]}); no = ninguno.\n"
            "# en_fondo_1a1: el gen entra al calculo de la firma solo si es 1:1 Tcas-Dmai y tiene GO (Rop no lo es).\n")
    f.write("simbolo\tespecie_experimento_RNAi\tid_estable\tgen_Tcas_RefSeq(GeneID)\tFBgn_ortologo_Dmel\t"
            "proteasomal\ten_fondo_1a1\treferencia\n")
    for r in E.itertuples():
        sim, esp, ide, ref = ENTR_INFO[r.nombre]
        f.write(f"{sim}\t{esp}\t{ide}\t{r.gen}\t{r.fbgn_dmel}\t{r.proteasomal}\t{r.en_fondo_1a1}\t{REF[ref]}\n")
    f.write(f"\n# proteasomales: {sorted(prot_si)} -> {len(prot_si)} de {len(E)}\n"
            f"# terminos de la firma (de {len(firma)}) cuyos genes de entrenamiento aportantes son TODOS proteasomales: "
            f"{len(solo_prot)} ({len(solo_prot) / len(firma):.1%})\n")
print("\n=== TAREA 0")
print(E[["nombre", "gen", "fbgn_dmel", "proteasomal", "en_fondo_1a1"]].to_string(index=False))
print(f"proteasomales: {len(prot_si)}/17 {sorted(prot_si)}; parcial: {sorted(E.nombre[E.proteasomal == 'parcial'])}")
print(f"terminos de la firma provenientes SOLO de genes proteasomales: {len(solo_prot)}/{len(firma)} "
      f"({len(solo_prot) / len(firma):.1%})")

# ---------------- TAREA 1: terminos metamorficos ----------------
SEMILLAS = {"GO:0007552": "metamorphosis", "GO:0002165": "instar larval or pupal development",
            "GO:0007444": "imaginal disc development", "GO:0007476": "imaginal disc-derived wing morphogenesis",
            "GO:0035220": "wing disc development", "GO:0007562": "eclosion", "GO:0035075": "response to ecdysone"}
CLAVES = ["metamorphosis", "pupa", "pupal", "puparium", "imaginal disc", "wing disc", "eclosion", "ecdysone",
          "ecdysteroid", "juvenile hormone", "instar"]
print("\n=== TAREA 1")
por_nombre = {v: k for k, v in nombre.items() if k not in obsoleto}
semillas = {}
for gid, n in SEMILLAS.items():
    if nombre.get(gid) == n and gid not in obsoleto:
        semillas[gid] = n; print(f"semilla {gid} resuelve a '{n}': OK")
    else:
        alt = por_nombre.get(n)
        print(f"semilla {gid}: en el .obo es '{nombre.get(gid)}'{' (obsoleto)' if gid in obsoleto else ''}; "
              f"resuelta por nombre -> {alt}")
        if alt: semillas[alt] = n
origen = defaultdict(set)
for s in semillas:
    for t in {s} | descendientes(s): origen[t].add(f"semilla:{s}")
for t, n in nombre.items():
    if t in obsoleto: continue
    for k in CLAVES:
        if re.search(r"\b" + re.escape(k), n, re.I): origen[t].add(f"nombre:{k}")
EXC_MUDA = re.compile(r"molt|moult|ecdysis|chitin|cuticl|embryo", re.I)
EXC_HEMI = re.compile(r"ecdysone|ecdysteroid|juvenile hormone|instar|larva", re.I)
ESTRICTO_OK = re.compile(r"metamorph|pupa|puparium|imaginal disc|wing disc|eclosion", re.I)
filas_t, filas_x = [], []
for t, o in origen.items():
    n = nombre[t]
    if EXC_MUDA.search(n):
        filas_x.append((t, n, ns[t], "muda/ecdisis/quitina/cuticula/embrion (transfiere a hemimetabolos)", ";".join(sorted(o))))
        continue
    hemi = bool(EXC_HEMI.search(n)) and not ESTRICTO_OK.search(n)
    if hemi:
        filas_x.append((t, n, ns[t], "aplica igual a un hemimetabolo (ecdisona/HJ/instar/larva sin marcador "
                        "metamorfico); queda solo en el conjunto amplio", ";".join(sorted(o))))
    filas_t.append((t, n, ns[t], "amplio" if hemi else "estricto", ";".join(sorted(o))))
TM = pd.DataFrame(filas_t, columns=["go", "nombre", "ontologia", "conjunto", "origen"]).sort_values(["conjunto", "go"],
                                                                                                   ascending=[False, True])
TM.to_csv(OUT / "terminos_metamorficos.tsv", sep="\t", index=False)
pd.DataFrame(filas_x, columns=["go", "nombre", "ontologia", "motivo", "origen"]).sort_values("go") \
  .to_csv(OUT / "exclusiones_aplicadas.tsv", sep="\t", index=False)
META_E = set(TM.go[TM.conjunto == "estricto"]); META_A = set(TM.go)
print(f"terminos metamorficos: estricto {len(META_E)}, amplio {len(META_A)}; excluidos {len(filas_x)} "
      f"(muda/quitina/cuticula/embrion: {sum(x[3].startswith('muda') for x in filas_x)}, aplica a hemimetabolo: "
      f"{sum(x[3].startswith('aplica') for x in filas_x)})")

def flag_tabla(ranking, gofile):
    U = pd.read_csv(ranking, sep="\t", index_col=0, dtype={0: str}, low_memory=False)
    U.index = U.index.astype(str)
    G = leer_go(gofile)
    te = [sorted(G.get(g, set()) & META_E) for g in U.index]
    ta = [sorted(G.get(g, set()) & META_A) for g in U.index]
    U["flag_metamorfosis"] = [bool(x) for x in te]
    U["terminos_metamorficos"] = [";".join(f"{t} {nombre[t]}" for t in x) for x in te]
    U["flag_metamorfosis_amplio"] = [bool(x) for x in ta]
    U["proteasomal"] = [proteasomal_g(G, g) for g in U.index]
    return U, G
def proteasomal_g(G, g):
    s = G.get(g, set())
    return "si" if PROT_CX in s else "parcial" if PROT_CAT in s else "no"

U, GU = flag_tabla("ranking/ranking_completo.tsv", "ranking/go_por_gen.tsv")
O, GO3 = flag_tabla("ranking_ogs3/ranking_completo.tsv", "ranking_ogs3/go_por_gen.tsv")
V20 = pd.read_csv("ranking/candidatos_20_verificacion.tsv", sep="\t", dtype={"tcas_geneid": str})
t100 = U[U["rank"] <= 100]
for n, X in (("RefSeq", U), ("OGS3", O)):
    print(f"{n}: marcados {int(X.flag_metamorfosis.sum())}/{len(X)} ({X.flag_metamorfosis.mean():.1%}); "
          f"amplio {int(X.flag_metamorfosis_amplio.sum())} ({X.flag_metamorfosis_amplio.mean():.1%}); "
          f"top100 {int(X[X['rank'] <= 100].flag_metamorfosis.sum())} (amplio {int(X[X['rank'] <= 100].flag_metamorfosis_amplio.sum())})")
print("top100 RefSeq marcados (estricto):")
print(t100[t100.flag_metamorfosis][["rank", "dmel", "terminos_metamorficos"]].to_string(max_colwidth=110))
f20 = U.loc[V20.tcas_geneid]
print(f"20 candidatos finales marcados: {int(f20.flag_metamorfosis.sum())}/20 "
      f"{list(f20[f20.flag_metamorfosis].dmel)} (amplio {int(f20.flag_metamorfosis_amplio.sum())})")
n_t = int(t100.flag_metamorfosis.sum())
print("CRITERIO DE LECTURA (brief): top100 <= 5 marcados -> objecion cuantificada, una linea en el poster; >= 20 -> "
      "discutir reponderar dali22; entre 5 y 20 -> inspeccion caso por caso. Resultado: "
      f"{n_t} -> " + ("<= 5" if n_t <= 5 else ">= 20" if n_t >= 20 else "entre 5 y 20"))

# ---------------- TAREA 2: auditoria de la firma ----------------
print("\n=== TAREA 2")
inter = sorted(set(firma.go) & META_A)
print(f"INTERSECCION firma GO x conjunto metamorfico (amplio): {len(inter)} {[(t, nombre[t]) for t in inter]}")
n_univ = Counter(t for g in U.index for t in GU.get(g, ()))
A = firma.copy()
A["nombre"] = A.go.map(nombre); A["ontologia"] = A.go.map(ns)
A["profundidad"] = A.go.map(profundidad)
A["genes_universo"] = A.go.map(n_univ).fillna(0).astype(int)
A["frac_universo"] = (A.genes_universo / len(U)).round(4)
A["genes_entrenamiento"] = A.go.map(lambda t: ",".join(aporta[t]))
A["solo_proteasomales"] = A.go.isin(solo_prot)
A["metamorfico"] = A.go.isin(META_A)
A["generico"] = (A.profundidad <= 3) | (A.frac_universo > 0.25)
A = A[["go", "nombre", "ontologia", "profundidad", "genes_universo", "frac_universo", "p", "k_entrenamiento",
       "k_fondo", "genes_entrenamiento", "solo_proteasomales", "metamorfico", "generico"]]
A.to_csv(OUT / "auditoria_firma_GO.tsv", sep="\t", index=False)
print(f"metamorficos: {int(A.metamorfico.sum())}/77; genericos: {int(A.generico.sum())}/77 "
      f"(prof <= 3: {int((A.profundidad <= 3).sum())}, > 25% universo: {int((A.frac_universo > 0.25).sum())})")
print(f"profundidad min/mediana/max: {A.profundidad.min()}/{A.profundidad.median()}/{A.profundidad.max()}")
print(f"solo de genes proteasomales: {int(A.solo_proteasomales.sum())}/77 ({A.solo_proteasomales.mean():.1%})")

# ---------------- TAREA 3: letal en ambos regimenes ----------------
print("\n=== TAREA 3")
def obs(X):
    d22 = X.dali22.where(~X.dali22_imputado.astype(bool))
    d11 = X.dali11.where(~X.dali11_imputado.astype(bool))
    return d22.fillna(d11), X.dapi11.where(~X.dapi11_imputado.astype(bool))
dali_o, dapi_o = obs(U)
def ambos(q, X=U):
    a, b = obs(X)
    ua, ub = a.quantile(q), b.quantile(q)
    return (a >= ua) & (b >= ub), ua, ub
lines = []
def rep(s): print(s); lines.append(s)
rep(f"observados: dali {dali_o.notna().sum()} genes (dali22 {int((~U.dali22_imputado.astype(bool)).sum())}, resto dali11); "
    f"dapi11 {dapi_o.notna().sum()}")
res = {}
for q in (0.75, 0.50, 0.90):
    m, ua, ub = ambos(q)
    res[q] = m
    top = m[U["rank"] <= 100]
    rho = spearmanr(U.score, m.astype(int)).correlation
    F = U[m].sort_values("score", ascending=False)
    sobrev = int(top.sum())
    ct = pd.crosstab(~m, U.flag_metamorfosis)
    no_y_meta = int(((~m) & U.flag_metamorfosis).sum()); meta = int(U.flag_metamorfosis.sum())
    rep(f"\nq = {q}{' (PRINCIPAL)' if q == 0.75 else ''}: umbral dali >= {ua:.1f}, dapi11 >= {ub:.1f}; "
        f"letal_ambos = {int(m.sum())}/{len(U)} ({m.mean():.1%})")
    rep(f"  top100 con letal_ambos: {sobrev}; Spearman(score, letal_ambos) = {rho:.3f}")
    rep(f"  letal_ambos=False y flag_metamorfosis=True: {no_y_meta} de {meta} marcados "
        f"({no_y_meta / max(meta, 1):.1%}); de los {int((~m).sum())} False, {no_y_meta / max((~m).sum(), 1):.1%} marcados "
        f"(vs {U[m].flag_metamorfosis.mean():.1%} entre los True)")
    rep(f"  ranking filtrado: {len(F)} genes; del top100 original sobreviven {sobrev}; el top100 filtrado llega hasta "
        f"el rank original {int(F['rank'].iloc[min(99, len(F) - 1)])}; Spearman rank filtrado vs original en los "
        f"sobrevivientes = 1 por construccion (el filtro no reordena)")
    if q == 0.75:
        F = F.assign(rank_filtro=np.arange(1, len(F) + 1))
        F[["rank_filtro", "rank", "score", "dmel", "dali11", "dali22", "dapi11", "flag_metamorfosis", "proteasomal"]] \
            .to_csv(OUT / "ranking_filtro_ambos_regimenes.tsv", sep="\t", float_format="%.4g")
U["letal_ambos_regimenes"] = res[0.75]
O["letal_ambos_regimenes"] = ambos(0.75, O)[0]
(OUT / "comparacion_filtro.txt").write_text(
    "Tarea 3 — letal_ambos_regimenes (diagnostico; NO reemplaza al ranking oficial)\n"
    "alto = >= percentil q de los valores observados (no imputados). dali = dali22 observado o, si es imputado, dali11 "
    "observado.\n" + "\n".join(lines) + "\n")

# ---------------- TAREA 6: ranking B sin GO ----------------
print("\n=== TAREA 6")
WB = {"letalidad": 0.35 / 0.70, "expresion": 0.20 / 0.70, "conservacion": 0.15 / 0.70}
print("pesos B (redistribucion proporcional):", {k: round(v, 4) for k, v in WB.items()})
U["score_B"] = sum(WB[c] * U[c + "_norm"] for c in WB)
U["rank_B"] = U.score_B.rank(ascending=False, method="first").astype(int)
tA, tB = U[U["rank"] <= 100], U[U.rank_B <= 100]
rhoAB = spearmanr(U["rank"], U.rank_B).correlation
print(f"top100 A en top100 B: {len(set(tA.index) & set(tB.index))}/100; Spearman A vs B = {rhoAB:.3f}")
print(f"proteasomales (GO:0000502) top100 A: {int((tA.proteasomal == 'si').sum())} | top100 B: "
      f"{int((tB.proteasomal == 'si').sum())}  (parcial: A {int((tA.proteasomal == 'parcial').sum())}, "
      f"B {int((tB.proteasomal == 'parcial').sum())})")
CTL = {"PSMB5": "Prosbeta5", "Rpt3": "Rpt3", "V-ATPasa D": "Vha36-1", "V-ATPasa B": "Vha55"}
for k, v in CTL.items():
    r = U[U.entrenamiento.fillna("").str.startswith(v)].iloc[0]
    print(f"  control {k}: A #{int(r['rank'])} -> B #{int(r.rank_B)}")
solo_B = tB[~tB.index.isin(tA.index)].sort_values("rank_B")
print(f"solo en top100 B: {len(solo_B)}")
print(solo_B[["rank_B", "rank", "dmel", "go_norm", "letalidad_norm", "proteasomal", "flag_metamorfosis"]].to_string())
print(f"top100 B con flag_metamorfosis: {int(tB.flag_metamorfosis.sum())} (A: {int(tA.flag_metamorfosis.sum())})")
U.sort_values("rank_B")[["rank_B", "rank", "score_B", "score", "dmel", "letalidad_norm", "expresion_norm",
                         "conservacion_norm", "go_norm", "proteasomal", "flag_metamorfosis", "letal_ambos_regimenes",
                         "entrenamiento"]].to_csv(OUT / "ranking_B_sin_GO.tsv", sep="\t", float_format="%.4g")

# ---------------- tabla maestra con columnas nuevas ----------------
for X, p in ((U, "ranking/ranking_completo.tsv"), (O, "ranking_ogs3/ranking_completo.tsv")):
    X.to_csv(p, sep="\t", float_format="%.4g")
print("\ncolumnas agregadas a la tabla maestra: flag_metamorfosis, terminos_metamorficos, flag_metamorfosis_amplio, "
      "proteasomal, letal_ambos_regimenes (+ score_B, rank_B en RefSeq)")

# ---------------- txt de la firma GO con el origen de cada termino ----------------
with open("resultados/firma_GO_entrenamiento.txt", "w") as f:
    f.write("FIRMA GO DEL SCORE (componente GO, 30%) Y SU ORIGEN\n"
            "=====================================================\n"
            "Construccion (scripts/score_candidatos.py): GO de Drosophila (FlyBase, todas las evidencias) transferidos al gen\n"
            "de Tribolium via ortologo Tcas-Dmel de OrthoFinder y propagados a ancestros (is_a, part_of; go-basic.obo\n"
            f"{version}). Test hipergeometrico (cola superior) de los genes de entrenamiento que estan en el fondo contra el\n"
            "FONDO = genes Tcas con ortologo 1:1 en D. maidis y con GO (5,739). Entra a la firma un termino con p < 0.01 y\n"
            ">= 2 genes de entrenamiento. Puntaje GO de un gen = suma de -log10(p) de los terminos de la firma que porta; los\n"
            "genes de entrenamiento se puntuan con la firma recalculada sin ellos (leave-one-out).\n\n"
            "1) GENES DE ENTRENAMIENTO (17; 16 en el fondo 1:1 — Rop no, por eso no aporta a ningun termino)\n")
    for r in E.itertuples():
        sim, esp, ide, ref = ENTR_INFO[r.nombre]
        f.write(f"  - {sim:22s} {ide:12s} especie del RNAi: {esp}; proteasomal: {r.proteasomal}; "
                f"en fondo: {r.en_fondo_1a1}\n      fuente: {REF[ref]}\n")
    f.write(f"\n2) LOS {len(A)} TERMINOS (ordenados por p). Columnas: GO | nombre | ontologia | p | genes de entrenamiento "
            "que lo portan (k) | genes del fondo que lo portan | profundidad | generico | solo proteasomales\n")
    for r in A.sort_values("p").itertuples():
        f.write(f"  {r.go} | {r.nombre} | {r.ontologia} | p={r.p:.2e} | {r.genes_entrenamiento} ({r.k_entrenamiento}) | "
                f"{r.k_fondo} | prof {r.profundidad} | {'GENERICO' if r.generico else '-'} | "
                f"{'SOLO_PROTEASOMA' if r.solo_proteasomales else '-'}\n")
    cnt = Counter(g for s in A.genes_entrenamiento for g in s.split(",") if g)
    f.write("\n3) EN CUANTOS TERMINOS DE LA FIRMA APARECE CADA GEN DE ENTRENAMIENTO\n")
    for g, n in cnt.most_common():
        f.write(f"  {g:24s} {n}\n")
    f.write(f"\nResumen: {len(prot_si)}/17 genes proteasomales ({', '.join(sorted(prot_si))}); terminos que provienen SOLO de "
            f"genes proteasomales: {len(solo_prot)}/{len(A)} ({len(solo_prot) / len(A):.1%}); terminos genericos "
            f"(profundidad <= 3 o > 25% del universo): {int(A.generico.sum())}/{len(A)}; metamorficos: "
            f"{int(A.metamorfico.sum())}.\n")
