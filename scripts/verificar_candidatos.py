#!/usr/bin/env python3
"""
verificar_candidatos.py — brief_rnaseq Tarea D sobre los 20 candidatos finales.

Candidatos finales = primeros 20 del ranking que cumplen:
  confianza_gen != baja; paralogo_resuelto != ambas_con_lecturas; loci_genoma == 1;
  sin paralogo Dmai con identidad proteica >= 70% (max_ident_paralogo < 70): riesgo de compensacion y de co-silenciamiento.

Para cada candidato se compara el modelo de Helixer con los transcriptos de StringTie de novo (BAM de 19 corridas,
sin anotacion guia) del mismo locus y hebra:
  intrones_helixer = intrones dentro del CDS de Helixer (los de UTR no afectan al dsRNA y Helixer los predice peor);
  intrones_confirmados = de esos, los que aparecen exactos (mismo donor/acceptor) en ALGUN transcripto StringTie del locus
  (= hay lecturas que saltan exactamente esa union);
  cds_cubierto (fraccion del CDS de Helixer contenida en exones del transcripto StringTie).
Transcripto StringTie elegido para la secuencia = el de mayor (intrones del CDS que comparte, cds_cubierto).
Estado: confirmado (todos los intrones y CDS >= 0.95) / discrepancia / sin_stringtie.
FASTA: secuencia del transcripto StringTie (ARNm respaldado por lecturas) si cubre >= 90% del CDS de Helixer, aunque
discrepe en intrones (las lecturas mandan; queda marcado para revision); si no, el ARNm de Helixer marcado como no verificado.
"""
import csv, re
from collections import defaultdict
import pandas as pd

RS = "rnaseq_insumos/"
U = pd.read_csv("ranking/ranking_completo.tsv", sep="\t", index_col=0)
ok = (U.confianza_gen != "baja") & (U.paralogo_resuelto.fillna("") != "ambas_con_lecturas") & \
     (U.loci_genoma == 1) & ~(U.max_ident_paralogo >= 70)
F = U[ok].head(20).copy()
excluidos = U.head(int(F["rank"].max()))[~ok.reindex(U.head(int(F["rank"].max())).index)]

def exones_gff(path, tipo_padre_re, feat):
    d = defaultdict(list); info = {}
    for line in open(path):
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != feat: continue
        m = re.search(tipo_padre_re, c[8])
        if m:
            d[m.group(1)].append((int(c[3]), int(c[4])))
            info[m.group(1)] = (c[0], c[6])
    return {k: sorted(v) for k, v in d.items()}, info

intr = lambda ex: {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}
hx_ex, hx_info = exones_gff(RS + "Dmaidis_helixer.gff3", r"Parent=([^;]+)", "exon")
hx_cds, _ = exones_gff(RS + "Dmaidis_helixer.gff3", r"Parent=([^;]+)", "CDS")
st_ex, st_info = exones_gff(RS + "stringtie_denovo.gtf", r'transcript_id "([^"]+)"', "exon")
st_por_cr = defaultdict(list)
for t, ex in st_ex.items():
    st_por_cr[st_info[t][0]].append((ex[0][0], ex[-1][1], t))

def leer_fa(path, ids):
    out, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]; cur = cur if cur in ids else None
            if cur: out[cur] = []
        elif cur: out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}

filas = []
for tg, r in F.iterrows():
    mrna = r.dmai_representativo.replace("Dmai|", "") + ".1"
    ex, cds = hx_ex[mrna], hx_cds[mrna]
    cr, hebra = hx_info[mrna]
    hi = intr(cds); ini, fin = ex[0][0], ex[-1][1]
    cds_len = sum(b - a + 1 for a, b in cds)
    mejor, union = (-1, -1.0, None), set()
    for a, b, t in st_por_cr[cr]:
        if b < ini or a > fin or st_info[t][1] not in (hebra, "."): continue
        sx = st_ex[t]
        union |= intr(sx)
        conf = len(hi & intr(sx))
        cub = sum(max(0, min(b2, e2) - max(a2, s2) + 1) for a2, b2 in cds for s2, e2 in sx) / cds_len
        if (conf, cub) > mejor[:2]: mejor = (conf, cub, t)
    _, cub, t = mejor
    conf = len(hi & union)
    estado = "sin_stringtie" if t is None else \
             "confirmado" if conf == len(hi) and cub >= 0.95 and len(hi & intr(st_ex[t])) == len(hi) else \
             "uniones_ok_secuencia_parcial" if conf == len(hi) else "discrepancia"
    lo, hi_c = min(a for a, b in cds), max(b for a, b in cds)
    extra = len({i for i in intr(st_ex[t]) if lo < i[0] < hi_c} - hi) if t else 0   # intrones de StringTie dentro del CDS que Helixer no tiene
    filas.append(dict(rank=int(r["rank"]), tcas_geneid=tg, dmel=r.dmel, helixer=mrna, stringtie=t or "",
                      intrones_helixer=len(hi), intrones_confirmados=max(conf, 0),
                      intrones_solo_stringtie=extra, cds_cubierto=round(max(cub, 0), 3), estado=estado,
                      id_GITV01=r.id_GITV01, confianza=r.confianza_gen, riesgo_paralogo=r.riesgo_paralogo,
                      paralogos_prot=r.paralogos_prot))
V = pd.DataFrame(filas)
V.to_csv("ranking/candidatos_20_verificacion.tsv", sep="\t", index=False)

st_seq = leer_fa(RS + "stringtie_tx.fa", set(V.stringtie))
hx_seq = leer_fa(RS + "helixer_mrna.fa", set(V.helixer))
with open("candidatos_20_verificados.fna", "w") as f:
    for v in V.itertuples():
        if v.estado == "confirmado":
            s, fuente = st_seq[v.stringtie], f"StringTie {v.stringtie} (coincide con Helixer)"
        elif v.stringtie and v.cds_cubierto >= 0.90:
            s, fuente = st_seq[v.stringtie], f"StringTie {v.stringtie} (REVISAR: difiere de Helixer, {v.estado})"
        else:
            s, fuente = hx_seq[v.helixer], f"Helixer {v.helixer} (NO verificado: {v.estado})"
        f.write(f">cand{v.rank:04d}_{v.dmel} tcas={v.tcas_geneid} {fuente} intrones={v.intrones_confirmados}/{v.intrones_helixer}"
                f" cds_cubierto={v.cds_cubierto}\n")
        for i in range(0, len(s), 70): f.write(s[i:i + 70] + "\n")

print(f"excluidos antes de completar 20 (rank <= {int(F['rank'].max())}): {len(excluidos)}")
for g, e in excluidos.iterrows():
    mot = []
    if e.paralogo_resuelto == "ambas_con_lecturas": mot.append("2 copias con lecturas")
    if e.max_ident_paralogo >= 70: mot.append(f"paralogo prot {e.max_ident_paralogo:.0f}%")
    if e.loci_genoma != 1: mot.append(f"loci={e.loci_genoma}")
    if e.confianza_gen == "baja": mot.append("confianza baja")
    print(f"  #{int(e['rank'])} {e.dmel}: {', '.join(mot)}")
print(V.drop(columns=["tcas_geneid", "helixer"]).to_string(index=False, max_colwidth=22))
print(V.estado.value_counts().to_dict())
