#!/usr/bin/env python3
"""
expresion_matrices.py — matrices de expresion por gen Helixer, separadas por BioProject.

Fuente: fragmentos por gen y corrida contados sobre los BAM de HISAT2 al genoma (sbatch r4_conteos_por_corrida:
read1 primario, MAPQ >= 10, cuenta si algun bloque cae en un exon del gen). Se dejo salmon-sobre-CDS porque daba TPM = 0
en algun estadio a ~2,000 genes con >= 1,000 lecturas de HISAT2 (artefacto del indice solo-CDS).
TPM = (fragmentos / largo exonico en kb) reescalado a suma 1e6 por muestra.

Salidas en expresion/:
  metadatos_corridas.tsv   por corrida: BioProject, BioSample, estadio, sexo, tratamiento, tiempo, temperatura, fragmentos
  tpm_corridas.tsv         TPM por gen x corrida (19 corridas; encabezado = run, cruzar con metadatos_corridas.tsv)
  estadios_318849.tsv      PRJNA318849 (ninfa, hembra, macho; 1 corrida c/u): TPM por estadio, TPM_minimo, CV_estadios,
                           sesgo_sexo = log2((TPM_hembra + 1) / (TPM_macho + 1))
  condiciones_631706.tsv   PRJNA631706: las 4 corridas de cada condicion son replicas TECNICAS de una misma BioSample ->
                           se colapsan sumando fragmentos antes de calcular TPM. CV_condiciones entre las 4 condiciones y
                           log2FC_MRFV_4h / _7d = log2((TPM_MRFV + 1) / (TPM_naive + 1)). Fold changes DESCRIPTIVOS:
                           sin test ni umbral (con replicas tecnicas la varianza esta subestimada).
CV = desvio estandar muestral / media (genes con media 0 -> NaN).
"""
import glob
from pathlib import Path
import numpy as np
import pandas as pd

CT = Path("rnaseq_insumos/conteos")
OUT = Path("expresion"); OUT.mkdir(exist_ok=True)

# BioSample por corrida: ENA filereport (2026-10-28); condiciones: titulos de BioSample/SRA
META = {
    "SRR5010611": ("PRJNA318849", "SAMN06009808", "ninfa", "mixto", "ninguno", ""),
    "SRR5010612": ("PRJNA318849", "SAMN06009807", "adulto", "hembra", "ninguno", ""),
    "SRR5010613": ("PRJNA318849", "SAMN06009806", "adulto", "macho", "ninguno", ""),
}
BS_631706 = {"SAMN14890446": ("naive", "4h"), "SAMN14890505": ("MRFV", "4h"),
             "SAMN14890635": ("naive", "7d"), "SAMN14890756": ("MRFV", "7d")}
RUN_BS = {"SRR11822339": "SAMN14890446", "SRR11822340": "SAMN14890446", "SRR11822347": "SAMN14890446",
          "SRR11822348": "SAMN14890446", "SRR11822335": "SAMN14890505", "SRR11822336": "SAMN14890505",
          "SRR11822337": "SAMN14890505", "SRR11822338": "SAMN14890505", "SRR11822333": "SAMN14890635",
          "SRR11822334": "SAMN14890635", "SRR11822345": "SAMN14890635", "SRR11822346": "SAMN14890635",
          "SRR11822341": "SAMN14890756", "SRR11822342": "SAMN14890756", "SRR11822343": "SAMN14890756",
          "SRR11822344": "SAMN14890756"}
for r, bs in RUN_BS.items():
    META[r] = ("PRJNA631706", bs, "adulto", "pool_hembra_macho", *BS_631706[bs])

bed = pd.read_csv(CT / "genes.bed12", sep="\t", header=None, index_col=3)
largo_kb = bed[10].map(lambda s: sum(int(x) for x in s.rstrip(",").split(","))) / 1000
N = pd.DataFrame({Path(f).stem: pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
                  for f in sorted(glob.glob(str(CT / "SRR*.tsv")))})
N.index = "Dmai|" + N.index
N.index.name = "gen_dmai"
largo_kb.index = "Dmai|" + largo_kb.index
largo_kb = largo_kb.reindex(N.index)
def tpm(cnt):
    rpk = cnt.div(largo_kb, axis=0)
    return rpk / rpk.sum() * 1e6

M = pd.DataFrame.from_dict(META, orient="index", columns=["bioproject", "biosample", "estadio", "sexo",
                                                          "tratamiento", "tiempo"])
M["temperatura"] = np.where(M.bioproject == "PRJNA631706", "30C", "")
M["fragmentos_en_genes"] = N.sum().astype(int)
M.index.name = "run"
M = M.sort_values(["bioproject", "tratamiento", "tiempo"])
M.to_csv(OUT / "metadatos_corridas.tsv", sep="\t")
TPM = tpm(N)[M.index]
TPM.to_csv(OUT / "tpm_corridas.tsv", sep="\t", float_format="%.4g")
cv = lambda d: (d.std(axis=1) / d.mean(axis=1)).where(d.mean(axis=1) > 0)

# ---- PRJNA318849: estadios ----
E = TPM[["SRR5010611", "SRR5010612", "SRR5010613"]].set_axis(["tpm_ninfa", "tpm_hembra", "tpm_macho"], axis=1)
E["TPM_minimo"] = E.min(axis=1)
E["CV_estadios"] = cv(E[["tpm_ninfa", "tpm_hembra", "tpm_macho"]])
E["sesgo_sexo"] = np.log2((E.tpm_hembra + 1) / (E.tpm_macho + 1))
E.to_csv(OUT / "estadios_318849.tsv", sep="\t", float_format="%.4g")

# ---- PRJNA631706: colapsar replicas tecnicas por BioSample (suma de fragmentos) ----
C = tpm(pd.DataFrame({f"{t}_{h}": N[[r for r, b in RUN_BS.items() if b == bs]].sum(axis=1)
                      for bs, (t, h) in BS_631706.items()}))[["naive_4h", "MRFV_4h", "naive_7d", "MRFV_7d"]]
C = C.add_prefix("tpm_")
C["CV_condiciones"] = cv(C)
C["log2FC_MRFV_4h"] = np.log2((C.tpm_MRFV_4h + 1) / (C.tpm_naive_4h + 1))
C["log2FC_MRFV_7d"] = np.log2((C.tpm_MRFV_7d + 1) / (C.tpm_naive_7d + 1))
C.to_csv(OUT / "condiciones_631706.tsv", sep="\t", float_format="%.4g")

print(M.groupby(["bioproject", "estadio", "sexo", "tratamiento", "tiempo"]).size().to_string())
print("fragmentos en genes por corrida (M):", (M.fragmentos_en_genes / 1e6).round(1).to_dict())
print(f"genes: {len(TPM)}; TPM_minimo = 0 en estadios: {(E.TPM_minimo == 0).sum()}; < 1: {(E.TPM_minimo < 1).sum()}")
print("CV_estadios mediana:", round(E.CV_estadios.median(), 3), "| CV_condiciones mediana:", round(C.CV_condiciones.median(), 3))
print("sesgo_sexo |log2| > 2:", int((E.sesgo_sexo.abs() > 2).sum()))
