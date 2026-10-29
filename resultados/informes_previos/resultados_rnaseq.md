# RNA-seq como evidencia independiente — resultados

Fecha: 2026-10-28. Scripts: `scripts/puente_gitv01.py`, `scripts/rnaseq_integracion.py`, `scripts/verificar_candidatos.py`.
Jobs de SDumont: `r1` (HISAT2 + bedtools + StringTie), `r2`/`r4` (minimap2), `r3` (DIAMOND Dmai vs Dmai). DIAMOND del
puente: cluster QB. El ranking previo quedó en `ranking/ranking_completo.pre_rnaseq.tsv`; las columnas nuevas están en
`ranking/ranking_completo.tsv`.

## Correcciones al brief (verificadas)

- El paper del transcriptoma y de MRFV **no es "Martin et al. 2021"**: es **Xu, Willman, Todd, Kim, Redinbaugh & Stewart
  (2021) *Microbiol. Spectrum* 9:e00612-21**, doi 10.1128/Spectrum.00612-21.
- Sus lecturas (1.7 mil millones, **PRJNA579843**) **no** son las 19 corridas que se usaron acá. Las corridas usadas son:
  - **PRJNA631706** (16): MRFV vs. no expuesto × 4 h / 7 d, **a 30 °C**. Son 4 corridas por condición de **una misma
    BioSample**, es decir réplicas técnicas.
  - **PRJNA318849** (3): **ninfas**, machos y hembras.

  O sea: no es solo RNA-seq de adultos, y no es el experimento de Xu et al.

## 1. Alineamiento al genoma (Tarea A)

HISAT2 contra `GCA_975739785.1`, por corrida:

| Grupo | Corridas | Al genoma | Salmon contra CDS de Helixer |
|---|---|---|---|
| PRJNA631706 (MRFV/naive) | 16 | 88.7–94.2% (mediana 93.1%) | 56–76% |
| PRJNA318849 (estadios) | 3 | 92.4–95.1% | **35–51%** |

Detalle por corrida en `rnaseq_insumos/alineamiento_por_corrida.tsv`.

**Diagnóstico: los datos están bien** (≥ 88.7% al genoma en todas las corridas, incluidas las tres que daban 35–51% en
Salmon). La pérdida en Salmon es del lado de la anotación y del índice, y tiene dos componentes:

1. El índice de Salmon era **solo de CDS** (`gffread -x`): las lecturas de UTR e intrones no tenían dónde caer.
2. **~19% de los alineamientos al genoma cae fuera de cualquier gen de Helixer.** De 938 M alineamientos, 764 M (81.4%) caen
   dentro de los genes. Es transcripción no anotada: genes que Helixer no predijo, ncRNA o TEs activos.

Es decir, aplica una versión moderada del caso "preocupante": Helixer no cubre una parte de lo transcripto. No afecta a los
candidatos (ver 2), pero sí a la completitud del proteoma (genes reales que no están). Que PRJNA318849 caiga más en Salmon
es compatible con que tenga más transcripción fuera de CDS (ninfas, otros estadios), pero no lo separé por corrida.

## 2. Validación de Helixer con lecturas (Tarea B)

`rnaseq_soporte` = **≥ 50% de las bases exónicas cubiertas y ≥ 10 lecturas** sobre el gen. Con ~390 M de pares, 10 lecturas
descarta lecturas sueltas (ADN contaminante, multimapeo en TEs) sin castigar a genes poco expresados; el 50% exige que la
transcripción recorra el modelo.

- **20,439 de 22,199 genes de Helixer (92.1%)** tienen soporte.
- **Top 100: 100/100 con soporte.** La cobertura mínima en el top 100 es 0.65 (mediana 1.0) y la mínima de lecturas es 9,947.
  Ningún candidato del top 100 parece un modelo falso.

## 3. Puente con GITV01 (Tarea G) → `mapeo_helixer_GITV01.tsv`

- **RBH de proteína** (DIAMOND, `--id 80 --query-cover 60 --subject-cover 60`, proteínas GITV01 de TransDecoder):
  **5,622 de 22,199** (25.3%). Top 100 con `id_GITV01`: **42**.
- Es bajo porque solo 23,848 de los 56,117 transcriptos tienen ORF ≥ 100 aa, y muchas están fragmentadas y no pasan la
  cobertura 60/60. Por eso agregué un **puente por locus** (minimap2 `-c -x splice:hq`, identidad ≥ 0.90, ≥ 50% del
  transcripto alineado, ≥ 100 pb sobre exones del gen): **10,407 genes**. RBH o locus: **11,236 genes (50.6%)**.
- La columna `otros_GITV01` del mapeo lista los demás transcriptos de GITV01 cuyo mejor hit es ese gen (isoformas o alelos
  del ensamblado). Sirve para buscar en la literatura.

## 4. Confianza del gen (Tarea H)

`soporte_GITV01` = RBH **o** locus.

| Confianza | Top 100 | Todo el ranking (7,881) |
|---|---|---|
| alta (lecturas + GITV01) | **68** | 5,828 |
| media (lecturas, sin GITV01) | **32** | 2,025 |
| media_revisar (GITV01, sin lecturas) | 0 | 2 |
| baja (ninguna) | **0** | 26 |

## 5. Top 100 con confianza baja

**Ninguno.** Los 26 genes con confianza baja de todo el ranking quedan fuera del top 100.

## 6. Flags de parálogo (Tarea C)

De los 896 genes con 2 copias en *D. maidis*:

| `paralogo_resuelto` | Genes |
|---|---|
| ambas_con_lecturas | **886** |
| solo_una | 9 |
| ninguna | 1 |

**El RNA-seq casi no rescata copias falsas:** en el 99% de los casos las dos copias tienen lecturas propias, así que son genes
reales distintos y el riesgo de compensación es real. En el top 100 hay **4** con el flag, y los 4 tienen ambas copias con
lecturas: #38 CG7483, #40 Vha16-2/3, #71 Rho1, #83 Vha100-2..5.

## 7. Loci y parálogos antiguos (Tarea F + complemento)

- **minimap2** (`-c -x splice:hq`, mRNA de Helixer contra el genoma, identidad ≥ 0.90, cobertura ≥ 0.5): **los 100 del top 100
  mapean a 1 solo locus**.
- **Esto no detecta parálogos antiguos:** a ≥ 90% de identidad nucleotídica solo aparecen copias recientes. Para cubrir lo
  que pedía el brief (parálogos que divergieron y quedaron en otro ortogrupo) agregué **DIAMOND *D. maidis* contra *D.
  maidis*** (proteína, ≥ 40% de identidad y ≥ 50% de cobertura de ambas). Columnas `paralogos_prot` y `max_ident_paralogo`.
  - Top 100 con parálogos proteicos: 0 → 72 genes; 1 → 9; 2 → 3; ≥ 3 → 16.
  - Con un parálogo ≥ 70% de identidad (riesgo de compensación **y** de co-silenciamiento): #6 eIF4A (73%), #18 Act5C (99%),
    #38 CG7483 (98.5%), #40 Vha16 (83%), #62 Arf1 (79%), #71 Rho1 (95%).
  - Muchos de los ≥ 3 son familias que no son funcionalmente redundantes a ~40% (las ATPasas *Rpt1–6* entre sí, por
    ejemplo); por eso no se usa como filtro duro salvo a ≥ 70%.

### Candidatos finales y secuencia verificada (Tarea D) → `candidatos_20_verificados.fna`

**Finales** = primeros 20 del ranking con confianza ≠ baja, 1 locus, sin dos copias expresadas y sin parálogo ≥ 70%.
Excluidos en el camino: #6 eIF4A y #18 Act5C.

StringTie de novo (sin anotación guía) contra Helixer, comparando **solo los intrones del CDS**. Un intrón se considera
confirmado si alguna lectura salta exactamente esa unión.

- **Confirmados (9):** todos los intrones del CDS confirmados y transcripto que cubre el CDS. Rpt4, Rpt3, snRNP-U1-70K,
  Prosα6, Rpt5, Rab1, VhaSFD, RpLP0, Vha14-1.
- **Discrepancia (11):** Rpt1 (7/9 intrones), Rpn11 (12/14), Rpt6 (6/7), Hel25E (7/11), Prosα2 (4/6), Prosα5 (5/8), Rpt2
  (7/10), Polr2E (3/5), CkIIα (6/8), hay (10/16), CG1703 (10/13). Helixer predice intrones que ninguna lectura salta, y en
  varios casos StringTie tiene intrones dentro del CDS que Helixer no. **Los límites de exones de Helixer están corridos en
  la mitad de los finales**, que es justo el riesgo que el brief quería medir.
- En el FASTA va el **transcripto de StringTie** (respaldado por lecturas) cuando cubre ≥ 90% del CDS: 18 de 20. Los que
  discrepan quedan marcados `REVISAR` en el encabezado. **Rpn11** (CDS cubierto 0.53) y **Prosα5** (0.87) van con la
  secuencia de Helixer marcada `NO verificado`.
- Detalle en `ranking/candidatos_20_verificacion.tsv`.

## 8. Competencia vectorial (Tarea I) → `ranking/candidatos_competencia_vectorial.tsv`

- **Fuente:** Tablas 3 y 4 de Xu et al. 2021 (DESeq2, p ajustado ≤ 0.1). Son **83 transcriptos con anotación** (4 h: 54 up y
  3 down; 7 d: 26 up), en `mrfv/mrfv_DE.tsv`. El paper reporta 72 + 67. Los que no tienen anotación no están en las tablas, y
  la Tabla S7 da fold changes pero **no p ajustados**, así que no los recuperé ni los inventé.
- **Mapeo a Helixer:** locus 7, RBH 10, mejor hit directo de proteína 10 → **27 de 83**; **24 en el ranking**. Lo verifiqué
  con la SOD (contig 4722 → gen de Helixer ortólogo de *Sod1* de *Drosophila*).
  - 38 de los 83 **no alinean al genoma**, mientras que en todo GITV01 alinea el 96%. No encontré la causa: la numeración de
    contigs está bien, como muestra la SOD. Queda como limitación.
- Lista aparte (sin mezclar con el ranking principal), con 9 de los 24 que tienen letalidad ≥ p75 o estabilidad ≥ p75. Los
  más interesantes:
  - **TER94** (#68, up 7 d, letalidad 0.86).
  - **VhaPPA1-1** (#99, up 4 h, 0.92).
  - **RpL6** (#240, up 7 d, letalidad 0.99, estabilidad 0.77).
  - **Slu7** (#263, up 7 d).
  - **Sod1** (#1137, up 4 h, estabilidad 0.80, letalidad baja: es candidato de competencia, no de letalidad).

## 9. Estabilidad de expresión (Tarea E, informativa, no entra al score)

`estabilidad_expresion` = 1 − percentil del rango log2 entre las medias de las 4 condiciones (MRFV/naive × 4 h/7 d).

| | Estabilidad |
|---|---|
| PSMB5 | 0.71 |
| V-ATPasa B | 0.74 |
| V-ATPasa D | 0.66 |
| Mediana del ranking | 0.44 |
| Mediana del top 100 | 0.68 |

Los controles salen por encima de la mediana general, pero **igual que el resto del top 100**. Correlación con el score:
Spearman 0.13. Discrimina algo, pero no separa a los controles del resto de los buenos candidatos. Además, con réplicas
técnicas no hay varianza biológica: el rango entre condiciones mezcla efecto real y ruido de muestreo. **No conviene sumarla
al score con estos datos.**

## Limitaciones

- La ausencia de cobertura no prueba que un gen sea falso: el RNA-seq es de adultos (y ninfas, en 1 corrida) de cuerpo entero.
- La ausencia de contraparte en GITV01 tampoco lo prueba (GITV01 tiene BUSCO C 60%). Por eso las dos evidencias se combinan.
- Las 16 corridas de MRFV son 4 réplicas técnicas por condición, sin réplicas biológicas.
- ~19% de las lecturas cae fuera de genes de Helixer: el proteoma no es completo.
- Helixer tiene límites de exones discordantes en 11 de 20 candidatos finales. Para el dsRNA usar la secuencia de StringTie y
  revisar los marcados.
- La lista de MRFV está incompleta (solo transcriptos anotados) y 38 de 83 no alinean al genoma.
- **Lo que este RNA-seq no puede dar:** expresión por tejido, en particular intestino medio (relevante para la vía oral).
  Requiere RNA-seq de tejidos disecados. Sigue siendo el experimento pendiente de mayor retorno.
