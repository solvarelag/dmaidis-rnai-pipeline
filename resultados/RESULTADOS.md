# Resultados — Priorización de blancos de RNAi en *Dalbulus maidis*

Versión: 2026-10-28, después de la regeneración completa (tarea 4 de `brief_metamorfosis_y_cierre.md`).
Este documento **reemplaza** a los informes parciales de `informes_previos/`: `resultados_rnaseq.md` y
`resultados_score.md` tienen cifras previas a las correcciones de esta versión.
El método está en `METODOLOGIA.md`; qué contiene cada archivo, en `LEEME.md`.

---

## 1. Insumos genómicos y anotación

| Ítem | Resultado |
|---|---|
| Genoma | GCA_975739785.1 (pseudohaplotipo principal, 972 Mb). GCA_975739835.1 es el alternativo y no se usa. |
| BUSCO del genoma (insecta) | C 99.6%, D 1.8% |
| Anotación oficial de *D. maidis* | No existe |
| Proteoma usado | **Helixer** (v0.3.6), 22,199 genes, BUSCO insecta **C 93.0%, D 1.6%** |
| Alternativa descartada | TransDecoder sobre GITV01: C 60%, D 18.6% (no pasó el chequeo 1) |
| *Macrosteles quadrilineatus* | Inflado por haplotipos (D 38% aun después de CD-HIT): no se usa para contar parálogos |

## 2. Ortología (OrthoFinder 2.5.5, 10 proteomas)

- Pasada principal (*Tribolium* RefSeq): 16,941 ortogrupos; 4,630 con las 10 especies; 6,528 pares 1:1 Dmai–Tcas.
  *D. maidis* tiene 89.6% de sus genes en ortogrupos.
- **Análisis separado con *Tribolium* OGS3** (IDs TC nativos, sin mapeo TC→GeneID; isoforma más larga, 16,593 genes;
  79.3% en ortogrupos). OGS3 no se sumó como especie número 11: dos proteomas casi idénticos de *Tribolium* distorsionan
  los ortogrupos.

| | RefSeq (oficial) | OGS3 (separado) |
|---|---|---|
| Universo con letalidad y ortólogo Dmai | 8,272 | 9,147 |
| Tras el filtro de copias | **7,881** | **8,651** |
| Copias 1 / 2 / ≥ 3 | 6,985 / 896 / 391 | 7,624 / 1,027 / 496 |
| iB letales (dali11 ≥ 80) en el ranking | 961 | 1,015 (69 solo en OGS3, 15 solo en RefSeq) |
| Top 100 compartido (por gen Dmai) | — | **96/100**; Spearman 0.992 en 7,171 genes comunes |
| PSMB5 / V-ATPasa D / V-ATPasa B | #25 / #125 / #216 | #27 / #128 / #219 |

- Solo en el top 100 de OGS3: **Rab6** (#26; #299 en RefSeq), **Cam** (#47, 2 copias; ausente en RefSeq), **atms** (#60,
  ausente) y **salsa** (#74; #195).
- Conclusión: el mapeo de IDs perdía ~54 iB letales netos y 770 genes del universo, pero casi no cambia la cabeza del
  ranking. El ranking oficial sigue siendo el de RefSeq.

## 3. Letalidad (iBeetle)

- Letales con dali11 ≥ 80: 1,499 iB → 1,142 con GeneID → 980 con ortólogo en Dmai.
  - Copias: 1 → 88.8%; 2 → 9.1%; ≥ 3 → 2.1%.
- Imputación de faltantes en el universo de 7,881:
  - dali22 4,658
  - dapi11 400
  - dali11 171
  - 4,707 genes tienen al menos una imputación.
- **dali22 se imputa como max(dali11, media)**, decisión de la usuaria. La mortalidad es acumulada: dali22 ≥ dali11 en el
  99.9% de las iB con ambas medidas. La media sola daba valores imposibles y castigaba a los genes que matan rápido.
- **Sensibilidad sin dali22** (letalidad = 0.6·dali11 + 0.4·dapi11): top 100 compartido **92/100**; Spearman 0.960.
  - Salen: Vha14-1 (#48 → 205; dali11 = 0, dali22 = 100), CG7483, Caf1-55, pit, Rpn12, Vha44, ATPsynδ y ND-ACP.
  - Entran: pont, Vha13, CkIIβ, LysRS, eIF4E, eEF2, Nfs1, mEFTu1.
  - Archivo: `tablas/ranking/sensibilidad_sin_dali22.tsv`.

## 4. Score y ranking oficial (A)

Score = 0.35 letalidad + 0.30 firma GO + 0.20 expresión constitutiva + 0.15 conservación (rango percentil).

**Top 30:**

| # | Gen (Dmel) | Proteasoma | Flag metam. | Confianza | # en B |
|---|---|---|---|---|---|
| 1 | Rpt1 | sí | | alta | 2 |
| 2 | Rpt5 | sí | | media | 4 |
| 3 | Vha16-2/3 | | | media | 6 |
| 4 | Rpt2 | sí | | alta | 8 |
| 5 | EloC | | ✔ | alta | 7 |
| 6 | Rpt6 | sí | | media | 9 |
| 7 | mts | | ✔ | alta | 10 |
| 8 | Rpn11 | sí | | alta | 14 |
| 9 | Rpn6 | sí | | media | 24 |
| 10 | Hsp83 | | | alta | 19 |
| 11 | snRNP-U1-70K | | | alta | 16 |
| 12 | Sar1 | | ✔ | alta | 5 |
| 13 | Act5C | | | alta | 23 |
| 14 | Arf1 | | ✔ | media | 15 |
| 15 | Rpt4 | sí | | alta | 30 |
| 16 | eIF4A | | | alta | 25 |
| 17 | Hel25E | | | alta | 31 |
| 18 | CG1703 | | | alta | 26 |
| 19 | Vha26 | | | media | 41 |
| 20 | hay | | | alta | 39 |
| 21 | Rpt3 | sí | | alta | 17 |
| 22 | TER94 | parcial | ✔ | alta | 42 |
| 23 | gskt/sgg | parcial | ✔ | alta | 51 |
| 24 | VhaSFD | | | alta | 54 |
| 25 | Prosβ5 (PSMB5) | sí | | alta | 34 |
| 26 | Vha68 | | ✔ | alta | 60 |
| 27 | crn | | | alta | 47 |
| 28 | lwr | | | media | 3 |
| 29 | Rho1 | | ✔ | media | 22 |
| 30 | Pop2 | | | media | 21 |

- **Controles positivos:**
  - PSMB5 (blanco de ledprona): #25
  - Rpt3: #21
  - V-ATPasa D (eficaz en *P. maidis*): #125
  - V-ATPasa B: #216
- Al inicio de la sesión eran #24, #5, #220 y #149. Los cambios vienen de la imputación de dali22 y de la expresión por
  estadios.
- **Resto del set de entrenamiento:**
  - Hsc70-3 #65
  - dre4 #93
  - Pp1α-96A #131
  - RpII140 #146
  - ncm #153
  - Rpn7 #229 (conservación 0.25)
  - Srp54k #230
  - αSnap #527
  - Rop #650 (riesgo_paralogo)
  - gw #1111
  - cact #1152
  - shi #1286
  - Inr-a #2736
- **Sesgo temático:** 22 del top 100 son subunidades del proteasoma (GO:0000502) y 7 más son parciales (catabolismo
  proteasomal).
- **Reparto del score en el top 50:** letalidad 0.341 de 0.35, GO 0.294 de 0.30, expresión 0.187 de 0.20,
  conservación 0.142 de 0.15. Ningún criterio domina solo, porque los cuatro saturan juntos.

## 5. Firma GO y circularidad (tareas 0, 2 y 6)

- Firma: **77 términos**, 17 genes de entrenamiento (16 en el fondo 1:1; Rop no está).
  - Detalle de cada término y su origen: `firma_GO_entrenamiento.txt`.
  - Ficha de cada gen de entrenamiento: `set_entrenamiento_GO.txt`.
- **Proteasomales en el entrenamiento: 3/17** (PSMB5, Rpn7, Rpt3); Hsc70-3 es parcial.
  - Términos de la firma que provienen **solo** de genes proteasomales: **6/77 (7.8%)**.
  - Los genes que más términos aportan son V-ATPasa B (33), V-ATPasa D (31), Rpt3 (29), Hsc70-3 (28) y shi (25).
- **Auditoría:**
  - 0/77 términos metamórficos.
  - **33/77 genéricos** (profundidad ≤ 3 en el DAG o anotados en > 25% del universo).
  - Profundidad: mínima 1, mediana 4, máxima 10.
- **Ranking B (sin GO; letalidad 50%, expresión 28.6%, conservación 21.4%):**
  - Top 100 A ∩ B: **63/100**; Spearman A vs. B: 0.898.
  - **Proteasomales en el top 100: A 22, B 15.**
  - Controles en B: PSMB5 #34, Rpt3 #17, V-ATPasa D #65, V-ATPasa B #184.
  - 37 genes entran al top 100 solo en B: Sec61α (#1 de B), Sec23/Sec24CD y el coatómero (γ, δ y β'COP), factores de
    traducción eIF3 (a, c, d, e), eRF1, splicing (Prp31, Sf3b2, Sf3b3, Bx42), ATP sintasa (γ, F), Gs1/2, AsnRS, Sams,
    Acly, Eno y Nop56, entre otros.
  - **Lectura** (criterio del brief): el proteasoma **sigue en el top 100 sin GO** (15 genes; Rpt3 sube a #17 y PSMB5 queda
    #34). La señal proteasomal no es solo un efecto de la firma. A sirve como ranking confirmatorio y B como exploratorio:
    B destapa secreción/COPII–COPI, traducción y splicing, que la firma GO dejaba más abajo.

## 6. Sesgo de metamorfosis (tareas 1 y 3)

- **Conjunto metamórfico:** 146 términos estrictos y 206 amplios, expandidos por DAG desde 7 semillas y por nombre. Se
  excluyeron muda, quitina, cutícula, ecdisis y embrión, y lo que aplica igual a un hemimetábolo. Todo es auditable en
  `metamorfosis/terminos_metamorficos.tsv` y `metamorfosis/exclusiones_aplicadas.tsv`.
- **Marcados:**
  - Universo: 625/7,881 (7.9%)
  - OGS3: 646/8,651 (7.5%)
  - **Top 100: 15**
  - **20 finales: 5** (EloC, mts, Sar1, TER94, gskt/sgg)
  - Ranking B: 18 en el top 100
- Criterio del brief: entre 5 y 20 → **inspección caso por caso**. Observación para esa inspección: los marcados del top
  100 son genes esenciales pleiotrópicos (PP2A/mts, Rho1, Arf1, Sar1, sqh, TER94, sgg, Cdk8). Se marcan porque en
  *Drosophila* sus mutantes afectan discos imaginales o el desarrollo larva-pupa, y no son genes de la red metamórfica
  propiamente dicha. **La firma GO no tiene ningún término metamórfico:** el sesgo no entra por el 30% de GO.
- **`letal_ambos_regimenes`** (alto en inyección larval y pupal; corte p75 de valores observados: dali ≥ 70, dapi11 ≥ 40):
  - 956 genes (12.1%); en el top 100, **93**.
  - Spearman con el score: 0.448.
  - Sensibilidad: p50 → 96 en el top 100; p90 → 76.
  - **No mide lo mismo que flag_metamorfosis.** Entre los letal_ambos = False, 7.5% están marcados; entre los True,
    10.9%.
  - El filtro deja 93 del top 100 original. Es diagnóstico: no se aplica al ranking oficial.

## 7. RNA-seq (dos BioProjects, tratados por separado)

**Corridas usadas** (metadatos en `tablas/expresion/metadatos_corridas.tsv`):

| BioProject | Contenido | Uso |
|---|---|---|
| PRJNA318849 | 3 corridas: ninfa (sexo mixto), hembra adulta, macho adulto | **Expresión constitutiva** (20% del score) y estabilidad entre estadios |
| PRJNA631706 | 16 corridas = 4 BioSamples × 4 réplicas **técnicas**: naive/MRFV × 4 h/7 d, 30 °C, adultos | Solo **descriptiva**: CV entre condiciones y log2FC; sin tests |

- **PRJNA579843** (Xu et al. 2021, 25 °C, la fuente de las Tablas 3/4 de MRFV) **no tiene lecturas públicas**: 0 corridas
  en SRA y en ENA. El registro TSA de GITV01 está bajo PRJNA631706 y lista 12 BioSamples: las 4 de 30 °C son las de las 16
  corridas usadas, y las de 25 °C no tienen corridas. El zip `Dalbulus_maidis_transcriptome_GITV01.zip` es el ensamblado
  (56,117 transcriptos), no lecturas.
- **Tarea A:** HISAT2 alinea **88.7–95.1%** de cada corrida al genoma: los datos están bien. La baja tasa de Salmon
  (35–76%) venía del índice solo de CDS. El "~19% de lecturas fuera de genes" del informe previo se midió por span, sin
  `-split`, y es una cota aproximada.
- **Expresión (cambio de método):** conteos de fragmentos exónicos con MAPQ ≥ 10 sobre los BAM de HISAT2, en lugar de
  Salmon.
  - 7,290 genes tienen TPM = 0 en algún estadio. Se verificó que los ceros son reales: las lecturas del span eran
    multimapeo o intrones de lecturas de otro gen.
  - En el universo hay 333 genes con TPM_minimo = 0; en el top 100, 0. La mediana del TPM mínimo en el top 100 es 104.
  - Medianas de CV: entre estadios 0.70 en todos los genes; entre condiciones 0.197.
  - |log2(hembra/macho)| > 2 en 1,620 genes del universo; en el top 100, 19 tienen |sesgo| > 1.
- **Tarea B (corregida con `-split`):** soporte (≥ 50% de bases exónicas cubiertas y ≥ 10 fragmentos) en **17,023 de
  22,199 genes (76.7%)**. Antes daba 92.1%, inflado porque los intrones contaban como cobertura. **Top 100: 100/100.**
- **Tarea G, puente con GITV01:**
  - RBH de proteína: 5,622 genes
  - Por locus (minimap2): 10,407
  - Cualquiera de los dos: 11,236
  - Top 100: 46 con RBH, 70 con algún soporte
- **Tarea H, confianza del top 100:** **70 alta, 30 media, 0 baja**. En todo el ranking: 5,818 alta, 2,000 media,
  12 media_revisar, 51 baja.
- **Tarea C, genes con 2 copias (896):** ambas con lecturas 857, solo una 35, ninguna 4. En el top 100 hay 3 y los 3 tienen
  ambas copias expresadas.
- **Tarea F:**
  - minimap2 con todos los mRNA: los 100 del top 100 tienen **1 locus**.
  - Parálogos proteicos (DIAMOND Dmai vs. Dmai, ≥ 40% de identidad, ≥ 50% de cobertura) en el top 100: 0 → 73 genes;
    1 → 10; 2 → 3; ≥ 3 → 14.
- **Tarea E, estabilidad (CV; ninguna entra al score):**

| | CV estadios (percentil) | CV condiciones (percentil) |
|---|---|---|
| PSMB5 | 0.114 (0.04) | 0.054 (0.10) |
| V-ATPasa D | 0.268 (0.19) | 0.080 (0.20) |
| V-ATPasa B | 0.330 (0.25) | 0.045 (0.07) |
| Mediana del ranking / del top 100 | 0.560 / 0.252 | 0.150 / 0.072 |

  Los controles caen en el cuarto más estable, igual que el top 100 en general. Spearman entre CV_condiciones y el score:
  −0.342; entre los dos CV: 0.439. estabilidad_estadios es el mismo CV_estadios que ya está dentro del subscore de
  expresión.

## 8. Candidatos finales y secuencias (tarea D)

**Criterio:** los 20 primeros del ranking A que cumplen todo esto:

- confianza distinta de baja
- 1 locus
- no tienen dos copias expresadas
- no tienen un parálogo proteico con ≥ 70% de identidad

Antes de completar los 20 se excluyeron #3 Vha16-2/3 (2 copias y parálogo 83%), #13 Act5C (parálogo 99%), #14 Arf1 (79%)
y #16 eIF4A (73%).

| # | Gen | StringTie | Intrones del CDS confirmados | CDS cubierto | Confianza | Proteasoma | Flag metam. | # en B |
|---|---|---|---|---|---|---|---|---|
| 1 | Rpt1 | discrepancia | 7/9 | 0.97 | alta | sí | | 2 |
| 2 | Rpt5 | confirmado | 7/7 | 1.00 | media | sí | | 4 |
| 4 | Rpt2 | discrepancia | 7/10 | 0.97 | alta | sí | | 8 |
| 5 | EloC | confirmado | 2/2 | 1.00 | alta | | ✔ | 7 |
| 6 | Rpt6 | discrepancia | 6/7 | 0.99 | media | sí | | 9 |
| 7 | mts | confirmado | 5/5 | 1.00 | alta | | ✔ | 10 |
| 8 | Rpn11 | discrepancia | 12/14 | **0.53** | alta | sí | | 14 |
| 9 | Rpn6 | confirmado | 0/0 | 1.00 | media | sí | | 24 |
| 10 | Hsp83 | confirmado | 2/2 | 1.00 | alta | | | 19 |
| 11 | snRNP-U1-70K | confirmado | 8/8 | 1.00 | alta | | | 16 |
| 12 | Sar1 | discrepancia | 4/6 | 1.00 | alta | | ✔ | 5 |
| 15 | Rpt4 | confirmado | 8/8 | 1.00 | alta | sí | | 30 |
| 17 | Hel25E | discrepancia | 7/11 | 0.96 | alta | | | 31 |
| 18 | CG1703 | discrepancia | 10/13 | 0.98 | alta | | | 26 |
| 19 | Vha26 | discrepancia | 3/5 | 1.00 | media | | | 41 |
| 20 | hay | discrepancia | 10/16 | 0.98 | alta | | | 39 |
| 21 | Rpt3 | confirmado | 8/8 | 1.00 | alta | sí | | 17 |
| 22 | TER94 | discrepancia | 12/14 | 0.93 | alta | parcial | ✔ | 42 |
| 23 | gskt/sgg | discrepancia | 7/8 | **0.86** | alta | parcial | ✔ | 51 |
| 24 | VhaSFD | confirmado | 9/9 | 1.00 | alta | | | 54 |

- **9 confirmados y 11 con discrepancia.** En las discrepancias, Helixer predice intrones que ninguna lectura salta.
- `tablas/candidatos/candidatos_20_verificados.fna`: lleva la secuencia de StringTie (respaldada por lecturas) cuando
  cubre ≥ 90% del CDS, en 18/20. Los que discrepan van marcados `REVISAR`. **Rpn11** y **gskt/sgg** van con la secuencia de
  Helixer marcada `NO verificado`.
- 19 de los 20 son `letal_ambos_regimenes`; la excepción es TER94.

## 9. Competencia vectorial (tarea I; lista aparte)

- **Fuente:** Tablas 3 y 4 de Xu et al. 2021, experimento a 25 °C, PRJNA579843 sin lecturas públicas.
  - Son 83 transcriptos **con anotación**, de los 72 + 67 significativos del paper.
  - 27 mapean a Helixer y **24 están en el ranking**: `tablas/candidatos/candidatos_competencia_vectorial.tsv`.
- **Sesgo declarado:** solo los transcriptos anotados aparecen en las tablas. Eso sesga hacia genes conservados y excluye a
  los específicos de linaje, que es justo donde estarían las proteínas salivales de hemípteros y los mejores candidatos por
  especificidad. **Esta lista cubre la fracción conservada de la respuesta a MRFV, no la respuesta completa.** La Tabla S7
  da fold changes sin p ajustados: no se reconstruyó significancia ni se inventaron umbrales.
- 38 de los 83 no alinean al genoma (contra el 96% de alineamiento de todo GITV01). No se encontró la causa.
- **Con letalidad ≥ p75 o CV_estadios en el cuarto más estable (14):**
  - TER94 (#22, up 7 d)
  - VhaPPA1-1 (#32, up 4 h)
  - Slu7 (#260)
  - RpL6 (#345)
  - UQCR-14 (#505)
  - Sod1 (#1182, up 4 h)
  - Karyβ3
  - msl-3
  - Acsl
  - CG1492
  - CG10211
  - Ugt36
  - Tsp29Fb
  - un gen sin ortólogo en *Drosophila* (#1487)
- **Los log2FC de PRJNA631706 (30 °C, réplicas técnicas colapsadas) no reproducen la dirección del paper (25 °C):**
  TER94 da −0.02 a 7 d y VhaPPA1-1 −0.18 a 4 h. Son experimentos distintos, y los log2FC se informan solo como
  descriptivos.

## 10. Limitaciones

1. La letalidad viene de un holometábolo (*Tribolium*). dali22 cuenta muertes en la pupa y el adulto.
   - El flag de metamorfosis marca 15 del top 100, todos genes pleiotrópicos.
   - La firma GO no tiene términos metamórficos.
   - Sin dali22, el top 100 conserva 92/100.
2. Circularidad de la firma GO: el entrenamiento incluye a los controles.
   - Solo 6/77 términos son exclusivamente proteasomales.
   - Sin GO (ranking B), el proteasoma sigue con 15 genes en el top 100.
   - 33/77 términos son genéricos.
3. La expresión constitutiva depende de 3 corridas sin réplicas (una por estadio). El CV entre 3 valores es una estimación
   ruidosa.
4. PRJNA631706 solo tiene réplicas técnicas: no admite tests de expresión diferencial. Los log2FC son descriptivos.
5. Helixer tiene límites de exones discordantes en 11 de los 20 finales. Para diseñar el dsRNA hay que usar la secuencia
   de StringTie y revisar los marcados. Rpn11 y gskt/sgg no están verificados.
6. La lista de MRFV está incompleta (solo transcriptos anotados) y viene de un experimento sin lecturas públicas.
7. Sin RNA-seq por tejido no se puede evaluar la expresión en el intestino medio, que es la vía oral del dsRNA.
8. El proteoma de Helixer no es completo: hay transcripción fuera de los genes anotados.
9. Errores propios corregidos en esta versión:
   - Cobertura sin `-split` (soporte inflado de 76.7% a 92.1%).
   - Imputación de dali22 por la media.
   - Expresión con Salmon solo-CDS.
   - Loci calculados solo para el top 100 anterior.
   - En la versión previa: identidades de minimap2 sin `-c` e intrones de UTR contados en la tarea D.
