# Resumen — sesgo de metamorfosis y cierre del ranking

Fecha: 2026-10-28. Script: `scripts/metamorfosis.py` (log completo en `resultados/metamorfosis/log_metamorfosis.txt`).
Regeneración: `resultados/logs/regeneracion_tarea4.log`. Ontología: `go-basic.obo` release **2026-07-26**, la misma con la
que se construyó la firma GO.

## Tarea 0 — Conjunto de entrenamiento de la firma GO

- Archivo: `set_entrenamiento_GO.txt`, junto al brief y en `resultados/`. No había plantilla adjunta en el directorio; se
  armó con los campos pedidos.
- Criterio proteasomal: el ortólogo de *Drosophila* está anotado a GO:0000502 *proteasome complex* → sí. Solo a
  GO:0010498 *proteasomal protein catabolic process* → parcial.
- **Proteasomales: 3 de 17** (Prosβ5/PSMB5, Rpn7, Rpt3). Parcial: 1 (Hsc70-3).
- **Términos de la firma que provienen solo de genes proteasomales: 6 de 77 (7.8%).**

## Tarea 1 — Flag de metamorfosis

- 7 semillas: todas resuelven al nombre esperado.
- Conjunto expandido por DAG (is_a, part_of) + nombre:
  - **estricto: 146 términos**
  - amplio: 206 términos
- Exclusiones aplicadas: 64.
  - 4 por muda, ecdisis, quitina, cutícula o embrión.
  - 60 por aplicar igual a un hemimetábolo (ecdisona, hormona juvenil, instar o larva sin marcador metamórfico). Estos
    quedan solo en el conjunto amplio.
- Genes marcados (conjunto estricto):

| | Marcados | % |
|---|---|---|
| Universo RefSeq | 625 / 7,881 | 7.9% |
| Universo OGS3 | 646 / 8,651 | 7.5% |
| Top 100 (RefSeq) | **15** | |
| 20 candidatos finales | **5** | |

- Con el conjunto amplio: RefSeq 658 (8.3%), OGS3 698 (8.1%); el top 100 sigue en 15 y los finales en 5.
- **Top 100 marcados (rank):** #5 EloC, #7 mts, #12 Sar1, #14 Arf1, #22 TER94, #23 gskt/sgg, #26 Vha68, #29 Rho1, #37 tsr,
  #50 wds, #52 EloB, #62 rept, #72 hyx, #90 sqh, #97 Cdk8.
- **20 finales marcados:** EloC, mts, Sar1, TER94, gskt/sgg.
- **Criterio de lectura:** ≤ 5 → una línea en el póster; ≥ 20 → discutir reponderar `dali22`; entre 5 y 20 → inspección
  caso por caso. **Resultado: 15 → entre 5 y 20, inspección caso por caso.**

## Tarea 2 — Auditoría de los 77 términos de la firma

- **Intersección de la firma con el conjunto metamórfico (amplio): 0.**
- Términos metamórficos: **0 / 77**.
- Genéricos (profundidad ≤ 3 o anotados en > 25% del universo): **33 / 77**.
  - 33 por profundidad.
  - 6 por fracción del universo.
- Profundidad: mínimo 1, mediana 4, máximo 10.
- Solo de genes proteasomales: **6 / 77 (7.8%)**.
- Archivo: `metamorfosis/auditoria_firma_GO.tsv`.

## Tarea 3 — `letal_ambos_regimenes`

"Alto" = ≥ percentil q de los valores observados. Para `dali` se usa `dali22` observado, o `dali11` observado cuando
`dali22` es imputado. Observados: dali 7,710 genes; dapi11 7,481.

| q | Umbral dali | Umbral dapi11 | `letal_ambos` | Top 100 con `letal_ambos` | Spearman(score, `letal_ambos`) | `False` ∩ flag_metamorfosis |
|---|---|---|---|---|---|---|
| **0.75 (principal)** | ≥ 70 | ≥ 40 | 956 (12.1%) | **93** | **0.448** | 521 de 625 marcados |
| 0.50 | ≥ 40 | ≥ 30 | 2,173 (27.6%) | 96 | 0.459 | 420 de 625 |
| 0.90 | ≥ 100 | ≥ 70 | 421 (5.3%) | 76 | 0.338 | 576 de 625 |

- Fracción marcada por metamorfosis (q = 0.75): **7.5%** entre los `False` y **10.9%** entre los `True`.
- Ranking con filtro (q = 0.75): **sobreviven 93 del top 100 original**. Spearman entre el rank filtrado y el original
  sobre los sobrevivientes = 1 por construcción, porque el filtro no reordena. El top 100 filtrado llega hasta el rank
  original 107.
- Archivos: `metamorfosis/ranking_filtro_ambos_regimenes.tsv` y `metamorfosis/comparacion_filtro.txt`. Es diagnóstico:
  no reemplaza al ranking oficial.

## Tarea 4 — Regeneración completa

- Orden: expresion_matrices → score_candidatos → score_candidatos ogs3 → comparar_ogs3 → rnaseq_integracion →
  verificar_candidatos. Los seis pasos terminaron con exit 0.
- Cobertura con `-split`: confirmado (`rnaseq_insumos/cobertura_exones_split.txt`). Ningún paso recalcula cobertura sin
  `-split`.
- Job 744221 (minimap2 con todos los mRNA): completo, 22,185 de 22,199 mRNA con algún alineamiento. En los 104 mRNA del
  top 100 previo coincide con el PAF anterior.

| Cifra | Corrida previa | Regenerado | Delta |
|---|---|---|---|
| Universo (RefSeq / OGS3) | 7,881 / 8,651 | 7,881 / 8,651 | 0 |
| Genes con soporte de lecturas | 17,023 / 22,199 | 17,023 / 22,199 | 0 |
| iB letales en el ranking (RefSeq / OGS3) | 961 / 1,015 | 961 / 1,015 | 0 |
| `paralogo_resuelto` (ambas / solo_una / ninguna) | 857 / 35 / 4 | 857 / 35 / 4 | 0 |
| Confianza top 100 (alta / media / baja) | 70 / 30 / 0 | 70 / 30 / 0 | 0 |
| `loci_genoma` top 100 | 85 con 1 locus + 15 sin dato | 100 con 1 locus | +15 con dato |
| Top 100 compartido | — | 100 / 100 | Spearman 0.9999 |

- Contra el ranking del inicio de la sesión (antes del cambio de `dali22` y de expresión): top 100 compartido 85/100,
  Spearman 0.988; confianza del top 100 68/32/0 → 70/30/0.
- Los 20 finales cambiaron contra la versión anterior del brief RNA-seq. Siguen 12; salen Prosα2, Prosα5, Prosα6, Rab1,
  Polr2E, CkIIα, RpLP0 y Vha14-1; entran EloC, mts, Rpn6, Hsp83, Sar1, Vha26, TER94 y gskt/sgg. De los 20: 9 confirmados
  y 11 con discrepancia con StringTie.

## Tarea 5 — Controles positivos

| Control | Rank previo | Rank regenerado | Desplazamiento |
|---|---|---|---|
| PSMB5 | #25 | #25 | 0 |
| Rpt3 | #21 | #21 | 0 |
| V-ATPasa D | #125 | #125 | 0 |
| V-ATPasa B | #216 | #216 | 0 |

PSMB5 y Rpt3 están en el top 50.

## Tarea 6 — Ranking B, sin firma GO

- Pesos redistribuidos manteniendo las proporciones relativas (cada peso / 0.70): letalidad 50%, expresión 28.6%,
  conservación 21.4%, como pide el brief.
- **Top 100 de A que siguen en el top 100 de B: 63/100.** Spearman A vs. B: **0.898**.
- **Genes proteasomales (GO:0000502) en el top 100 de B: 15** (en A: 22). Parciales: B 6, A 7.
- Controles en B:

| Control | A | B |
|---|---|---|
| PSMB5 | #25 | #34 |
| Rpt3 | #21 | #17 |
| V-ATPasa D | #125 | #65 |
| V-ATPasa B | #216 | #184 |

- **Entran al top 100 solo en B (37):**
  - Sec61α (B #1)
  - Sec23
  - Prp31
  - Gs1/Gs2
  - cpb
  - l(1)10Bb
  - Bx42
  - Snx3
  - ATPsynγ
  - AsnRS
  - Pp1 (Pp1α-96A)
  - Acly
  - TfIIB
  - RpLP1
  - γCOP
  - Arp2
  - Sec24CD
  - eIF3e
  - CG10333
  - Cchl
  - Sams
  - c12.1
  - Vha36-1
  - δCOP
  - eIF3c
  - eRF1
  - Chmp1
  - Sf3b2
  - Eno (×2)
  - β'COP
  - eIF3a
  - Gβ13F
  - ATPsynF
  - Nop56
  - eIF3d
  - Sf3b3
  - Lista completa en `metamorfosis/ranking_B_sin_GO.tsv`.
- **Top 100 de B con flag_metamorfosis: 18** (A: 15).
- **Lectura (criterio del brief):** los proteasomales siguen en el top 100 de B sin usar GO (15).
