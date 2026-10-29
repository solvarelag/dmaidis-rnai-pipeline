# Score de candidatos — letalidad continua, copias como filtro

Fecha: 2026-10-27/28. Script: `scripts/score_candidatos.py` (en el docstring está el método de normalización).
Salidas: `ranking/ranking_completo.tsv` (score desglosado: aporte, valor normalizado y valor crudo de cada criterio, más
flags e IDs), `ranking/firma_go.tsv` y `ranking/reporte_score.txt`.

## 1. Universo

Genes de *Tribolium* con dato de letalidad en iBeetle (sin corte) **y** con ortólogo en *D. maidis*: **8,272**.

| Copias en *D. maidis* | Genes | Acción |
|---|---|---|
| 1 | 6,985 (84.4%) | pasa |
| 2 | 896 (10.8%) | pasa con `riesgo_paralogo` |
| ≥ 3 | 391 (4.7%) | descartado |

**Quedan 7,881 en el ranking.**

## 2. Letalidad continua

`0.5·dali22 + 0.3·dali11 + 0.2·dapi11`. Si un gen tiene varios iB, se promedia por columna. Los faltantes se imputan con la
media del universo.

- **iBeetle codifica "no medido" como `-1`**, además de las celdas vacías: hay 317 en `dali22` y 913 en `dapi11`. Se trataron
  como faltantes.
- Imputaciones: `dali22` **4,658**, `dapi11` 400, `dali11` 171. **4,707 de los 7,881 genes (60%)** tienen al menos una.
- ⚠ **El faltante de `dali22` no es al azar:** la fase II del screen casi nunca lo midió (falta en 9,846 de 9,868 iB). Además,
  la mortalidad es acumulada: en las 4,425 iB con ambas mediciones, `dali22 ≥ dali11` en el **99.9%**. La media (46.95)
  produce valores imposibles: un gen con `dali11 = 100` queda con `dali22 = 47` y baja en letalidad (p. ej. *Rpt2*, *Rpt5* y
  *Rpn6* tienen `L = 0.93` en vez de 0.99). **Propuesta:** imputar `dali22 = max(dali11, media)`. Se mantiene el espíritu de
  no castigar al gen y se respeta la monotonía. Es un cambio de una línea; no lo apliqué porque la instrucción fue explícita.

## 3. Firma GO

- GO de *Drosophila* transferidos vía el ortólogo *Tribolium*–*Drosophila* de OrthoFinder y propagados a ancestros
  (go-basic.obo, 2026-07-26).
- **Fondo:** genes de *Tribolium* con ortólogo 1:1 en *D. maidis* y con GO: **5,739**.
- **Entrenamiento:** 17 genes, todos mapeados; 16 están en el fondo 1:1.
  - Los 11 de Ulrich et al. 2015 (Fig. 1 y Tabla S1): *cact*, *Srp54k*, *Rop*, *αSnap*, *shi*, *Pp1α-96A*, *Inr-a*,
    *Hsc70-3*, *Rpn7*, *gw*, *Rpt3*.
  - *dre4*, *ncm* y *RpII140* (= *Polr2B*) de Knorr et al. 2018.
  - *Prosβ5* (PSMB5).
  - V-ATPasa **B (*Vha55*) y D (*Vha36-1*)**: las dos subunidades usadas en *P. maidis* (Yao et al. 2013).
- Hipergeométrica, p < 0.01 y ≥ 2 genes de entrenamiento: **77 términos**. El puntaje de un gen es la suma de −log10(p) de los
  términos de la firma que tiene.
- Para que los propios genes de entrenamiento no se puntúen con una firma que los incluye, cada uno se puntúa con la firma
  recalculada sin él (*leave-one-out*).

## 4. Expresión y conservación

- **Expresión:** salmon 1.10.3 sobre las CDS de Helixer, **19 corridas** (PRJNA631706 + PRJNA318849), con 35–76% de lecturas
  mapeadas. El subscore promedia el percentil de TPM mínimo y el percentil inverso del CV. 12 genes sin dato.
- **Conservación:** % de identidad del mejor HSP de DIAMOND (`--more-sensitive`) entre el gen de *Tribolium* y su ortólogo de
  *D. maidis*. Cuando hay 2 copias se usa la de mayor identidad, y la expresión se toma de esa misma copia. **259 genes sin
  hit** quedan con 0 en esa columna.

## 5. Normalización

Rango percentil (0–1) por columna, explicado en el docstring del script. Los genes sin dato en una columna reciben 0 en esa
columna.

## 6. Controles positivos (de 7,881)

| Control | Posición | Percentil | L | GO | E | C |
|---|---|---|---|---|---|---|
| **PSMB5** (*Prosβ5*, blanco de ledprona) | **#24** | top 0.3% | 0.99 | 0.97 | 0.99 | 0.84 |
| **V-ATPasa B** (*Vha55*) | **#149** | top 1.9% | 0.91 | 0.82 | 0.94 | 0.98 |
| **V-ATPasa D** (*Vha36-1*) | **#220** | top 2.8% | 0.92 | 0.83 | 0.85 | 0.93 |

Los dos controles quedan arriba. **Pero la validación es parcial:** los dos están en el conjunto de entrenamiento de la firma
GO. El *leave-one-out* evita que se puntúen con su propio aporte, pero otros miembros de su complejo (*Rpt3*, *Rpn7*; la otra
subunidad de V-ATPasa) siguen en el entrenamiento y empujan los términos de proteasoma y V-ATPasa. Sin la columna GO,
PSMB5 igual sería top: L 0.99, E 0.99, C 0.84.

Resto del entrenamiento: *Rpt3* #5, *Hsc70-3* #45, *RpII140* #105, *ncm* #117, *Srp54k* #146, *Pp1α-96A* #215, *Rpn7* #241,
*dre4* #271, *αSnap* #471, *Rop* #528 (`riesgo_paralogo`), *shi* #955, *gw* #1404, *cact* #1599, *Inr-a* #2193.
Los peor rankeados lo son por conservación baja (*cact* C = 0.08, *Inr-a* C = 0.05, *gw* C = 0.39): son genes divergentes
entre escarabajo y cicadélido.

## 7. Recuperación de IDs (DIAMOND recíproco OGS3 ↔ RefSeq)

- RBH a nivel de gen (`--id 90 --query-cover 70 --subject-cover 70`, isoformas colapsadas al máximo bitscore): **10,511 pares
  TC–GeneID** (`tabla_maestra/Tcas_TC_GeneID_rbh.tsv`).
- De los 402 letales (`dali11 ≥ 80`) perdidos, se recuperaron **48**. **Quedan afuera 354 iB (363 TC):** todos están en OGS3,
  pero solo 2 tienen algún hit directo que pase los umbrales. Son modelos de genes que cambiaron entre OGS3 y RefSeq Tcas5.2.
- La ortología que publica iBeetle-Base (`1.OrthoFinder_orthology.tsv`) es solo *Tribolium* ↔ *Drosophila* (7,505 pares), así
  que **no** resuelve el mapeo hacia *D. maidis*. La solución de raíz es **usar el proteoma OGS3 de *Tribolium* en OrthoFinder**
  en lugar de RefSeq: los IDs ya serían TC y no habría pérdida de mapeo. La pasada 1 tarda ~6 min.
- Conteo de copias de los letales `dali11 ≥ 80` rehecho con RBH: 1,142 con GeneID → 980 con ortólogo → **1 copia 870
  (88.8%) · 2 copias 89 (9.1%) · ≥3 copias 21 (2.1%)**. Los 48 recuperados no cambian la distribución.

## 8. ¿Domina algún criterio?

- Aporte medio en el top 50: letalidad 0.336/0.35, GO 0.293/0.30, expresión 0.192/0.20, conservación 0.142/0.15. Los cuatro
  saturan juntos y ninguno domina solo.
- **Sesgo temático:** 21 del top 100 son subunidades de proteasoma (*Pros*/*Rpn*/*Rpt*), y el top 30 está dominado por
  proteasoma, V-ATPasa, ribosoma y splicing. Es esperable porque el entrenamiento pesa hacia esas funciones (Ulrich 2015 ya lo
  reportó), pero la firma GO tiende a reproducir el entrenamiento. Si se buscan blancos *nuevos*, convendría mirar también un
  ranking sin GO, o limitar los candidatos por complejo.
- Hay 105 genes empatados en la letalidad máxima. En el tope, la letalidad casi no discrimina y el orden lo deciden GO,
  expresión y conservación.
- Ejemplo de lo que aporta `dali22`: *Vha14-1* (#22) tiene `dali11 = 0` pero `dali22 = 100`. Mata lento, y con el corte
  anterior habría quedado afuera.

## Pendientes

- Decidir la imputación de `dali22` (punto 2).
- Opcional: rehacer OrthoFinder con OGS3 de *Tribolium* (punto 7).
- La pasada MSA se relanzó después del job de expresión (job `743624`, corriendo).
