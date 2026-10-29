# Resultados lote 1 — camino rápido a OrthoFinder

Fecha: 2026-10-27, 23:30. Cómputo: SDumont II, nodo reservado `sdumont2nd2015` (`lncc-h100`, 4× H100, 96 CPUs).
Scripts en `sbatch_sdumont/`; resultados en el cluster, en `/scratch/PROYECTO/USUARIO/dmaidis/` (sin backup).
La copia local de lo esencial de la pasada 1 está en `orthofinder_pasada1/`.

## Cambio respecto del brief 3: el proteoma de *D. maidis* sale de Helixer, no de TransDecoder

El CHEQUEO 1 **falló** con TransDecoder: BUSCO insecta **C 60.0% / S 41.4% / D 18.6% / F 6.1% / M 33.9%** (23,848 proteínas,
corrido con `--single_best_only`). El problema es el transcriptoma: GITV01 tiene 56,117 contigs (N50 2.4 kb), pero en ~32 mil
no hay ORF ≥ 100 aa y falta un tercio de los BUSCO. CD-HIT no arregla eso, así que no se siguió por ese camino.

Se decidió (con tu OK) predecir genes directamente sobre el genoma (`GCA_975739785.1`, BUSCO genoma C 99.6% / D 1.8%) con
**Helixer v0.3.6**, modelo `invertebrate_v0.3_m_0100`. El genoma se partió en 4 bloques de cromosomas (222–266 Mb), uno por
GPU; tardó 1 h 12 min en total. Proteínas extraídas con `gffread`, una por gen, con prefijo `Dmai|`.

## 1. CHEQUEO 1 — BUSCO del proteoma de *D. maidis* (Helixer)

| Linaje | C | S | D | F | M | n |
|---|---|---|---|---|---|---|
| insecta_odb10 | **93.0%** | 91.4% | **1.6%** | 4.1% | 2.9% | 1367 |
| hemiptera_odb10 | 94.1% | 91.9% | 2.3% | 2.4% | 3.4% | 2510 |

**Veredicto: pasa** (C > 85%, D < 10%). Además, al venir del genoma, cada proteína es un locus: no hay isoformas ni alelos
mezclados, así que el conteo de parálogos es más confiable que con TransDecoder y CD-HIT no hace falta.
22,199 genes predichos.

## 2. Los 10 proteomas (una isoforma por gen)

| Especie | Fuente | Proteínas |
|---|---|---|
| *Dalbulus maidis* | Helixer sobre `GCA_975739785.1` | 22,199 |
| *Macrosteles quadrilineatus* | `GCF_028750875.1` | 21,979 |
| *Homalodisca vitripennis* | `GCF_021130785.1` | 19,904 |
| *Nilaparvata lugens* | `GCF_014356525.2` (2020) | 18,139 |
| *Acyrthosiphon pisum* | `GCF_005508785.2` | 17,680 |
| *Daphnia pulex* | `GCF_021134715.1` | 15,295 |
| *Drosophila melanogaster* | `GCF_000001215.4` (RefSeq) | 13,953 |
| *Bemisia tabaci* | `GCF_918797505.1` | 13,438 |
| *Tribolium castaneum* | `GCF_000002335.3` (Tcas5.2) | 12,873 |
| *Apis mellifera* | `GCF_003254395.2` | 9,934 |

Filtrado con `una_isoforma_por_gen.py` (la proteína más larga por `GeneID`, ≥ 30 aa).
*Macrosteles* va con el proteoma filtrado, sin colapsar: tiene haplotipos (ver `resultados_pipeline.md`), así que **no se usa
para contar parálogos**.

## 3. CHEQUEO 2 — `Statistics_PerSpecies.tsv` (OrthoFinder 2.5.5, pasada 1, DIAMOND, sin MSA)

| Especie | Genes | En ortogrupos | % en OGs | Sin asignar |
|---|---|---|---|---|
| *Macrosteles* | 21,979 | 21,265 | 96.8 | 714 |
| *Tribolium* | 12,873 | 11,971 | 93.0 | 902 |
| *Homalodisca* | 19,904 | 18,439 | 92.6 | 1,465 |
| *Acyrthosiphon* | 17,680 | 16,345 | 92.4 | 1,335 |
| *Apis* | 9,934 | 9,178 | 92.4 | 756 |
| *Bemisia* | 13,438 | 12,331 | 91.8 | 1,107 |
| *Nilaparvata* | 18,139 | 16,634 | 91.7 | 1,505 |
| ***Dalbulus maidis*** | 22,199 | 19,884 | **89.6** | **2,315** |
| *Daphnia* | 15,295 | 13,856 | 90.6 | 1,439 |
| *Drosophila* | 13,953 | 11,781 | 84.4 | 2,172 |

**Veredicto: pasa.** *D. maidis* (89.6%) queda dentro del rango de las demás (84.4–96.8%), a ~2 puntos de la mediana.
Tiene la mayor cantidad de genes sin asignar (2,315, 10.4%), lo esperable en una predicción ab initio: algunos modelos
espurios o fragmentados que no encuentran homólogos. No es anómalo.

## 4. Ortogrupos y ortología con *Tribolium*

- Ortogrupos totales: **16,941** (4,317 específicos de especie). 91.7% de todos los genes quedaron asignados.
- Ortogrupos con las 10 especies: **4,630**. De copia única en las 10: 1,071.
- *D. maidis* aparece en 11,465 ortogrupos (67.7%).
- *D. maidis* ↔ *Tribolium* (`Orthologues/Orthologues_Dalbulus_maidis/Dalbulus_maidis__v__Tribolium_castaneum.tsv`):
  8,029 relaciones ortólogas, de las cuales **6,528 son 1:1**.

## 5. Pasada MSA

Lanzada el 2026-10-27 a las ~23:30 (job `743584`, `orthofinder -b .../WorkingDirectory -M msa -t 96 -a 8`), **en curso**.
Según el brief: si a media mañana no terminó, se mata y se trabaja con la pasada 1.

## 6. BRAKER3 en segundo plano

**Sin lanzar** (decisión tuya: "no por ahora"). Todo está listo para lanzarlo: el genoma, las 19 corridas de RNA-seq en `.sra`
(PRJNA631706 + PRJNA318849), OrthoDB12 Arthropoda (16 M de proteínas) y los contenedores BRAKER 3.0.8 y RepeatModeler 2.0.7.
Scripts: `lote2_t4_t5a.sbatch` y `lote3_t5b_braker.sbatch`. Si se lanza, bajar RepeatModeler a 24–32 hilos.
Con Helixer en C 93%, BRAKER pasa a ser una mejora opcional, no un reemplazo necesario.

## Insumos de la tabla maestra (Tarea 4, ya descargados en `tabla_maestra/`)

- `iBeetle_lethality.csv` (15,213 iB) e `iB_TC.csv` (iB → TC; 509 iB mapean a más de un TC). Columnas: `dali11`/`dali22` =
  % de larvas muertas 11/22 días después de la inyección en larvas; `dapi11` = 11 días después de la inyección en pupas.
  Con **`dali11 ≥ 80%`** (umbral propuesto) quedan 1,499 iB (1,496 con TC). Con `dali11 ≥ 80` y `dapi11 ≥ 80`, 398.
- ⚠ **IDs de *Tribolium*:** el GFF de `GCF_000002335.3` trae el ID TC (BEETLEBASE) solo en **632 genes**, así que el supuesto
  del brief ("IDs compatibles con iBeetle") no se cumple directamente. Se resolvió con NCBI Gene (`locus_tag` =
  `TcasGA2_TCxxxxxx`): `Tcas_GeneID_TC.tsv` mapea **10,729 de 12,480** genes codificantes. Los headers del proteoma son
  `Tcas|NP_… gene=GeneID:…`.
- GO de *Drosophila*: `Dmel_FBgn_GO.tsv` (GAF de FlyBase del GO Consortium, sin anotaciones NOT) y `Dmel_GeneID_FBgn.tsv`
  (desde el GFF de RefSeq). 12,736 de los 13,953 genes del proteoma tienen GO. El FTP de FlyBase no respondía.

## Limitaciones

- El proteoma de *D. maidis* es ab initio (Helixer, sin evidencia de RNA-seq). Modelos de genes con bordes de exones
  aproximados; F 4.1%. Para el dsRNA conviene validar cada candidato contra el transcriptoma GITV01 (blastn/minimap2) antes de
  diseñar.
- Parálogos: con Helixer, cada gen es un locus. La Tarea 5 del brief (minimap2 de candidatos sobre el genoma) sigue siendo
  barata y útil para confirmar que no hay copias muy similares, pero ya no es imprescindible.
- *Drosophila* desde RefSeq y no desde FlyBase (el brief 3 lo acepta).
- BUSCO de los 9 proteomas de referencia no se corrió (el brief 3 lo deja para después).
