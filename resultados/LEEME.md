# Cómo leer los resultados

Carpeta generada el 2026-10-28 con `scripts/armar_resultados.sh`, después de la regeneración completa. Las tablas son TSV
(separadas por tab, con encabezado) y se abren con Excel/LibreOffice, `pandas.read_csv(sep="\t")` o R `read.delim`.

## Por dónde empezar

1. **`RESULTADOS.md`**: todos los resultados con números, el top 30, los 20 finales y las limitaciones. Reemplaza a
   los informes de `informes_previos/`.
2. **`resumen_metamorfosis_y_cierre.md`**: los números del último brief, sin interpretación, para el póster.
3. **`METODOLOGIA.md`**: cómo se calculó cada cosa: parámetros, umbrales y criterios.
4. **`tablas/ranking/ranking_completo.tsv`**: la tabla maestra, una fila por gen candidato. Todo lo demás se puede
   reconstruir a partir de ella.

## Estructura

```
resultados/
├── LEEME.md, RESULTADOS.md, METODOLOGIA.md, resumen_metamorfosis_y_cierre.md
├── firma_GO_entrenamiento.txt     los 77 términos GO de la firma, de qué genes de entrenamiento sale cada uno y de dónde salió cada gen
├── set_entrenamiento_GO.txt       ficha de los 17 genes de entrenamiento (tarea 0)
├── tablas/
│   ├── ranking/                   ranking oficial (A) y sensibilidad sin dali22
│   ├── ranking_ogs3/              análisis separado con Tribolium OGS3
│   ├── candidatos/                20 finales verificados (+ FASTA) y lista de competencia vectorial
│   ├── expresion/                 matrices de TPM y metadatos por corrida, separadas por BioProject
│   ├── gitv01_mrfv/               puente Helixer↔GITV01 y transcriptos de MRFV del paper
│   └── firma_go/                  firma GO cruda y set de entrenamiento mapeado
├── metamorfosis/                  tareas 1–3 y 6 del brief de metamorfosis
├── logs/                          salida de cada script en la última regeneración
├── informes_previos/              informes parciales anteriores (históricos; algunas cifras quedaron viejas)
(los scripts y sbatch están en `scripts/`, `sbatch/` y `preparacion/` en la raíz del repositorio)
```

## Archivo por archivo

### `tablas/ranking/ranking_completo.tsv` — tabla maestra (7,881 genes, ordenada por `rank`)

La clave es `tcas_geneid` (GeneID de NCBI del gen de *Tribolium*). El gen blanco en *D. maidis* es
`dmai_representativo`.

| Columna | Qué es |
|---|---|
| `rank`, `score` | Posición y score del ranking oficial A (0–1). |
| `*_aporte` | Contribución de cada criterio al score (peso × normalizado). Suman `score`. |
| `*_norm` | Valor normalizado por rango percentil (0–1) de letalidad, go, expresion y conservacion. |
| `letalidad_raw` | 0.5·dali22 + 0.3·dali11 + 0.2·dapi11, en % de mortalidad. |
| `dali11`, `dali22`, `dapi11` | Mortalidad de iBeetle (inyección larval a 11/22 días; pupal a 11 días). Ya incluye las imputaciones. |
| `dali11_imputado`, `dali22_imputado`, `dapi11_imputado`, `n_imputadas` | True si el valor fue imputado. dali22 se imputa como max(dali11, media). |
| `go_raw` | Σ −log10 p de los términos de la firma que porta el gen (leave-one-out para el entrenamiento). |
| `conservacion_raw` | % de identidad DIAMOND entre Tcas y el Dmai representativo. |
| `TPM_minimo`, `CV_estadios`, `sesgo_sexo` | Expresión constitutiva (PRJNA318849): mínimo entre ninfa, hembra y macho; CV entre los tres; log2((H+1)/(M+1)). |
| `tpm_ninfa`, `tpm_hembra`, `tpm_macho` | TPM por estadio. |
| `rank_sin_dali22` | Posición si la letalidad no usara dali22 (sensibilidad). |
| `copias_dmai`, `riesgo_paralogo` | Número de ortólogos en Dmai (1 o 2) y flag cuando son 2. |
| `entrenamiento` | Nombre si el gen pertenece al set de entrenamiento de la firma GO. |
| `tc`, `iB`, `dmel` | IDs TC de OGS3, IDs iBeetle y símbolos de los ortólogos de *Drosophila*. |
| `dmai_representativo`, `dmai_todos` | Gen(es) Helixer de *D. maidis*. El representativo es el de mayor identidad. |
| `rnaseq_cobertura`, `rnaseq_lecturas`, `rnaseq_soporte` | Fracción exónica cubierta (`-split`), fragmentos exónicos y soporte (≥ 0.5 y ≥ 10). |
| `paralogo_resuelto` | Solo con 2 copias: ambas_con_lecturas / solo_una / ninguna. |
| `estabilidad_estadios`, `estabilidad_condiciones` | CV entre estadios (318849) y entre las 4 condiciones de MRFV (631706). Menor = más estable. No entran al score. |
| `log2FC_MRFV_4h`, `log2FC_MRFV_7d` | log2((TPM_MRFV+1)/(TPM_naive+1)) con las réplicas técnicas colapsadas. Descriptivos: **no** indican significancia. |
| `id_GITV01`, `soporte_GITV01` | Mejor hit recíproco en GITV01 y si hay respaldo de GITV01 (RBH o locus). |
| `confianza_gen` | alta (lecturas + GITV01) / media (solo lecturas) / media_revisar (solo GITV01) / baja. |
| `respuesta_MRFV` | up/down × 4h/7d según las Tablas 3/4 de Xu et al. 2021 (solo transcriptos anotados). |
| `loci_genoma` | Loci del mRNA en el genoma (minimap2, ≥ 90% de identidad). 1 = copia única. |
| `paralogos_prot`, `max_ident_paralogo` | Parálogos proteicos en Dmai (≥ 40% de identidad, ≥ 50% de cobertura) y la identidad máxima. ≥ 70% excluye de los finales. |
| `flag_metamorfosis`, `terminos_metamorficos` | True si el gen tiene GO del conjunto metamórfico estricto, y cuáles. |
| `flag_metamorfosis_amplio` | Lo mismo con el conjunto amplio (incluye ecdisona, hormona juvenil e instar). |
| `proteasomal` | si (GO:0000502), parcial (solo GO:0010498) o no. |
| `letal_ambos_regimenes` | Alto en inyección larval y pupal (p75 de los valores observados). Diagnóstico. |
| `score_B`, `rank_B` | Ranking B, sin firma GO (50/28.6/21.4). |

**Cómo usarla:**

- Filtrar `rank <= 100` para el top 100.
- Ordenar por `rank_B` para ver el ranking exploratorio.
- `confianza_gen == "baja"` o `max_ident_paralogo >= 70` son señales de alarma.
- `rnaseq_soporte == False` indica un modelo dudoso.

### `tablas/ranking/sensibilidad_sin_dali22.tsv`

Genes que están en el top 100 con o sin dali22. Compara `rank` con `rank_sin_dali22`. Si un gen tiene rank ≤ 100 y
rank_sin_dali22 > 100, entra al top solo por la mortalidad tardía.

### `tablas/ranking_ogs3/`

- `ranking_completo.tsv`: la misma estructura que la tabla maestra, con clave `tcas_tc` (ID de OGS3), 8,651 genes. Es un
  análisis separado y **no** es el ranking oficial.
- `top100_solo_ogs3.tsv`: genes que entran al top 100 solo con OGS3. `rank_refseq` vacío = el gen no estaba en el ranking
  RefSeq.

### `tablas/candidatos/`

- `candidatos_20_verificacion.tsv`: los 20 finales.
  - `intrones_helixer`: intrones del CDS predichos.
  - `intrones_confirmados`: los que alguna lectura salta exactamente.
  - `intrones_solo_stringtie`: intrones dentro del CDS que Helixer no tiene.
  - `cds_cubierto`: fracción del CDS dentro del transcripto de StringTie.
  - `estado`: confirmado o discrepancia.
- `candidatos_20_verificados.fna`: la secuencia para diseñar el dsRNA. El encabezado dice la fuente:
  - `StringTie ... (coincide con Helixer)`: lista para usar.
  - `REVISAR`: StringTie difiere de Helixer; se usa la de StringTie, respaldada por lecturas.
  - `NO verificado`: secuencia de Helixer sin respaldo completo (Rpn11, gskt/sgg).
- `candidatos_competencia_vectorial.tsv`: los 24 genes del ranking con respuesta a MRFV. Es una lista **aparte**, sesgada
  a genes conservados (ver RESULTADOS §9). `buena_letalidad` = letalidad_norm ≥ 0.75; `expresion_estable` = CV_estadios en
  el cuarto más bajo.

### `tablas/expresion/`

- `metadatos_corridas.tsv`: una fila por corrida (run): BioProject, BioSample, estadio, sexo, tratamiento, tiempo,
  temperatura y fragmentos asignados a genes. Es la clave para leer las matrices.
- `tpm_corridas.tsv`: TPM por gen (`gen_dmai`) × 19 corridas. Las columnas son runs; se cruzan con los metadatos.
  **No mezclar los dos BioProjects en una misma comparación.**
- `estadios_318849.tsv`: TPM por estadio, TPM_minimo, CV_estadios y sesgo_sexo, para todos los genes de Helixer.
- `condiciones_631706.tsv`: TPM de las 4 condiciones (réplicas técnicas sumadas), CV_condiciones y log2FC. Solo
  descriptivos.

### `tablas/gitv01_mrfv/`

- `mapeo_helixer_GITV01.tsv`: gen Helixer → mejor hit recíproco en GITV01 (`pident`) y otros transcriptos que caen en ese
  gen.
- `mrfv_DE.tsv`: los 83 transcriptos de las Tablas 3/4 de Xu et al. 2021 (LFC y padj del paper).
- `mrfv_DE_mapeado.tsv`: lo mismo con `gen` (Helixer) y `via` (locus, rbh o hit_directo; vacío = no mapeó).

### `tablas/firma_go/`

- `firma_go.tsv`: los 77 términos con p, `k_entrenamiento` y `k_fondo`.
- `entrenamiento.tsv`: los 17 genes de entrenamiento con su GeneID, FBgn, y si están en el fondo 1:1 y en el ranking.

### `metamorfosis/`

| Archivo | Qué da |
|---|---|
| `terminos_metamorficos.tsv` | Conjunto de términos (estricto/amplio) y de qué semilla o palabra clave sale cada uno. |
| `exclusiones_aplicadas.tsv` | Términos descartados y el motivo (muda/quitina/cutícula/embrión, o aplica a hemimetábolo). |
| `auditoria_firma_GO.tsv` | Los 77 términos de la firma: profundidad, genes del universo, genes de entrenamiento aportantes, `solo_proteasomales`, `metamorfico`, `generico`. |
| `ranking_filtro_ambos_regimenes.tsv` | Ranking A restringido a `letal_ambos_regimenes` (diagnóstico). |
| `comparacion_filtro.txt` | Números de la tarea 3 con los tres cortes de percentil. |
| `ranking_B_sin_GO.tsv` | Ranking B completo, ordenado por `rank_B`, con el rank A al lado. |
| `log_metamorfosis.txt` | Salida completa del script: todos los números de las tareas 0–3 y 6. |

### `logs/`

- `regeneracion_tarea4.log`: la última regeneración completa, paso por paso, con su exit code.
- `reporte_score.txt` y `reporte_score_ogs3.txt`: universo, imputaciones, firma, sensibilidad y controles.
- `reporte_comparacion_ogs3.txt`: RefSeq contra OGS3.
- `reporte_rnaseq.txt`: tareas B, C, E, F, G, H e I.
- `reporte_candidatos20.txt`: los excluidos antes de llegar a 20 y la verificación de cada final.

## Advertencias para la interpretación

- El score es **relativo** (percentiles): 0.96 no es una probabilidad, sino que el gen está arriba en los cuatro criterios.
- `log2FC_MRFV_*` y `estabilidad_condiciones` vienen de réplicas técnicas: no sirven para afirmar expresión diferencial.
- `flag_metamorfosis` marca genes cuyo ortólogo de *Drosophila* tiene anotaciones de disco imaginal o
  de desarrollo larva-pupa. Muchos son esenciales pleiotrópicos: el flag es una alerta, no una exclusión.
- En los IDs de *D. maidis* (`Dmai|Dalbulus_maidis_<cromosoma>_<n>`), el prefijo indica el cromosoma o contig del genoma
  GCA_975739785.1.
