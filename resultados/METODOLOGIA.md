# Metodología — Priorización de blancos de RNAi en *Dalbulus maidis*

Versión: 2026-10-28. Cada paso indica el script o sbatch que lo implementa (en `scripts/`, `sbatch/` y `preparacion/` del repositorio). Los sbatch
corrieron en SDumont (LNCC): partición `lncc-h100`, una cuenta de proyecto, nodo `sdumont2nd2015` (96 CPU, 4×H100), con
herramientas en contenedores Singularity. Lo demás corrió en local con Python 3.13 (pandas, numpy, scipy).

---

## 1. Genoma y control de calidad (`chequeos_minimos.md`, `t2_busco_genoma`)

- *D. maidis* GCA_975739785.1 (ihDalMaid1.1, 972 Mb, nivel cromosoma). Es el pseudohaplotipo principal; GCA_975739835.1,
  el alternativo, se descartó para no duplicar loci.
- BUSCO 6.1.0 (linaje insecta) sobre el genoma: C 99.6%, D 1.8%.

## 2. Anotación de *D. maidis* (`h1_helixer.sbatch`)

- No hay anotación oficial. Se probaron dos alternativas:
  - **TransDecoder** 5.7.1 sobre el transcriptoma GITV01: BUSCO C 60%, D 18.6%. No pasa el chequeo 1.
  - **Helixer v0.3.6**: modelo `invertebrate`, `--subsequence-length 213840 --overlap-offset 106920
    --overlap-core-length 160380`, el genoma en 4 bloques en paralelo sobre 4 GPU.
- Proteínas con gffread `-y`. Se toma una por gen, se quitan los stops terminales y se descartan las < 30 aa.
- Resultado: 22,199 genes; BUSCO insecta **C 93.0%, D 1.6%** → `proteomas/Dalbulus_maidis.faa`.
- BRAKER3 y RepeatModeler quedaron preparados (`lote2`, `lote3`) y no se lanzaron, por decisión de la usuaria.

## 3. Proteomas y ortología (`h2_orthofinder.sbatch`, `o1_orthofinder_ogs3.sbatch`)

- **10 proteomas**, una isoforma por gen (`una_isoforma_por_gen.py`; la más larga):
  - Hemípteros: *D. maidis* (Helixer), *Macrosteles quadrilineatus*, *Homalodisca vitripennis*, *Nilaparvata lugens*,
    *Acyrthosiphon pisum*, *Bemisia tabaci*.
  - Holometábolos: *Tribolium castaneum* (RefSeq GCF_000002335.3), *Drosophila melanogaster* (RefSeq GCF_000001215.4;
    FlyBase estaba caído), *Apis mellifera*.
  - Grupo externo: *Daphnia pulex*.
- *Macrosteles* tiene duplicación por haplotipos (D 38%, colapsa 23.6% con CD-HIT): se mantiene en la ortología, pero no se
  usa para contar parálogos.
- **OrthoFinder 2.5.5:** `-S diamond -a 8 -t 96`, árbol de especies fijo (`arbol_especies.nwk`), sin MSA en la pasada
  principal. Hubo una pasada MSA aparte, no usada en el score.
- **Análisis separado con OGS3:** el proteoma *Tribolium* OGS3 (`tabla_maestra/OGS3_proteins.fasta.gz`, 18,536 proteínas →
  16,593 genes, isoforma más larga) reemplaza a RefSeq y se corre OrthoFinder con los mismos parámetros y las mismas otras
  9 especies. OGS3 no se agrega como especie número 11, para no tener dos *Tribolium* casi idénticos que distorsionen los
  ortogrupos. Además se corre DIAMOND `--more-sensitive -e 1e-5 -k 25` de Dmai contra OGS3 para la conservación.

## 4. Letalidad (iBeetle)

- Fuente: `tabla_maestra/iBeetle_lethality.csv` (iBeetle-Base).
  - dali11/dali22 = % de mortalidad tras inyección en **larva**, puntuada a 11 y 22 días.
  - dapi11 = inyección en **pupa**, a 11 días.
  - **El valor −1 significa "no medido"** y se trata como faltante, igual que las celdas vacías.
- Mapeo iB → TC (`iB_TC.csv`) → GeneID de RefSeq (`Tcas_GeneID_TC.tsv`, el locus_tag de NCBI Gene). Se completa con:
  - mapeo por secuencia (`Tcas_TC_GeneID_secuencia.tsv`)
  - **RBH OGS3↔RefSeq** (DIAMOND `--id 90 --query-cover 70 --subject-cover 70`, isoformas colapsadas al máximo bitscore;
    `scripts/rbh_ogs3_refseq.py`): 10,511 pares
- Si un gen tiene varios iB, se promedia por columna. En el modo OGS3 el gen es directamente el TC.
- **Universo** = genes con dato de letalidad (sin corte) y ≥ 1 ortólogo en Dmai.
- **Filtro de copias en Dmai:**
  - 1 → pasa
  - 2 → pasa con `riesgo_paralogo`
  - ≥ 3 → se descarta
- **Imputación** (sobre el universo filtrado):
  - dali11 y dapi11 → la media de la columna.
  - **dali22 → max(dali11 ya imputado, media de dali22).** La mortalidad es acumulada: dali22 ≥ dali11 en el 99.9% de las
    4,425 iB con ambas medidas.
  - Las columnas `*_imputado` registran qué se imputó.
- **Letalidad** = 0.5·dali22 + 0.3·dali11 + 0.2·dapi11.
- **Sensibilidad:** 0.6·dali11 + 0.4·dapi11, sin dali22 → `rank_sin_dali22`.
- Conteo de letales (`scripts/conteo_letales.py`): iB con dali11 ≥ 80.

## 5. Score (`scripts/score_candidatos.py`; modo `ogs3` para el análisis separado)

**Score = 0.35·letalidad + 0.30·GO + 0.20·expresión + 0.15·conservación.**

Cada columna cruda se normaliza por **rango percentil**: `(rank − 1)/(n − 1)` sobre el universo filtrado, con rango
promedio en los empates. Se eligió percentil y no min-max porque TPM y −log10 p tienen colas largas; el costo es que el score
expresa posiciones relativas. Un gen sin dato en una columna recibe 0 en esa columna.

### 5.1 Firma GO (30%)

- Anotaciones de *Drosophila* (`Dmel_FBgn_GO.tsv`, FlyBase, todas las evidencias), propagadas a ancestros por is_a y part_of
  (`go-basic.obo` releases/2026-07-26).
- Se transfieren al gen de *Tribolium* a través del ortólogo Tcas–Dmel de OrthoFinder (unión si hay varios FBgn).
- **Entrenamiento (17 genes):**
  - 11 de Ulrich et al. 2015 (Fig. 1 y Tabla S1)
  - dre4, ncm y RpII140 de Knorr et al. 2018
  - Prosβ5/PSMB5 (ledprona; Rodrigues et al. 2021)
  - V-ATPasa B y D (Yao et al. 2013, *P. maidis*)
  - Ficha completa en `set_entrenamiento_GO.txt`.
- **Fondo** = genes Tcas con ortólogo 1:1 en Dmai y con GO (5,739); 16 de los 17 están en el fondo.
- **Test:** hipergeométrico de cola superior. Un término entra a la firma con p < 0.01 y ≥ 2 genes de entrenamiento: quedan
  77 términos.
- **Puntaje** = Σ −log10(p) de los términos de la firma que porta el gen.
- Los genes de entrenamiento se puntúan con una firma **recalculada sin ellos** (leave-one-out).

### 5.2 Expresión constitutiva (20%) — solo PRJNA318849 (`scripts/expresion_matrices.py`)

- Corridas: SRR5010611 (ninfa, sexo mixto), SRR5010612 (hembra adulta), SRR5010613 (macho adulto).
- **Cuantificación:** fragmentos por gen y corrida sobre los BAM de HISAT2 (`r4_conteos_por_corrida.sbatch`):
  - read1, primario, mapeado, MAPQ ≥ 10 (`samtools view -q 10 -f 0x40 -F 0x904`)
  - cuenta una vez por gen si algún bloque del alineamiento cae en un exón (`bedtools coverage -split -counts` contra un
    BED12 con los exones de Helixer)
  - TPM = (fragmentos / kb exónicas), reescalado a una suma de 10⁶ por muestra
- Reemplaza a Salmon sobre CDS de Helixer, que perdía las lecturas de UTR y daba ceros espurios.
- **TPM_minimo** = mínimo entre los 3 estadios. **CV_estadios** = desvío muestral / media.
- **sesgo_sexo** = log2((TPM_hembra + 1)/(TPM_macho + 1)). Se informa y **no** entra al score.
- **Subscore** = 0.5·pct(log10(TPM_minimo + 0.01)) + 0.5·(1 − pct(CV_estadios)); después se renormaliza por percentil.
- Se usa el gen Dmai representativo: el de mayor identidad con el gen Tcas si hay 2 copias.

### 5.3 Conservación (15%)

- % de identidad del mejor HSP de DIAMOND `--more-sensitive -e 1e-5` entre el gen Tcas y su ortólogo Dmai representativo
  (`s4_expr_rbh_ident.sbatch`).

### 5.4 Ranking B (`scripts/metamorfosis.py`, tarea 6)

- Sin GO. Los pesos se redistribuyen manteniendo las proporciones relativas (peso / 0.70): letalidad 50%, expresión 28.6%,
  conservación 21.4%. Se reutilizan las mismas columnas normalizadas.

## 6. RNA-seq (`r1`–`r6`, `scripts/rnaseq_integracion.py`, `scripts/verificar_candidatos.py`)

### 6.1 Datos y separación por BioProject

- Las 19 corridas se bajaron del espejo AWS de SRA con sra-tools.
- Metadatos (BioSample, estadio, sexo, tratamiento, tiempo) de ENA/BioSample → `tablas/expresion/metadatos_corridas.tsv`.
- **PRJNA318849** → expresión constitutiva y estabilidad entre estadios.
- **PRJNA631706** → 4 BioSamples × 4 corridas = **réplicas técnicas**: naive/MRFV × 4 h/7 d, 30 °C.
  - No se aplica ningún test de expresión diferencial: las réplicas técnicas subestiman la varianza.
  - Las 4 corridas de cada BioSample se **colapsan sumando fragmentos** antes de calcular TPM.
  - log2FC_MRFV_4h/_7d = log2((TPM_MRFV + 1)/(TPM_naive + 1)): descriptivos, como variable continua y sin umbral.
- PRJNA579843 (Xu et al. 2021) no tiene lecturas públicas. Se verificó en ENA, en SRA (esearch/elink) y en las BioSamples
  del registro TSA de GITV01.

### 6.2 Alineamiento (tarea A; `r1_rnaseq_genoma.sbatch`)

- HISAT2 (`--dta`) contra el genoma, por corrida; samtools sort/index; merge de las 19.
- StringTie 3.0.3 **de novo** (sin `-G`) sobre el BAM combinado → `stringtie_denovo.gtf` y `stringtie_tx.fa`.

### 6.3 Soporte por lecturas (tarea B; `r5_cobertura_split.sbatch`)

- `bedtools coverage -split` de los exones de Helixer contra el BAM combinado, un proceso por cromosoma →
  `cobertura_exones_split.txt`.
- `rnaseq_cobertura` = fracción de bases exónicas (UTR incluidas) con ≥ 1 lectura.
- `rnaseq_lecturas` = suma de los fragmentos exónicos de 6.1.
- `rnaseq_soporte` = cobertura ≥ 0.5 **y** ≥ 10 fragmentos.
- La versión anterior no usaba `-split` y contaba los intrones (N del CIGAR) como cobertura: 92.1% contra 76.7% con
  `-split`.

### 6.4 Puente con GITV01 (tarea G; `scripts/puente_gitv01.py`, `r2`/`r3`)

- RBH de proteína Helixer ↔ proteínas TransDecoder de GITV01: DIAMOND `--id 80 --query-cover 60 --subject-cover 60`.
- **Por locus:** minimap2 2.31 `-c -x splice:hq` de GITV01 al genoma. Se acepta un transcripto con identidad
  matches/(qe − qs) ≥ 0.90, ≥ 50% de la query alineada y ≥ 100 pb sobre exones del gen.
- `soporte_GITV01` = RBH o locus.

### 6.5 Confianza (tarea H)

- alta = lecturas + GITV01
- media = solo lecturas
- media_revisar = solo GITV01
- baja = ninguna

### 6.6 Parálogos (tareas C y F)

- `paralogo_resuelto` (genes con 2 copias): cuántas copias tienen `rnaseq_soporte`.
- **Loci** (`r6_loci_todos.sbatch`): minimap2 `-c -x splice:hq --secondary=yes -N 20 -p 0.5` con todos los mRNA de
  Helixer. Se cuentan los loci no solapados con identidad ≥ 0.90 y cobertura ≥ 0.5.
- **Parálogos proteicos:** DIAMOND Dmai vs. Dmai `--more-sensitive -e 1e-10`; otro gen con ≥ 40% de identidad y ≥ 50% de
  cobertura de ambas proteínas → `paralogos_prot` y `max_ident_paralogo`.

### 6.7 Estabilidad (tarea E)

- `estabilidad_estadios` = CV_estadios.
- `estabilidad_condiciones` = CV entre las 4 condiciones colapsadas de PRJNA631706.
- Ninguna entra al score; un CV menor indica más estabilidad.

### 6.8 Candidatos finales y secuencia (tarea D; `scripts/verificar_candidatos.py`)

- **Finales** = los primeros 20 del ranking A que cumplen todo esto:
  - confianza distinta de baja
  - `paralogo_resuelto` distinto de ambas_con_lecturas
  - loci = 1
  - `max_ident_paralogo` < 70%
- **Comparación con StringTie** del mismo locus y hebra:
  - Se usan solo los intrones del **CDS** de Helixer. Un intrón está confirmado si aparece con exactamente el mismo donor
    y acceptor en algún transcripto de StringTie.
  - `cds_cubierto` = fracción del CDS de Helixer dentro de exones del transcripto elegido.
- **Estado:**
  - confirmado = todos los intrones y CDS ≥ 0.95
  - discrepancia = el resto
- **FASTA:** la secuencia de StringTie si cubre ≥ 90% del CDS (las lecturas mandan); si no, el ARNm de Helixer marcado
  `NO verificado`.

### 6.9 Competencia vectorial (tarea I)

- Tablas 3 y 4 de Xu et al. 2021 (DESeq2, p ajustado ≤ 0.1): 83 transcriptos anotados → `mrfv/mrfv_DE.tsv`.
- Mapeo de GITV01 a gen Helixer, en orden de prioridad: locus, RBH, mejor hit directo de proteína (≥ 80% de identidad,
  cobertura 60/60).
- Lista aparte, sin mezclar con el ranking, que marca letalidad_norm ≥ 0.75 y CV_estadios en el cuarto más bajo.
- Se declara el sesgo hacia genes conservados. No se reconstruye significancia a partir de la Tabla S7.

## 7. Sesgo de metamorfosis (`scripts/metamorfosis.py`)

### 7.1 Términos metamórficos

- Ontología: `go-basic.obo` releases/2026-07-26, la misma de la firma. No se volvió a descargar, para no desacoplar el
  flag de la firma.
- **Semillas** (las 7 resuelven al nombre esperado): GO:0007552, GO:0002165, GO:0007444, GO:0007476, GO:0035220,
  GO:0007562, GO:0035075.
- **Expansión:** descendientes por is_a y part_of, más los términos no obsoletos cuyo nombre contiene (sin distinguir
  mayúsculas, desde límite de palabra) metamorphosis, pupa, pupal, puparium, imaginal disc, wing disc, eclosion, ecdysone,
  ecdysteroid, juvenile hormone o instar.
- **Exclusión** (`exclusiones_aplicadas.tsv`):
  - (a) nombres con molt, moult, ecdysis, chitin, cuticl o embryo
  - (b) nombres que aplican igual a un hemimetábolo (ecdysone, ecdysteroid, juvenile hormone, instar, larva) sin
    marcador metamórfico (metamorph, pupa, puparium, imaginal disc, wing disc, eclosion)
- **Conjuntos:**
  - estricto = sin (a) ni (b): 146 términos
  - amplio = sin (a): 206 términos
- **Flag del gen:** sus GO propagados intersecan el conjunto estricto. `flag_metamorfosis_amplio` usa el conjunto amplio.

### 7.2 Proteasomal

- sí = anotado a GO:0000502 *proteasome complex*
- parcial = solo a GO:0010498 *proteasomal protein catabolic process*

### 7.3 Auditoría de la firma

- Profundidad = camino mínimo a la raíz de su ontología por is_a y part_of.
- Genérico = profundidad ≤ 3 o anotado en > 25% del universo.
- Genes aportantes = genes de entrenamiento del fondo que portan el término.
- "Solo proteasomales" = todos los aportantes son proteasomales.

### 7.4 `letal_ambos_regimenes`

- dali = dali22 observado; si dali22 es imputado, dali11 observado.
- dapi = dapi11 observado.
- Alto = ≥ percentil 0.75 de los observados. Cortes alternativos: 0.50 y 0.90.
- Ranking filtrado = solo los True, en el mismo orden. Es diagnóstico.

## 8. Orden de regeneración

```
python3 scripts/expresion_matrices.py
python3 scripts/score_candidatos.py
python3 scripts/score_candidatos.py ogs3
python3 scripts/comparar_ogs3.py
python3 scripts/rnaseq_integracion.py
python3 scripts/verificar_candidatos.py
python3 scripts/metamorfosis.py > resultados/metamorfosis/log_metamorfosis.txt
bash scripts/armar_resultados.sh
```

Las entradas pesadas (BAM, conteos, PAF, cobertura) están en `rnaseq_insumos/` y `score_insumos/`, y vienen de los sbatch.

## 9. Referencias

- Ulrich J. et al. (2015) *BMC Genomics* 16:674 — screen iBeetle y genes de entrenamiento (PMC4559001).
- Knorr E. et al. (2018) *Sci Rep* 8:2061 — dre4, ncm, RpII140 (PMC5794766).
- Rodrigues T.B. et al. (2021) *Front Plant Sci* 12:728652 — ledprona, PSMB5 en *Leptinotarsa decemlineata*.
- Yao J. et al. (2013) *PLoS ONE* — V-ATPasa B y D en *Peregrinus maidis* (PMC3737204).
- Xu Y., Willman M.R., Todd J.C., Kim K.-H., Redinbaugh M.G., Stewart L.R. (2021) *Microbiol Spectrum*
  doi 10.1128/Spectrum.00612-21 — transcriptoma GITV01 y respuesta a MRFV.
- Herramientas:
  - Helixer v0.3.6
  - OrthoFinder 2.5.5
  - DIAMOND 2.2.8 (SDumont) / 2.1.18 (cluster QB)
  - HISAT2
  - samtools
  - bedtools
  - StringTie 3.0.3
  - minimap2 2.31
  - gffread
  - BUSCO 6.1.0
  - Salmon 1.10.3 (versión previa de la expresión)
