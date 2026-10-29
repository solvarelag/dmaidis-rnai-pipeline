# Chequeos mínimos previos al pipeline de ortología

Fecha: 2026-10-27. Fuente: API REST de NCBI Datasets v2 (`/genome/accession/.../dataset_report` y `/download`).
El CLI `datasets` no estaba instalado; todo se hizo con `curl` + `jq`. La v2alpha también respondió (200), pero se usó la v2.

---

## 1. Ensamblados de *Dalbulus maidis*

| Campo | GCA_975739785.1 | GCA_975739835.1 | GCA_054756405.1 |
|---|---|---|---|
| assembly_name | ihDalMaid1.1 | ihDalMaid1.1 alternate haplotype | CIAP v1 |
| assembly_type | haploid | **alternate-pseudohaplotype** | haploid |
| diploid_role | **principal_pseudohaplotype_of_diploid** | **alternate_pseudohaplotype_of_diploid** | — |
| linked_assemblies | GCA_975739835.1 (alternate) | GCA_975739785.1 (principal) | — |
| Nivel | Chromosome | Scaffold | Contig |
| BioProject | PRJEB97791 | PRJEB97790 | PRJNA1120083 |
| BioProject lineage | PRJEB96548 (*D. maidis*) → Sanger ToL / AEGIS / EBP | ídem | — |
| BioSample | SAMEA117578701 | SAMEA117578701 (el mismo) | SAMN41674197 |
| Submitter | Wellcome Sanger Institute | Wellcome Sanger Institute | Instituto de Patología Vegetal CIAP-INTA |
| atypical | (ausente) | (ausente) | (ausente) |
| infraspecific_names | (ausente) | (ausente) | isolate: CIAP |
| Scaffolds | 140 | 7,864 | 22,333 |
| Scaffold N50 | 93.7 Mb | 493 kb | 50.5 kb |
| Longitud total | 972.2 Mb | 488.5 Mb | 585.2 Mb |
| Anotación | no | no | no |

**Conclusión:** `GCA_975739785.1` es el **pseudohaplotipo principal** (un solo juego cromosómico, 972 Mb) de un individuo diploide.
`GCA_975739835.1` es su **haplotipo alternativo** (mismo BioSample, parcial y fragmentado). Es decir, 972 Mb **no** es "2 × 488":
el alternativo solo trae las regiones que se pudieron separar en fase.
→ Para la arquitectura del pipeline, el ensamblado a anotar (BRAKER3) es `GCA_975739785.1` solo; el alternativo no se debe
sumar al proteoma. (El ensamblado CIAP, de 585 Mb, es de otro individuo/grupo y su tamaño difiere bastante del de Sanger;
no se investigó por qué.)

Metadatos de los cicadélidos/delfácidos RefSeq: todos `assembly_type: haploid`, sin `diploid_role` ni ensamblado alternativo
vinculado, así que los metadatos no informan ploidía. *Macrosteles* (`isolate MER2022`) y *Homalodisca* (`AUS2020`) dicen
`sex: pooled male and female`, o sea que **se secuenciaron varios individuos juntos**. Eso aumenta la heterocigosis efectiva
y el riesgo de que haya haplotipos sin colapsar. *Nilaparvata* 2014: `isolate NLH13, sex female`.

---

## 2. Ratio proteínas/genes (pregunta 2)

Genes únicos = `GeneID` distintos en las líneas CDS del GFF. Genes con más de una isoforma = genes con más de un `protein_id` en CDS.

| Especie | Accession | Proteínas | Genes únicos | Ratio | Genes multi-isoforma | Genes prot.-coding NCBI | Ensamblado | BUSCO D (NCBI) | Veredicto |
|---|---|---|---|---|---|---|---|---|---|
| *Macrosteles quadrilineatus* | GCF_028750875.1 | 34,269 | 21,979 | **1.55** | 25.6% | 21,979 ✔ | 1.32 Gb, 1,164 scaff. | **38.1%** | D alto **no** explicado por isoformas (ver abajo) → probable duplicación haplotípica |
| *Homalodisca vitripennis* | GCF_021130785.1 | 31,163 | 19,904 | 1.56 | 25.0% | 19,904 ✔ | 2.30 Gb, 14,903 scaff. | 10.8% | Isoformas presentes; D moderado |
| *Nilaparvata lugens* (2014) | GCF_000757685.1 | 24,319 | 19,819 | 1.22 | 12.8% | 19,819 ✔ | 1.14 Gb, 46,558 scaff. | 9.6% | Pocas isoformas; D moderado |

El conteo de genes coincide exactamente con el de NCBI en las tres especies, así que está bien hecho.

**Interpretación (ojo: el criterio del brief no alcanza para decidir en este caso):**

- Siguiendo el criterio literal del brief, *Macrosteles* tiene ratio > 1.5, lo que daría "isoformas → filtrar baja el D".
  **Pero** *Homalodisca* tiene **el mismo ratio (1.56) y la misma fracción de genes multi-isoforma (25%)** con un D de solo 10.8%.
  Si las isoformas explicaran el 38%, *Homalodisca* también tendría el D alto. Por lo tanto, el exceso de D de
  *Macrosteles* (~27 puntos por encima de los otros dos) no se explica por isoformas.
- Además, según la documentación del pipeline de anotación de NCBI (EGAP), el BUSCO del reporte se calcula sobre la
  proteína más larga de cada gen. Si es así, ese 38% ya está libre de isoformas. **No lo verifiqué en esta sesión**, es
  conocimiento previo; si importa, conviene confirmarlo en la página de la anotación.
- El número de genes de *Macrosteles* (21,979) no es extremadamente alto, así que el dato "adicional gratis" del brief
  (>20k genes) es débil. El ensamblado pooled y el BUSCO D de 38% son más indicativos.

**Veredicto por especie:**
- *Macrosteles*: **probablemente haplotipos sin colapsar** → no es confiable para contar parálogos, aunque se filtre a una
  isoforma por gen. Hay que confirmarlo con BUSCO en modo proteína sobre el proteoma filtrado (pendiente 2).
- *Homalodisca* y *Nilaparvata* 2014: filtrar a una isoforma por gen es necesario para OrthoFinder; el D de ~10% es
  compatible con duplicación génica real o artefacto menor. Usables, con cautela.

---

## 3. *Nilaparvata lugens* 2020 — GCF_014356525.2

Confirmado que tiene anotación y que es descargable (`PROT_FASTA` + `GENOME_GFF`, ~15 MB comprimidos):

- ASM1435652v1, nivel Chromosome, 1.09 Gb, 6,874 scaffolds, scaffold N50 77.6 Mb, haploid, Huazhong Agricultural University, 2020-09-09.
- Anotación RefSeq Gnomon: 20,891 genes (18,139 codificantes de proteína).
- `protein.faa`: **32,453 proteínas** sobre 18,139 genes, ratio **1.78** (29.9% de los genes con más de una isoforma) → hay
  que filtrar a una isoforma por gen.
- NCBI **no** reporta BUSCO para esta anotación (campo `busco` nulo), así que queda en el pendiente 6.

---

## 4. Limitaciones y chequeos pendientes (no ejecutados por falta de cómputo)

| # | Chequeo | Qué decisión del pipeline quedaría mejor fundada |
|---|---|---|
| 1 | **BUSCO modo genoma** sobre `GCA_975739785.1` | La ploidía ya está resuelta por metadatos (pseudohaplotipo principal). Quedaría confirmada la completitud real del genoma y su D bajo, lo que justifica anotarlo con BRAKER3 en vez de usar el transcriptoma. *Nota:* en la carpeta hay un intento previo (`busco_39836.log`) que falló porque no existía `GCA_975739785.1_genomic.fna`; no llegó a correr. |
| 2 | **BUSCO modo proteína** sobre los proteomas filtrados (una isoforma por gen) | Confirmaría empíricamente que el D de *Macrosteles* sigue alto después del filtro, o sea que son haplotipos y que debe excluirse del conteo de parálogos. Es el test decisivo para la pregunta 2. |
| 3 | **`purge_dups`** sobre los ensamblados con duplicación haplotípica (*Macrosteles*) | Permitiría recuperar *Macrosteles* para el anillo 1 en vez de descartarlo. Requiere reanotar o reasignar genes, así que es caro. |
| 4 | **Anotación con BRAKER3** de *D. maidis* (`GCA_975739785.1`) usando el transcriptoma (56,116 transcriptos, Martin et al. 2021) como evidencia | Daría un proteoma de la especie blanco con un gen por locus. Sin esto, el conteo de parálogos en *D. maidis* depende de TransDecoder, que mezcla isoformas y alelos. |
| 5 | **Anotación de *Nephotettix cincticeps*** (`GCA_023375725.1`) | Reforzaría el anillo 1 (Deltocephalinae) si *Macrosteles* queda descartado. Hoy el anillo 1 depende de una especie dudosa. |
| 6 | **BUSCO local** de los ensamblados sin BUSCO de NCBI (Drosophila, *Nilaparvata* 2020) | Justificaría reemplazar *Nilaparvata* 2014 por la versión 2020 y confirmaría la calidad de las referencias externas. |

---

### Notas de reproducibilidad
- Los ZIPs descargados y `report_v2.json` quedaron en el scratchpad de la sesión (temporal), no en el proyecto.
- Comandos: `curl "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<ACC>/download?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF"`,
  `grep -c "^>" protein.faa`, y `awk '$3=="CDS"' genomic.gff | grep -oE "GeneID:[0-9]+" | sort -u | wc -l`.
