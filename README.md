# Priorización de blancos de RNAi en *Dalbulus maidis*

Pipeline para elegir genes blanco de RNAi en la chicharrita del maíz, *Dalbulus maidis* (Hemiptera: Cicadellidae). La
letalidad medida por RNAi en *Tribolium castaneum* (screen **iBeetle**) se transfiere por ortología a un proteoma de
*D. maidis* anotado con Helixer. Los candidatos se priorizan con un score que combina letalidad, una firma GO de blancos
validados, expresión constitutiva y conservación. Después se validan con RNA-seq: soporte del modelo génico, parálogos,
loci y secuencia verificada para diseñar el dsRNA.

> **Resultados:** empezar por [`resultados/LEEME.md`](resultados/LEEME.md), que explica cómo leer cada archivo, y por
> [`resultados/RESULTADOS.md`](resultados/RESULTADOS.md). El método completo está en
> [`resultados/METODOLOGIA.md`](resultados/METODOLOGIA.md).

## Qué produce

- **Ranking oficial (A)**, 7,881 genes:
  `score = 0.35·letalidad + 0.30·firma GO + 0.20·expresión constitutiva + 0.15·conservación` (rango percentil).
- **Ranking B**, sin firma GO, para controlar la circularidad del set de entrenamiento.
- **Análisis separado con *Tribolium* OGS3**, para descartar pérdidas por el mapeo de IDs.
- **20 candidatos finales** con secuencia verificada por RNA-seq (`candidatos_20_verificados.fna`).
- **Diagnósticos de sesgo:**
  - letalidad sin `dali22`
  - flag de metamorfosis (holometábolo → hemimetábolo)
  - `letal_ambos_regimenes` (inyección larval y pupal)
  - auditoría de la firma GO
- **Lista aparte de competencia vectorial** (respuesta a MRFV).

## Estructura del repositorio

```
├── README.md
├── requirements.txt
├── arbol_especies.nwk        árbol de especies fijo para OrthoFinder
├── preparacion/              búsqueda y descarga de proteomas de NCBI; una isoforma por gen
├── sbatch/                   jobs Slurm (cómputo pesado: Helixer, OrthoFinder, HISAT2, StringTie, minimap2, DIAMOND)
├── scripts/                  análisis en Python (score, integración RNA-seq, verificación, metamorfosis)
└── resultados/               tablas finales, documentos y logs de la última corrida
```

Los datos crudos e intermedios (genoma, BAM, GFF, proteomas, tablas de iBeetle, transcriptoma GITV01) **no** están en el
repositorio por tamaño y licencia. Abajo se indica de dónde sale cada uno.

## Requisitos

- **Python ≥ 3.10** con `pandas`, `numpy` y `scipy` (`pip install -r requirements.txt`).
- **Cluster con Slurm** y contenedores **Singularity/Apptainer**, que los sbatch esperan en `<BASE>/sif/`:
  - Helixer v0.3.6 (GPU)
  - OrthoFinder 2.5.5
  - DIAMOND 2.2.8
  - HISAT2
  - samtools
  - bedtools
  - StringTie 3.0.3
  - minimap2 2.31
  - gffread
  - BUSCO 6.1.0
  - sra-tools
  - salmon 1.10.3 (opcional)
- `ncbi-datasets-cli`, para bajar proteomas (`preparacion/descargar_proteomas.sh`).

Los sbatch se escribieron para SDumont II (LNCC; partición `lncc-h100`, 96 CPU por nodo). Antes de usarlos hay que
reemplazar los marcadores `PROYECTO`, `USUARIO` y `RESERVA`, y la ruta base `/scratch/PROYECTO/USUARIO/dmaidis`. Si no hay
reserva, se borra la línea `#SBATCH --reservation`.

## Datos de entrada

| Dato | Fuente | Dónde lo espera el pipeline |
|---|---|---|
| Genoma *D. maidis* GCA_975739785.1 (pseudohaplotipo principal) | NCBI Datasets | `<BASE>/data/Dmaidis_genome/…` (sbatch) |
| Proteomas de las otras 9 especies (RefSeq) | `preparacion/descargar_proteomas.sh` + `una_isoforma_por_gen.py` | `proteomas/*.faa` |
| Letalidad iBeetle (`iB;dali11;dali22;dapi11;phase`) | iBeetle-Base (https://ibeetle-base.uni-goettingen.de/) | `tabla_maestra/iBeetle_lethality.csv` |
| iB → TC | iBeetle-Base | `tabla_maestra/iB_TC.csv` |
| Proteoma *Tribolium* OGS3 | iBeetle-Base (`download/species/Tcas/OGS3_proteins.fasta.gz`) | `tabla_maestra/OGS3_proteins.fasta.gz` |
| GeneID ↔ TC de *Tribolium* | NCBI `gene_info` (Invertebrates), locus_tag | `tabla_maestra/Tcas_GeneID_TC.tsv` |
| GeneID ↔ FBgn y símbolos de *Drosophila* | NCBI `gene_info` (dbXrefs) | `tabla_maestra/Dmel_GeneID_FBgn.tsv`, `Dmel_FBgn_simbolo.tsv` |
| GO de *Drosophila* (`id, GO, aspecto`) | anotaciones GO de FlyBase (GAF) | `tabla_maestra/Dmel_FBgn_GO.tsv` |
| Ontología GO | `go-basic.obo` (release usada: 2026-07-26) | `tabla_maestra/go-basic.obo` |
| RNA-seq | SRA: PRJNA318849 (ninfa, hembra, macho) y PRJNA631706 (MRFV, 16 corridas) | se baja en `sbatch/r1_rnaseq_genoma.sbatch` |
| Transcriptoma GITV01 | NCBI TSA GITV01 (Xu et al. 2021) | `<BASE>/data/Dalbulus_maidis_transcriptome_GITV01.fas` |
| Transcriptos DE de MRFV | Tablas 3 y 4 de Xu et al. 2021 | `mrfv/mrfv_DE.tsv` (está en `resultados/tablas/gitv01_mrfv/`) |

## Cómo se usa

Los scripts de Python se corren desde la raíz del directorio de trabajo, con las carpetas de entrada de la tabla anterior.
Las salidas del cluster se copian a `score_insumos/` y `rnaseq_insumos/`.

### 1. Anotación y ortología (cluster)

| Paso | Job | Salida principal |
|---|---|---|
| BUSCO del genoma | `sbatch/t2_busco_genoma.sbatch` | completitud del genoma |
| Proteoma de *D. maidis* | `sbatch/h1_helixer.sbatch` (4 GPU) | `Dalbulus_maidis_helixer.faa`, `Dmaidis_helixer.gff3` |
| Ortología, 10 especies | `sbatch/h2_orthofinder.sbatch` | `Orthologues_Dalbulus_maidis/…`, `Tribolium_castaneum__v__Drosophila_melanogaster.tsv` → `orthofinder_pasada1/` |
| Ortología con OGS3 (separada) | `sbatch/o1_orthofinder_ogs3.sbatch` | → `orthofinder_ogs3/` + `score_insumos/dmai_vs_ogs3.tsv` |
| Mapeo OGS3↔RefSeq e identidad Dmai–Tcas | `sbatch/s4_expr_rbh_ident.sbatch` | `score_insumos/` → `python3 scripts/rbh_ogs3_refseq.py` |

### 2. RNA-seq (cluster)

| Paso | Job | Salida |
|---|---|---|
| HISAT2 por corrida, merge y StringTie de novo | `sbatch/r1_rnaseq_genoma.sbatch` | BAM, `stringtie_denovo.gtf` |
| GITV01 al genoma | `sbatch/r2_minimap2.sbatch` | `gitv01_genoma_c.paf`, `stringtie_tx.fa` |
| DIAMOND Dmai vs. Dmai | `sbatch/r3_gitv01_paralogos.sbatch` | `dmai_vs_dmai.tsv` |
| Fragmentos por gen y corrida | `sbatch/r4_conteos_por_corrida.sbatch` | `conteos/SRR*.tsv` |
| Cobertura exónica con `-split` | `sbatch/r5_cobertura_split.sbatch` | `cobertura_exones_split.txt` |
| Loci de todos los mRNA | `sbatch/r6_loci_todos.sbatch` | `todos_loci_c.paf` |

El puente con GITV01 se arma con `python3 scripts/puente_gitv01.py` (RBH de proteína con DIAMOND), que genera
`mapeo_helixer_GITV01.tsv`.

### 3. Score, validación y diagnósticos (local)

Hay que correrlos **en este orden**: cada paso usa la salida del anterior.

```bash
python3 scripts/expresion_matrices.py        # TPM por BioProject -> expresion/
python3 scripts/score_candidatos.py          # ranking A -> ranking/
python3 scripts/score_candidatos.py ogs3     # análisis separado -> ranking_ogs3/
python3 scripts/comparar_ogs3.py             # RefSeq vs. OGS3
python3 scripts/rnaseq_integracion.py        # soporte, confianza, parálogos, loci, estabilidad, MRFV
python3 scripts/verificar_candidatos.py      # 20 finales + candidatos_20_verificados.fna
python3 scripts/metamorfosis.py > resultados/metamorfosis/log_metamorfosis.txt
bash scripts/armar_resultados.sh             # copia todo a resultados/
```

Opcional: `python3 scripts/conteo_letales.py` (conteo de letales con dali11 ≥ 80 y sus copias en *D. maidis*).

## Resultados principales (corrida del 2026-10-28)

- **Proteoma Helixer:** 22,199 genes, BUSCO insecta C 93.0%.
- **Ranking A:** el top 10 tiene 6 subunidades del proteasoma (Rpt1, Rpt5, Rpt2, Rpt6, Rpn11, Rpn6). Los controles
  positivos quedan arriba: PSMB5 (blanco de ledprona) #25, Rpt3 #21, V-ATPasa D #125, V-ATPasa B #216.
- **Sin firma GO (ranking B):** comparte 63/100 del top 100 y el proteasoma sigue con 15 genes en su top 100.
- **Sesgo de metamorfosis:** 0 de los 77 términos de la firma son metamórficos. El flag marca 15 del top 100, todos genes
  esenciales pleiotrópicos.
- **RNA-seq:** 100/100 del top 100 con soporte de lecturas y un solo locus. De los 20 finales, 9 tienen todos los intrones
  confirmados y 11 discrepan con Helixer, así que para el dsRNA se usa la secuencia de StringTie.

El detalle y las limitaciones están en [`resultados/RESULTADOS.md`](resultados/RESULTADOS.md).

## Referencias

- Ulrich J. et al. (2015) *BMC Genomics* 16:674 — screen iBeetle.
- Knorr E. et al. (2018) *Sci Rep* 8:2061.
- Rodrigues T.B. et al. (2021) *Front Plant Sci* 12:728652 — ledprona (PSMB5).
- Yao J. et al. (2013) *PLoS ONE* — V-ATPasa en *Peregrinus maidis*.
- Xu Y. et al. (2021) *Microbiol Spectrum*, doi 10.1128/Spectrum.00612-21 — GITV01 y MRFV.
