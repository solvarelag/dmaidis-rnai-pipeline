# Resultados del pipeline — *Dalbulus maidis*

Estado: **parcial** (tareas 1–3 terminadas; 6 plan B corriendo; 4, 5, 7 sin lanzar). Fecha: 2026-10-27.

Cómputo: SDumont II, partición `lncc-h100` con una reserva (nodo `sdumont2nd2015`, 96 CPUs, ~1 TB de RAM).
La partición asigna el nodo completo, así que T2, T3 y T6 se corrieron empaquetadas en un único job (`743259`, 64+24+8 CPUs).
Ninguna herramienta usa GPU. Contenedores de Singularity de depot.galaxyproject.org: BUSCO 6.1.0, CD-HIT 4.8.1, TransDecoder 5.7.1.
Directorio en el cluster: `/scratch/PROYECTO/USUARIO/dmaidis` (sin backup). Los scripts están en `sbatch_sdumont/`.

---

## 1. ¿Existe anotación oficial de *D. maidis*?

**No se encontró.**
- Ensembl Rapid Release (`ftp.ensembl.org/pub/rapid-release/species/Dalbulus_maidis/`) y Ensembl organisms
  (`ftp.ebi.ac.uk/pub/ensemblorganisms/Dalbulus_maidis/`) → 404.
- API de búsqueda de beta.ensembl.org → sin resultados.
- ENA: ningún *analysis* (anotación) para el taxid 74065. El BioProject `PRJEB96548` solo tiene lecturas de PacBio Revio (WGS) y Hi-C,
  **no tiene RNA-seq ni Iso-Seq** de Sanger.
- Portal de DToL (API `data_portal?search=Dalbulus`) → 0 resultados. El genoma viene de AEGIS/ToL; es posible que la API del portal
  no lo indexe, así que vale un vistazo manual en el navegador.
- NCBI: los tres ensamblados están sin anotar (ver `chequeos_minimos.md`).

→ Las tareas 4–5 **no** se cancelan.

## 2. BUSCO modo genoma — `GCA_975739785.1`

BUSCO 6.1.0, `insecta_odb10`, modo `euk_genome_min` (miniprot), 64 CPUs, 232 s.

| C | S | D | F | M | n |
|---|---|---|---|---|---|
| **99.6%** | 97.8% | **1.8%** | 0.2% | 0.2% | 1367 |

Ensamblado: 141 scaffolds, 500 contigs, 0.004% de gaps, contig N50 5 Mb.

**Veredicto:** genoma haploide limpio y casi completo (C > 90%, D < 10%) → **anotarlo (tareas 4–5)**.
Con esto, la hipótesis "972 Mb = dos haplotipos" queda descartada empíricamente, además de por los metadatos.

⚠ BUSCO advierte que 331 de 1361 (24.3%) de los completos tienen codones stop internos en las predicciones de miniprot
(columna `E:`). Esto no afecta C ni D, pero para un ensamblado HiFi es alto. Puede ser un artefacto de miniprot (alineamientos
con frameshift o exones mal delimitados) o errores de indel en el ensamblado. Hay que mirar si BRAKER produce muchas proteínas
truncadas; si pasa, conviene volver sobre esto.

## 3. *Macrosteles*: ¿haplotipos o isoformas?

### 3a. Colapso CD-HIT (c = 0.98, aS = 0.9) sobre proteomas ya filtrados a una isoforma por gen

| Especie | Proteínas antes | Después | % colapsado |
|---|---|---|---|
| *Macrosteles quadrilineatus* | 21,979 | 16,789 | **23.6%** |
| *Homalodisca vitripennis* | 19,904 | 18,012 | 9.5% |

### 3b. BUSCO de proteínas sobre los proteomas filtrados

| Especie | Linaje | C | S | D | F | M |
|---|---|---|---|---|---|---|
| *Macrosteles* | insecta_odb10 | 97.3% | 59.8% | **37.5%** | 0.6% | 2.1% |
| *Macrosteles* | hemiptera_odb10 | 97.4% | 59.3% | **38.0%** | 0.2% | 2.4% |
| *Homalodisca* | insecta_odb10 | 91.4% | 80.5% | 10.9% | 2.3% | 6.4% |
| *Homalodisca* | hemiptera_odb10 | 92.4% | 81.7% | 10.7% | 1.3% | 6.3% |

(Agregué `hemiptera_odb10` porque es el linaje del reporte de NCBI, así el D se puede comparar directamente con el 38.1% de NCBI.)

**Veredicto:** filtrar a una isoforma por gen **no** cambia el D de *Macrosteles* (38.1% → 38.0%), y casi una cuarta parte de
su proteoma tiene una copia ≥ 98% idéntica (2.5× lo de *Homalodisca*). Son **haplotipos sin colapsar**.
→ *Macrosteles* **no se usa para contar parálogos**. Sí sirve para la topología y la transferencia de ortología, con cautela.
Opcional: usar `out/cdhit/Macrosteles_quadrilineatus_98.faa` (colapsado) en OrthoFinder en vez del filtrado, para reducir
falsos 1:2. Ver Limitaciones.

*Homalodisca* colapsa 9.5% y tiene D ≈ 11%. Es compatible con duplicaciones reales más algo de artefacto; usable con cautela.

---

## 4–5. Enmascarado + BRAKER3

Sin lanzar todavía (pendiente de decisión; ver resumen de la sesión).

Evidencia de RNA-seq disponible en SRA para *D. maidis* (taxid 74065), por ejemplo `PRJNA631706` (16 corridas Illumina,
~20 M de lecturas cada una, 2021), `PRJNA318849`, `PRJNA1282137` (71 M de lecturas, 2025). **No verifiqué** cuál corresponde a
Martin et al. 2021.

## 6. Set de proteomas

Plan B de *D. maidis* (TransDecoder `--single_best_only` sobre GITV01, 56,117 transcriptos) corriendo. Filtrados hasta ahora (local, `proteomas/`):
*Macrosteles* 21,979 · *Homalodisca* 19,904 · *Nilaparvata* 2020 18,139.

## 7. OrthoFinder

Sin lanzar.

## Limitaciones

(se completa al final)
