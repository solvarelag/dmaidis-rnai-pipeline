#!/bin/bash
# armar_resultados.sh — copia las tablas finales y los reportes a resultados/ (correr desde la raiz del proyecto,
# despues de la regeneracion completa y de scripts/metamorfosis.py). Los documentos .md de resultados/ se escriben a mano
# y este script no los toca.
set -euo pipefail
R=resultados
mkdir -p $R/tablas/ranking $R/tablas/ranking_ogs3 $R/tablas/candidatos $R/tablas/expresion $R/tablas/gitv01_mrfv \
         $R/tablas/firma_go $R/logs $R/informes_previos $R/scripts
cp ranking/ranking_completo.tsv ranking/sensibilidad_sin_dali22.tsv $R/tablas/ranking/
cp ranking_ogs3/ranking_completo.tsv ranking_ogs3/top100_solo_ogs3.tsv $R/tablas/ranking_ogs3/
cp ranking/candidatos_20_verificacion.tsv candidatos_20_verificados.fna ranking/candidatos_competencia_vectorial.tsv \
   $R/tablas/candidatos/
cp expresion/*.tsv $R/tablas/expresion/
cp mapeo_helixer_GITV01.tsv mrfv/mrfv_DE.tsv mrfv/mrfv_DE_mapeado.tsv $R/tablas/gitv01_mrfv/
cp ranking/firma_go.tsv ranking/entrenamiento.tsv $R/tablas/firma_go/
cp set_entrenamiento_GO.txt $R/
cp ranking/reporte_score.txt ranking/reporte_rnaseq.txt ranking/reporte_candidatos20.txt $R/logs/
cp ranking_ogs3/reporte_score.txt $R/logs/reporte_score_ogs3.txt
cp ranking_ogs3/reporte_comparacion.txt $R/logs/reporte_comparacion_ogs3.txt
cp chequeos_minimos.md resultados_pipeline.md resultados_lote1.md resultados_score.md resultados_rnaseq.md \
   $R/informes_previos/
cp scripts/*.py scripts/*.sh $R/scripts/
cp sbatch_sdumont/*.sbatch $R/scripts/ 2>/dev/null || true
find $R -type f | sort
