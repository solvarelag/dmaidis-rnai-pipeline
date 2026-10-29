#!/bin/bash
#SBATCH --job-name=r5_cob_split
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=03:00:00
#SBATCH --no-requeue
#SBATCH -o /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.out
#SBATCH -e /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis/out/rnaseq_genoma
# Rehace la cobertura de exones de r1 con -split: sin -split, el intron (N del CIGAR) de una lectura empalmada contaba
# como cobertura y daba soporte falso a genes dentro de intrones de otros genes. Un proceso por cromosoma.
set -euo pipefail
R=/scratch/PROYECTO/USUARIO/dmaidis; SX="singularity exec -B /scratch:/scratch"
mkdir -p cob_split
for ch in $(cut -f1 exones.s.bed | uniq); do
  ( awk -v c=$ch '$1==c' exones.s.bed > cob_split/$ch.bed
    $SX $R/sif/samtools.sif samtools view -u rnaseq_merged.bam $ch \
      | $SX $R/sif/bedtools.sif bedtools coverage -split -sorted -g genoma.txt -a cob_split/$ch.bed -b stdin > cob_split/$ch.cov ) &
  while [ $(jobs -r | wc -l) -ge 90 ]; do sleep 2; done
done
wait
for ch in $(cut -f1 exones.s.bed | uniq); do cat cob_split/$ch.cov; done > cobertura_exones_split.txt
wc -l exones.s.bed cobertura_exones_split.txt
rm -rf cob_split
date +"fin %F %T"
