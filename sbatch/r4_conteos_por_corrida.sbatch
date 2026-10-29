#!/bin/bash
#SBATCH --job-name=r4_conteos
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
# Fragmentos por gen y por corrida desde los BAM de HISAT2 (reemplaza a salmon-sobre-CDS, que da ceros espurios).
# Fragmento = read1 (0x40), primario, mapeado, MAPQ >= 10 (unico); cuenta 1 vez por gen si algun bloque cae en un exon (-split).
set -euo pipefail
R=/scratch/PROYECTO/USUARIO/dmaidis; SX="singularity exec -B /scratch:/scratch"
mkdir -p conteos
$SX $R/sif/samtools.sif samtools idxstats bam/SRR5010611.bam | awk '$1!="*"{print $1"\t"$2}' > genoma_bam.txt
for b in bam/*.bam; do
  r=$(basename $b .bam)
  ( $SX $R/sif/samtools.sif samtools view -u -q 10 -f 0x40 -F 0x904 $b \
    | $SX $R/sif/bedtools.sif bedtools coverage -split -counts -sorted -g genoma_bam.txt -a genes.bed12 -b stdin \
    | cut -f4,13 > conteos/$r.tsv ) &
done
wait
wc -l conteos/*.tsv | tail -1
date +"fin %F %T"
