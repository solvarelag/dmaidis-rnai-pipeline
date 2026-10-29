#!/bin/bash
#SBATCH --job-name=r3_gitv01_paralogos
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH -o /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.out
#SBATCH -e /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis/out/rnaseq_genoma
# GITV01 -> genoma con -ub (orientacion desconocida) y busqueda de paralogos a nivel proteina Dmai vs Dmai.
set -euo pipefail
R=/scratch/PROYECTO/USUARIO/dmaidis
SX="singularity exec -B /scratch:/scratch"
$SX $R/sif/minimap2.sif minimap2 -t 48 -x splice:hq -ub --secondary=no genoma.mmi \
    $R/data/Dalbulus_maidis_transcriptome_GITV01.fas > gitv01_genoma_ub.paf 2> /dev/null &
P=$R/orthofinder/proteomas/Dalbulus_maidis.faa
$SX $R/sif/diamond.sif diamond makedb --in $P -d dmai --quiet
$SX $R/sif/diamond.sif diamond blastp -q $P -d dmai -o dmai_vs_dmai.tsv --more-sensitive --evalue 1e-10 \
    --max-target-seqs 50 --threads 48 --quiet --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore &
wait
wc -l gitv01_genoma_ub.paf dmai_vs_dmai.tsv
