#!/bin/bash
#SBATCH --job-name=r6_loci_todos
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
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# Tarea F para TODOS los mRNA de Helixer (el ranking cambio y top100_loci.paf solo cubria el top 100 anterior).
set -euo pipefail
R=/scratch/PROYECTO/USUARIO/dmaidis; SX="singularity exec -B /scratch:/scratch"
cd $R/out/rnaseq_genoma
$SX $R/sif/minimap2.sif minimap2 -t 96 -c -x splice:hq -uf --secondary=yes -N 20 -p 0.5 genoma.mmi \
    helixer_mrna.fa > todos_loci_c.paf 2> /dev/null
wc -l todos_loci_c.paf
date +"fin %F %T"
