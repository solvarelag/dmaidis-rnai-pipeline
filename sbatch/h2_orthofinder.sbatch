#!/bin/bash
#SBATCH --job-name=h2_orthofinder
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH -o /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.out
#SBATCH -e /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis/orthofinder
# OrthoFinder pasada 1 (sin MSA). La pasada 2 (-b ... -M msa) se lanza aparte.
set -euo pipefail
date +"inicio %F %T"
ls proteomas/*.faa | wc -l
singularity exec -B /scratch:/scratch ../sif/orthofinder.sif orthofinder -f proteomas/ -S diamond \
    -t 96 -a 8 -s arbol_especies.nwk -n pasada1
S=proteomas/OrthoFinder/Results_pasada1/Comparative_Genomics_Statistics
cat $S/Statistics_PerSpecies.tsv | head -12
head -20 $S/Statistics_Overall.tsv
date +"fin %F %T"
