#!/bin/bash
#SBATCH --job-name=h3_orthofinder_msa
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
# OrthoFinder pasada 2: reusa los DIAMOND de la pasada 1 y arma arboles por MSA.
set -euo pipefail
date +"inicio %F %T"
ls proteomas/*.faa | wc -l
singularity exec -B /scratch:/scratch ../sif/orthofinder.sif orthofinder -b proteomas/OrthoFinder/Results_pasada1/WorkingDirectory/ -M msa \
    -t 96 -a 8 -s arbol_especies.nwk -n pasada2_msa
S=$(ls -d proteomas/OrthoFinder/Results_pasada1*/WorkingDirectory/OrthoFinder/Results_pasada2_msa* proteomas/OrthoFinder/Results_pasada2_msa* 2>/dev/null | head -1)/Comparative_Genomics_Statistics
cat $S/Statistics_PerSpecies.tsv | head -12
head -20 $S/Statistics_Overall.tsv
date +"fin %F %T"
