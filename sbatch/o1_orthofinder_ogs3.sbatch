#!/bin/bash
#SBATCH --job-name=o1_of_ogs3
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH -o /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.out
#SBATCH -e /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis/orthofinder_ogs3
# OrthoFinder pasada 1 (sin MSA) con Tribolium OGS3 (IDs TC) en lugar de RefSeq; las otras 9 especies iguales.
# Ademas DIAMOND Dmai vs OGS3 para la columna de conservacion.
set -euo pipefail
date +"inicio %F %T"
ls proteomas/*.faa | wc -l
singularity exec -B /scratch:/scratch ../sif/orthofinder.sif orthofinder -f proteomas/ -S diamond \
    -t 96 -a 8 -s ../orthofinder/arbol_especies.nwk -n ogs3
singularity exec -B /scratch:/scratch ../sif/diamond.sif diamond makedb --in proteomas/Tribolium_castaneum.faa -d tcas_ogs3 --quiet
singularity exec -B /scratch:/scratch ../sif/diamond.sif diamond blastp -q proteomas/Dalbulus_maidis.faa -d tcas_ogs3 \
    --more-sensitive -e 1e-5 -k 25 -p 96 --quiet -o dmai_vs_ogs3.tsv
S=proteomas/OrthoFinder/Results_ogs3/Comparative_Genomics_Statistics
head -12 $S/Statistics_PerSpecies.tsv
date +"fin %F %T"
