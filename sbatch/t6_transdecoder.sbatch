#!/bin/bash
#SBATCH --job-name=t6_transdecoder
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
set -euo pipefail
date +"inicio %F %T"
SX="singularity exec -B /scratch:/scratch"
N=$SLURM_CPUS_PER_TASK
mkdir -p out/transdecoder && cd out/transdecoder
T=../../data/Dalbulus_maidis_transcriptome_GITV01.fas
$SX ../../sif/transdecoder.sif TransDecoder.LongOrfs -t $T
$SX ../../sif/transdecoder.sif TransDecoder.Predict -t $T --single_best_only
P=$(ls *.transdecoder.pep); echo "proteinas: $(grep -c '^>' $P)"
cd ../..
$SX sif/busco.sif busco -i out/transdecoder/$P -l busco_downloads/lineages/insecta_odb10 -m proteins -c $N \
    --offline --download_path busco_downloads --out_path out -o busco_Dmaidis_transdecoder_insecta -f > /dev/null
grep -m1 'C:' out/busco_Dmaidis_transdecoder_insecta/short_summary*.txt
date +"fin %F %T"
