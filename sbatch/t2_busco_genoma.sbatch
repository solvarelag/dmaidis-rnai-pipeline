#!/bin/bash
#SBATCH --job-name=t2_busco_genoma
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=16:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
set -euo pipefail
date +"inicio %F %T"
SX="singularity exec -B /scratch:/scratch"
N=$SLURM_CPUS_PER_TASK
G=data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna
$SX sif/busco.sif busco -i $G -l busco_downloads/lineages/insecta_odb10 -m genome -c $N \
    --offline --download_path busco_downloads --out_path out -o busco_Dmaidis_genoma_insecta -f
cat out/busco_Dmaidis_genoma_insecta/short_summary*.txt
date +"fin %F %T"
