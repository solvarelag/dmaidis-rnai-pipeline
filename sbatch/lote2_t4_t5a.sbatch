#!/bin/bash
#SBATCH --job-name=lote2_t4_t5a
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=5-00:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# Empaqueta T4 (80 CPUs) y T5a (16 CPUs) en el nodo exclusivo.
set -uo pipefail
date +"inicio %F %T"
SLURM_CPUS_PER_TASK=80 bash sbatch/t4_repeats.sbatch > logs/t4_${SLURM_JOB_ID}.log 2>&1 & p4=$!
SLURM_CPUS_PER_TASK=16 bash sbatch/t5a_rnaseq.sbatch > logs/t5a_${SLURM_JOB_ID}.log 2>&1 & p5=$!
fail=0
for p in $p4 $p5; do wait $p || { echo "pid $p fallo"; fail=1; }; done
date +"fin %F %T"
exit $fail
