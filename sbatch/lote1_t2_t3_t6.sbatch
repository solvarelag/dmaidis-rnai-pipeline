#!/bin/bash
#SBATCH --job-name=lote1_t2_t3_t6
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=16:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# Empaqueta tres tareas independientes en el nodo exclusivo (64+24+8 = 96 CPUs).
set -uo pipefail
date +"inicio %F %T"
SLURM_CPUS_PER_TASK=64 bash sbatch/t2_busco_genoma.sbatch > logs/t2_${SLURM_JOB_ID}.log 2>&1 & p2=$!
SLURM_CPUS_PER_TASK=24 bash sbatch/t3_macrosteles.sbatch  > logs/t3_${SLURM_JOB_ID}.log 2>&1 & p3=$!
SLURM_CPUS_PER_TASK=8  bash sbatch/t6_transdecoder.sbatch > logs/t6_${SLURM_JOB_ID}.log 2>&1 & p6=$!
fail=0
for p in $p2 $p3 $p6; do wait $p || { echo "pid $p fallo"; fail=1; }; done
date +"fin %F %T"
exit $fail
