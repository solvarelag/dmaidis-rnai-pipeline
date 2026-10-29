#!/bin/bash
#SBATCH --job-name=t3_macrosteles
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
set -euo pipefail
date +"inicio %F %T"
SX="singularity exec -B /scratch:/scratch"
N=$SLURM_CPUS_PER_TASK
mkdir -p out/cdhit
for SP in Macrosteles_quadrilineatus Homalodisca_vitripennis; do
  $SX sif/cd-hit.sif cd-hit -i data/$SP.faa -o out/cdhit/${SP}_98.faa \
      -c 0.98 -n 5 -aS 0.9 -T $N -M 64000 > out/cdhit/${SP}.log
  a=$(grep -c "^>" data/$SP.faa); b=$(grep -c "^>" out/cdhit/${SP}_98.faa)
  echo -e "CDHIT\t$SP\t$a\t$b"
  for L in insecta_odb10 hemiptera_odb10; do
    $SX sif/busco.sif busco -i data/$SP.faa -l busco_downloads/lineages/$L -m proteins -c $N \
        --offline --download_path busco_downloads --out_path out -o busco_${SP}_filtrado_${L%_odb10} -f > /dev/null
    echo "BUSCO $SP $L: $(grep -m1 'C:' out/busco_${SP}_filtrado_${L%_odb10}/short_summary*.txt)"
  done
done
date +"fin %F %T"
