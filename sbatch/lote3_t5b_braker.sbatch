#!/bin/bash
#SBATCH --job-name=lote3_t5b_braker
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=4-00:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# T5b: BRAKER3 (RNA-seq + OrthoDB12 Arthropoda) sobre el genoma enmascarado, y BUSCO del resultado.
set -euo pipefail
date +"inicio %F %T"
R=/scratch/PROYECTO/USUARIO/dmaidis
SX="singularity exec -B /scratch:/scratch"
N=$SLURM_CPUS_PER_TASK
# AUGUSTUS necesita una copia escribible de su config
if [ ! -d out/augustus_config ]; then
  $SX sif/braker3.sif bash -c 'cp -r ${AUGUSTUS_CONFIG_PATH:-/usr/share/augustus/config} '$R'/out/augustus_config'
fi
export SINGULARITYENV_AUGUSTUS_CONFIG_PATH=$R/out/augustus_config
$SX sif/braker3.sif braker.pl --genome=out/repeats/rmask/genome.fna.masked \
    --bam=out/hisat2/rnaseq.bam --prot_seq=data/Arthropoda_odb12.fa \
    --species=Dalbulus_maidis_ihDalMaid1 --threads=$N --softmasking --gff3 \
    --workingdir=out/braker_dmaidis
B=out/braker_dmaidis
echo "genes: $(awk -F'\t' '$3=="gene"' $B/braker.gff3 | wc -l)  proteinas: $(grep -c '^>' $B/braker.aa)"
# una isoforma por gen (la mas larga); IDs de BRAKER: g123.t1
awk '/^>/{if(id)print id"\t"len"\t"seq; id=substr($1,2); seq=""; len=0; next}{seq=seq $0; len+=length($0)} END{print id"\t"len"\t"seq}' $B/braker.aa \
  | awk -F'\t' '{g=$1; sub(/\.t[0-9]+$/,"",g); if($2>L[g]){L[g]=$2; S[g]=$3}} END{for(g in S) print ">Dmai|"g"\n"S[g]}' > data/Dalbulus_maidis.faa
echo "una isoforma por gen: $(grep -c '^>' data/Dalbulus_maidis.faa)"
for L in insecta_odb10 hemiptera_odb10; do
  $SX sif/busco.sif busco -i data/Dalbulus_maidis.faa -l busco_downloads/lineages/$L -m proteins -c $N \
      --offline --download_path busco_downloads --out_path out -o busco_Dmaidis_braker_${L%_odb10} -f > /dev/null
  echo "BUSCO $L: $(grep -m1 'C:' out/busco_Dmaidis_braker_${L%_odb10}/short_summary*.txt)"
done
date +"fin %F %T"
