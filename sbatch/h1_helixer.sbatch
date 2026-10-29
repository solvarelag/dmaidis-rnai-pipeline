#!/bin/bash
#SBATCH --job-name=h1_helixer
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=10:00:00
#SBATCH --no-requeue
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# Proteoma rapido de D. maidis: Helixer (invertebrate) sobre el genoma, 4 bloques de cromosomas -> 4 GPUs.
# Termina con CHEQUEO 1 (BUSCO); sale con error si no pasa, para que OrthoFinder (afterok) no arranque.
set -euo pipefail
date +"inicio %F %T"
R=/scratch/PROYECTO/USUARIO/dmaidis
G=$R/data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna
SX="singularity exec -B /scratch:/scratch"
MODEL=$(find $R/helixer_models -name "invertebrate*.h5" | sort | tail -1); echo "modelo: $MODEL"
W=$R/out/helixer; mkdir -p $W && cd $W
# reparto goloso por longitud en 4 bloques
awk '/^>/{if(id)print id"\t"l; id=substr($1,2); l=0; next}{l+=length($0)} END{print id"\t"l}' $G | sort -k2,2nr \
  | awk '{b=0; for(i=1;i<4;i++) if(t[i]<t[b]) b=i; t[b]+=$2; print $1"\t"b} END{for(i=0;i<4;i++) printf "bloque %d: %.0f Mb\n",i,t[i]/1e6 > "/dev/stderr"}' > bloques.tsv
awk 'NR==FNR{b[$1]=$2; next} /^>/{id=substr($1,2); f="bloque_"b[id]".fa"} {print > f}' bloques.tsv $G
pids=()
for i in 0 1 2 3; do
  mkdir -p $W/tmp_$i
  CUDA_VISIBLE_DEVICES=$i $SX --nv $R/sif/helixer.sif Helixer.py --fasta-path bloque_$i.fa --lineage invertebrate \
      --model-filepath $MODEL --subsequence-length 213840 --overlap-offset 106920 --overlap-core-length 160380 --species Dalbulus_maidis --gff-output-path bloque_$i.gff3 \
      --temporary-dir $W/tmp_$i > helixer_$i.log 2>&1 & pids+=($!)
done
sleep 120; nvidia-smi || true
for p in "${pids[@]}"; do wait $p; done
cat bloque_[0-3].gff3 | grep -v "^##" | sed '1i ##gff-version 3' > Dmaidis_helixer.gff3
$SX $R/sif/gffread.sif gffread Dmaidis_helixer.gff3 -g $G -y Dmaidis_helixer.pep.raw -S
# una proteina por gen (Helixer da un transcripto por gen; se asegura igual), sin '.' de stop, >= 30 aa
awk '/^>/{if(id)print id"\t"seq; id=substr($1,2); seq=""; next}{seq=seq $0} END{print id"\t"seq}' Dmaidis_helixer.pep.raw \
  | awk -F'\t' '{g=$1; sub(/\.[0-9]+$/,"",g); gsub(/[.*]$/,"",$2); if(length($2)>length(S[g])) S[g]=$2} END{for(g in S) if(length(S[g])>=30) print ">Dmai|"g"\n"S[g]}' > $R/data/Dalbulus_maidis_helixer.faa
echo "genes: $(awk -F'\t' '$3=="gene"' Dmaidis_helixer.gff3 | wc -l)  proteinas finales: $(grep -c '^>' $R/data/Dalbulus_maidis_helixer.faa)"
cd $R
for L in insecta_odb10 hemiptera_odb10; do
  $SX sif/busco.sif busco -i data/Dalbulus_maidis_helixer.faa -l busco_downloads/lineages/$L -m proteins -c 96 \
      --offline --download_path busco_downloads --out_path out -o busco_Dmaidis_helixer_${L%_odb10} -f > /dev/null
  echo "BUSCO $L: $(grep -m1 'C:' out/busco_Dmaidis_helixer_${L%_odb10}/short_summary*.txt)"
done
# CHEQUEO 1: C > 85% y D < 15% (insecta) para seguir
read C D < <(grep -m1 'C:' out/busco_Dmaidis_helixer_insecta/short_summary*.txt | sed -E 's/.*C:([0-9.]+)%.*D:([0-9.]+)%.*/\1 \2/')
date +"fin %F %T"
if awk -v c=$C -v d=$D 'BEGIN{exit !(c>85 && d<15)}'; then
  cp data/Dalbulus_maidis_helixer.faa orthofinder/proteomas/Dalbulus_maidis.faa; echo "CHEQUEO1 OK C=$C D=$D"
else
  echo "CHEQUEO1 FALLA C=$C D=$D"; exit 3
fi
