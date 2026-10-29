#!/bin/bash
#SBATCH --job-name=r1_rnaseq_genoma
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
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# brief_rnaseq Tarea A (+ insumos de B, C, D): HISAT2 por corrida al genoma, cobertura por gen, StringTie de novo.
set -euo pipefail
date +"inicio %F %T"
R=/scratch/PROYECTO/USUARIO/dmaidis
SX="singularity exec -B /scratch:/scratch"
G=$R/data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna
H=$R/out/helixer/Dmaidis_helixer.gff3
W=$R/out/rnaseq_genoma; mkdir -p $W/fq $W/bam && cd $W

# indice (32 hilos) en paralelo con la extraccion de FASTQ (8 x 8 hilos)
[ -s idx.8.ht2 ] || $SX $R/sif/hisat2.sif hisat2-build -p 32 $G idx > build.log 2>&1 &
pb=$!
ls $R/data/rnaseq/*.sra | xargs -P8 -I{} sh -c '[ -s fq/$(basename {} .sra)_2.fastq ] || '"$SX"' '$R'/sif/sra-tools.sif fasterq-dump --split-files -e 8 -t fq/tmp_$(basename {} .sra) -O fq {} > /dev/null 2>&1'
wait $pb
date +"indice y fastq listos %F %T"

alinear() {
  r=$1
  [ -s bam/$r.bam.bai ] && return 0
  $SX $R/sif/hisat2.sif hisat2 -p 20 --dta -x idx -1 fq/${r}_1.fastq -2 fq/${r}_2.fastq 2> bam/$r.hisat2.log \
    | $SX $R/sif/samtools.sif samtools sort -@ 4 -m 3G -T bam/tmp_$r -o bam/$r.bam
  $SX $R/sif/samtools.sif samtools index bam/$r.bam
  rm -f fq/${r}_1.fastq fq/${r}_2.fastq
  echo -e "$r\t$(tail -n 1 bam/$r.hisat2.log)"
}
export -f alinear; export SX R
ls $R/data/rnaseq/*.sra | xargs -n1 basename | sed 's/.sra$//' | xargs -P4 -I{} bash -c 'alinear {}' > alineamiento_por_corrida.tsv
cat alineamiento_por_corrida.tsv
$SX $R/sif/samtools.sif samtools merge -@ 32 -f rnaseq_merged.bam bam/*.bam
$SX $R/sif/samtools.sif samtools index -@ 16 rnaseq_merged.bam

# StringTie de novo (sin -G) sobre todo el genoma, en paralelo con la cobertura
$SX $R/sif/stringtie.sif stringtie rnaseq_merged.bam -p 32 -o stringtie_denovo.gtf -l STRG > stringtie.log 2>&1 &
ps=$!

# cobertura: exones por gen (fraccion cubierta) y span del gen (n de lecturas)
$SX $R/sif/samtools.sif samtools faidx $G 2>/dev/null || true
cut -f1,2 $G.fai > genoma.txt
awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon"{match($9,/Parent=[^;]+/); p=substr($9,RSTART+7,RLENGTH-7); sub(/\.[0-9]+$/,"",p); print $1,$4-1,$5,p}' $H \
  | sort -k1,1 -k2,2n > exones.bed
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{match($9,/ID=[^;]+/); print $1,$4-1,$5,substr($9,RSTART+3,RLENGTH-3)}' $H \
  | sort -k1,1 -k2,2n > genes.bed
# ordenar el bed segun el orden del genoma (requisito de -sorted)
for b in exones genes; do
  awk 'NR==FNR{o[$1]=NR;next}{print o[$1]"\t"$0}' genoma.txt $b.bed | sort -k1,1n -k3,3n | cut -f2- > $b.s.bed
done
$SX $R/sif/bedtools.sif bedtools coverage -sorted -g genoma.txt -a exones.s.bed -b rnaseq_merged.bam > cobertura_exones.txt &
$SX $R/sif/bedtools.sif bedtools coverage -sorted -g genoma.txt -counts -a genes.s.bed -b rnaseq_merged.bam > lecturas_genes.txt &
wait
wc -l cobertura_exones.txt lecturas_genes.txt stringtie_denovo.gtf
rm -rf fq
date +"fin %F %T"
