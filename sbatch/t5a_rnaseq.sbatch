#!/bin/bash
# T5a: SRA -> FASTQ -> HISAT2 -> BAM ordenado. No depende del enmascarado. Se corre desde lote2.
set -euo pipefail
date +"inicio %F %T"
SX="singularity exec -B /scratch:/scratch"
N=$SLURM_CPUS_PER_TASK
R=/scratch/PROYECTO/USUARIO/dmaidis
G=$R/data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna
mkdir -p $R/out/hisat2/fq $R/out/hisat2/tmp && cd $R/out/hisat2
for s in $R/data/rnaseq/*.sra; do
  r=$(basename $s .sra)
  [ -s fq/${r}_2.fastq ] || $SX $R/sif/sra-tools.sif fasterq-dump --split-files -e $N -t tmp -O fq $s
done
[ -s idx.1.ht2 ] || $SX $R/sif/hisat2.sif hisat2-build -p $N $G idx > build.log 2>&1
R1=$(ls fq/*_1.fastq | paste -sd,); R2=$(ls fq/*_2.fastq | paste -sd,)
$SX $R/sif/hisat2.sif hisat2 -p $N --dta -x idx -1 $R1 -2 $R2 2> hisat2.log \
  | $SX $R/sif/samtools.sif samtools sort -@ 4 -m 4G -T tmp/sort -o rnaseq.bam
$SX $R/sif/samtools.sif samtools index rnaseq.bam
tail -n 1 hisat2.log
rm -rf fq tmp
date +"fin %F %T"
