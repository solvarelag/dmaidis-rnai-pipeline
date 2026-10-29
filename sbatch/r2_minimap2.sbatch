#!/bin/bash
#SBATCH --job-name=r2_minimap2
#SBATCH --partition=lncc-h100
#SBATCH --account=PROYECTO
#SBATCH --reservation=RESERVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH -o /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.out
#SBATCH -e /scratch/PROYECTO/USUARIO/dmaidis/logs/%x_%j.err
#SBATCH --chdir=/scratch/PROYECTO/USUARIO/dmaidis
# brief_rnaseq: puente por locus GITV01 -> genoma, Tarea F (loci del top 100) y secuencias para Tarea D.
set -euo pipefail
date +"inicio %F %T"
R=/scratch/PROYECTO/USUARIO/dmaidis
SX="singularity exec -B /scratch:/scratch"
G=$R/data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna
W=$R/out/rnaseq_genoma; cd $W
$SX $R/sif/minimap2.sif minimap2 -t 32 -d genoma.mmi $G 2> /dev/null
# GITV01 completo al genoma (un mejor alineamiento por transcripto)
$SX $R/sif/minimap2.sif minimap2 -t 48 -x splice:hq -uf --secondary=no genoma.mmi \
    $R/data/Dalbulus_maidis_transcriptome_GITV01.fas > gitv01_genoma.paf 2> /dev/null &
# Tarea F: mRNA de Helixer del top 100, hasta 20 loci por secuencia
$SX $R/sif/gffread.sif gffread $R/out/helixer/Dmaidis_helixer.gff3 -g $G -w helixer_mrna.fa
awk 'NR==FNR{k[$1".1"]=1; next} /^>/{p=(substr($1,2) in k)} p' $R/out/top100_helixer_ids.txt helixer_mrna.fa > top100_mrna.fa
grep -c "^>" top100_mrna.fa
$SX $R/sif/minimap2.sif minimap2 -t 32 -x splice:hq -uf --secondary=yes -N 20 -p 0.5 genoma.mmi \
    top100_mrna.fa > top100_loci.paf 2> /dev/null
wait
# lecturas totales alineadas (diagnostico Tarea A)
$SX $R/sif/samtools.sif samtools idxstats rnaseq_merged.bam > idxstats.txt
# transcriptos StringTie como FASTA (Tarea D)
$SX $R/sif/gffread.sif gffread stringtie_denovo.gtf -g $G -w stringtie_tx.fa
wc -l gitv01_genoma.paf top100_loci.paf; grep -c "^>" stringtie_tx.fa
date +"fin %F %T"
