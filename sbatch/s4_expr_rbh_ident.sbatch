#!/bin/bash
#SBATCH --job-name=s4_expr_rbh_ident
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
# Insumos del score: (a) DIAMOND reciproco OGS3<->RefSeq Tcas, (b) identidad Dmai vs Tcas, (c) salmon sobre 19 corridas.
set -euo pipefail
date +"inicio %F %T"
R=/scratch/PROYECTO/USUARIO/dmaidis
SX="singularity exec -B /scratch:/scratch"
DM="$SX $R/sif/diamond.sif diamond"
mkdir -p out/score && cd out/score

# (a) mejor hit reciproco OGS3 <-> RefSeq (todas las isoformas; se colapsa a gen despues)
zcat $R/data/tcas_map/OGS3_proteins.fasta.gz > ogs3.faa
$DM makedb --in ogs3.faa -d ogs3 --quiet
$DM makedb --in $R/data/tcas_map/protein.faa -d refseq --quiet
F="6 qseqid sseqid pident length qlen slen evalue bitscore"
$DM blastp -q ogs3.faa -d refseq -o ogs3_vs_refseq.tsv --id 90 --query-cover 70 --subject-cover 70 \
    --max-target-seqs 5 --threads 32 --quiet --outfmt $F &
$DM blastp -q $R/data/tcas_map/protein.faa -d ogs3 -o refseq_vs_ogs3.tsv --id 90 --query-cover 70 --subject-cover 70 \
    --max-target-seqs 5 --threads 32 --quiet --outfmt $F &
# (b) identidad Dmai vs Tcas (proteomas usados en OrthoFinder)
$DM makedb --in $R/orthofinder/proteomas/Tribolium_castaneum.faa -d tcas_of --quiet
$DM blastp -q $R/orthofinder/proteomas/Dalbulus_maidis.faa -d tcas_of -o dmai_vs_tcas.tsv --more-sensitive \
    --max-target-seqs 50 --evalue 1e-5 --threads 32 --quiet --outfmt $F &
wait
wc -l ogs3_vs_refseq.tsv refseq_vs_ogs3.tsv dmai_vs_tcas.tsv

# (c) expresion: CDS de Helixer -> salmon, 4 corridas en paralelo x 24 hilos
$SX $R/sif/gffread.sif gffread $R/out/helixer/Dmaidis_helixer.gff3 \
    -g $R/data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna -x dmai_cds.fa
$SX $R/sif/salmon.sif salmon index -t dmai_cds.fa -i salmon_idx -p 32 > salmon_index.log 2>&1
mkdir -p quant fq
cuant() {
  r=$1
  [ -s quant/$r/quant.sf ] && return 0
  $SX $R/sif/sra-tools.sif fasterq-dump --split-files -e 8 -t fq/tmp_$r -O fq $R/data/rnaseq/$r.sra > /dev/null 2>&1
  $SX $R/sif/salmon.sif salmon quant -i salmon_idx -l A -1 fq/${r}_1.fastq -2 fq/${r}_2.fastq -p 16 \
      --validateMappings -o quant/$r > quant_$r.log 2>&1
  rm -rf fq/${r}_*.fastq fq/tmp_$r
  echo "$r $(grep -oE '"percent_mapped": [0-9.]+' quant/$r/aux_info/meta_info.json)"
}
export -f cuant; export SX R
ls $R/data/rnaseq/*.sra | xargs -n1 basename | sed 's/.sra$//' | xargs -P4 -I{} bash -c 'cuant {}'
ls quant/*/quant.sf | wc -l
date +"fin %F %T"
