#!/bin/bash
# T4: biblioteca de repeticiones de novo + enmascarado suave. Se corre desde lote2.
set -euo pipefail
date +"inicio %F %T"
SX="singularity exec -B /scratch:/scratch"
N=$SLURM_CPUS_PER_TASK
R=/scratch/PROYECTO/USUARIO/dmaidis
mkdir -p $R/out/repeats && cd $R/out/repeats
ln -sf $R/data/Dmaidis_genome/ncbi_dataset/data/GCA_975739785.1/GCA_975739785.1_ihDalMaid1.1_genomic.fna genome.fna
[ -s dmaidis.nsq ] || $SX $R/sif/tetools.sif BuildDatabase -name dmaidis genome.fna
[ -s dmaidis-families.fa ] || $SX $R/sif/tetools.sif RepeatModeler -database dmaidis -threads $N -LTRStruct
date +"repeatmodeler fin %F %T"
# -pa cuenta procesos de ~4 hilos cada uno
$SX $R/sif/tetools.sif RepeatMasker -pa $((N/4)) -lib dmaidis-families.fa -xsmall -gff -dir rmask genome.fna
grep -E "total length|bases masked" rmask/genome.fna.tbl
date +"fin %F %T"
