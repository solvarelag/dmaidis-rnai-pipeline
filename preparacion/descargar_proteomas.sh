#!/usr/bin/env bash
# Generado por check_proteomes.py v2
# Requiere: conda install -c conda-forge ncbi-datasets-cli
set -euo pipefail
mkdir -p proteomas_crudos proteomas

# Apis mellifera [calibrador no-blanco] BUSCO C=98.78 D=0.15
datasets download genome accession GCF_003254395.2 \
    --include protein,gff3 \
    --filename proteomas_crudos/Apis_mellifera.zip
unzip -o proteomas_crudos/Apis_mellifera.zip -d proteomas_crudos/Apis_mellifera
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Apis_mellifera/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Apis_mellifera/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Apis_mellifera.faa --etiqueta Amel

# Daphnia pulex [outgroup] BUSCO C=97.24 D=0.89
datasets download genome accession GCF_021134715.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Daphnia_pulex.zip
unzip -o proteomas_crudos/Daphnia_pulex.zip -d proteomas_crudos/Daphnia_pulex
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Daphnia_pulex/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Daphnia_pulex/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Daphnia_pulex.faa --etiqueta Dpul

# Tribolium castaneum [fenotipo (iBeetle)] BUSCO C=99.15 D=0.52
datasets download genome accession GCF_031307605.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Tribolium_castaneum.zip
unzip -o proteomas_crudos/Tribolium_castaneum.zip -d proteomas_crudos/Tribolium_castaneum
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Tribolium_castaneum/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Tribolium_castaneum/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Tribolium_castaneum.faa --etiqueta Tcas

# Homalodisca vitripennis [anillo 1 Cicadellinae] BUSCO C=92.75 D=10.8
datasets download genome accession GCF_021130785.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Homalodisca_vitripennis.zip
unzip -o proteomas_crudos/Homalodisca_vitripennis.zip -d proteomas_crudos/Homalodisca_vitripennis
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Homalodisca_vitripennis/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Homalodisca_vitripennis/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Homalodisca_vitripennis.faa --etiqueta Hvit

# Macrosteles quadrilineatus [anillo 1 Deltocephalinae] BUSCO C=97.49 D=38.09
datasets download genome accession GCF_028750875.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Macrosteles_quadrilineatus.zip
unzip -o proteomas_crudos/Macrosteles_quadrilineatus.zip -d proteomas_crudos/Macrosteles_quadrilineatus
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Macrosteles_quadrilineatus/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Macrosteles_quadrilineatus/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Macrosteles_quadrilineatus.faa --etiqueta Mqua

# Laodelphax striatellus [anillo 2 Delphacidae] BUSCO C=nan D=nan
datasets download genome accession GCA_003335185.2 \
    --include protein,gff3 \
    --filename proteomas_crudos/Laodelphax_striatellus.zip
unzip -o proteomas_crudos/Laodelphax_striatellus.zip -d proteomas_crudos/Laodelphax_striatellus
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Laodelphax_striatellus/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Laodelphax_striatellus/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Laodelphax_striatellus.faa --etiqueta Lstr

# Nilaparvata lugens [anillo 2 Delphacidae] BUSCO C=97.61 D=9.64
datasets download genome accession GCF_000757685.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Nilaparvata_lugens.zip
unzip -o proteomas_crudos/Nilaparvata_lugens.zip -d proteomas_crudos/Nilaparvata_lugens
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Nilaparvata_lugens/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Nilaparvata_lugens/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Nilaparvata_lugens.faa --etiqueta Nlug

# Acyrthosiphon pisum [anillo 3 afido] BUSCO C=98.65 D=2.63
datasets download genome accession GCF_005508785.2 \
    --include protein,gff3 \
    --filename proteomas_crudos/Acyrthosiphon_pisum.zip
unzip -o proteomas_crudos/Acyrthosiphon_pisum.zip -d proteomas_crudos/Acyrthosiphon_pisum
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Acyrthosiphon_pisum/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Acyrthosiphon_pisum/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Acyrthosiphon_pisum.faa --etiqueta Apis

# Bemisia tabaci [anillo 3 mosca blanca] BUSCO C=98.76 D=1.79
datasets download genome accession GCF_001854935.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Bemisia_tabaci.zip
unzip -o proteomas_crudos/Bemisia_tabaci.zip -d proteomas_crudos/Bemisia_tabaci
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Bemisia_tabaci/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Bemisia_tabaci/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Bemisia_tabaci.faa --etiqueta Btab

# Diaphorina citri [anillo 3 psilido] BUSCO C=86.77 D=4.86
datasets download genome accession GCF_000475195.1 \
    --include protein,gff3 \
    --filename proteomas_crudos/Diaphorina_citri.zip
unzip -o proteomas_crudos/Diaphorina_citri.zip -d proteomas_crudos/Diaphorina_citri
python una_isoforma_por_gen.py \
    --faa proteomas_crudos/Diaphorina_citri/ncbi_dataset/data/*/protein.faa \
    --gff proteomas_crudos/Diaphorina_citri/ncbi_dataset/data/*/genomic.gff \
    --out proteomas/Diaphorina_citri.faa --etiqueta Dcit

# Drosophila melanogaster: NO desde NCBI -> FlyBase (mejor anotacion y fuente de los GO)
#   descargar manualmente y correr una_isoforma_por_gen.py

