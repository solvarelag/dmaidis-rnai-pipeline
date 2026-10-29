#!/usr/bin/env python3
"""
una_isoforma_por_gen.py

Reduce un proteoma de NCBI a una isoforma por gen (la mas larga).

Por que hace falta: un proteoma con todas las isoformas infla
artificialmente los ortogrupos. OrthoFinder va a reportar relaciones
1:muchos que no son parasitos biologicos sino variantes de splicing, y el
filtro "1:1" deja de significar algo.

Usa el GFF3 de NCBI para mapear proteina -> gen (via Dbxref=GeneID:).

Uso:
    pip install biopython
    python una_isoforma_por_gen.py \
        --faa proteomas_crudos/Tribolium_castaneum/ncbi_dataset/data/GCF_*/protein.faa \
        --gff proteomas_crudos/Tribolium_castaneum/ncbi_dataset/data/GCF_*/genomic.gff \
        --out proteomas/Tribolium_castaneum.faa \
        --etiqueta Tcas

La etiqueta se antepone a cada ID (Tcas|XP_12345) para que despues sea
trivial saber de que especie viene cada secuencia en los arboles de genes.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path


def abrir(ruta: Path):
    """Abre texto plano o .gz indistintamente."""
    if str(ruta).endswith(".gz"):
        return gzip.open(ruta, "rt")
    return open(ruta, "r")


def mapear_proteina_a_gen(gff: Path) -> dict:
    """Lee un GFF3 de NCBI y devuelve {protein_id: gene_id}.

    En los GFF de NCBI, las features CDS llevan en la columna 9:
        protein_id=XP_008191234.1
        Dbxref=GeneID:103312345,Genbank:XP_008191234.1
    Se usa GeneID porque es estable; si falta, se cae a 'gene=' y por
    ultimo al 'Parent='.
    """
    prot2gen = {}
    re_prot = re.compile(r"protein_id=([^;]+)")
    re_geneid = re.compile(r"GeneID:(\d+)")
    re_gene = re.compile(r"(?:^|;)gene=([^;]+)")
    re_parent = re.compile(r"(?:^|;)Parent=([^;]+)")

    with abrir(gff) as fh:
        for linea in fh:
            if linea.startswith("#"):
                continue
            campos = linea.rstrip("\n").split("\t")
            if len(campos) < 9 or campos[2] != "CDS":
                continue
            attrs = campos[8]

            m = re_prot.search(attrs)
            if not m:
                continue
            pid = m.group(1)
            if pid in prot2gen:
                continue  # ya mapeado por otro exon del mismo CDS

            g = re_geneid.search(attrs)
            if g:
                prot2gen[pid] = "GeneID:" + g.group(1)
                continue
            g = re_gene.search(attrs)
            if g:
                prot2gen[pid] = g.group(1)
                continue
            g = re_parent.search(attrs)
            if g:
                prot2gen[pid] = g.group(1)

    return prot2gen


def leer_fasta(faa: Path):
    """Generador de (id, descripcion, secuencia). Evita dependencias."""
    sid, desc, partes = None, "", []
    with abrir(faa) as fh:
        for linea in fh:
            linea = linea.rstrip("\n")
            if linea.startswith(">"):
                if sid is not None:
                    yield sid, desc, "".join(partes)
                cabecera = linea[1:]
                sid = cabecera.split()[0] if cabecera else ""
                desc = cabecera[len(sid):].strip()
                partes = []
            else:
                partes.append(linea.strip())
    if sid is not None:
        yield sid, desc, "".join(partes)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", required=True, type=Path, help="protein.faa de NCBI")
    ap.add_argument("--gff", required=True, type=Path, help="genomic.gff de NCBI")
    ap.add_argument("--out", required=True, type=Path, help="FASTA de salida")
    ap.add_argument("--etiqueta", default=None,
                    help="prefijo de especie para los IDs, ej. Tcas")
    ap.add_argument("--min-aa", type=int, default=30,
                    help="descartar proteinas mas cortas que esto (default 30)")
    args = ap.parse_args()

    if not args.faa.exists():
        sys.exit(f"No existe: {args.faa}")
    if not args.gff.exists():
        sys.exit(f"No existe: {args.gff}")

    print(f"Mapeando proteina -> gen desde {args.gff.name} ...")
    prot2gen = mapear_proteina_a_gen(args.gff)
    print(f"  {len(prot2gen):,} proteinas mapeadas a gen")

    # Elegir la isoforma mas larga por gen
    mejor = {}            # gen -> (largo, sid, secuencia)
    sin_mapeo = 0
    total = 0
    cortas = 0

    for sid, _desc, seq in leer_fasta(args.faa):
        total += 1
        seq = seq.replace("*", "")
        if len(seq) < args.min_aa:
            cortas += 1
            continue
        gen = prot2gen.get(sid)
        if gen is None:
            # Sin mapeo: se conserva como gen propio para no perder datos,
            # pero se cuenta aparte porque es señal de GFF/FAA desparejos.
            gen = f"__sinmapeo__{sid}"
            sin_mapeo += 1
        largo = len(seq)
        if gen not in mejor or largo > mejor[gen][0]:
            mejor[gen] = (largo, sid, seq)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        for gen, (_largo, sid, seq) in sorted(mejor.items()):
            nuevo_id = f"{args.etiqueta}|{sid}" if args.etiqueta else sid
            fh.write(f">{nuevo_id} gene={gen}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # ---- Reporte ----
    print(f"\n  proteinas leidas      : {total:,}")
    print(f"  descartadas por cortas: {cortas:,}")
    print(f"  sin mapeo a gen       : {sin_mapeo:,}")
    print(f"  genes de salida       : {len(mejor):,}")
    if total:
        reduccion = 100 * (1 - len(mejor) / total)
        print(f"  reduccion             : {reduccion:.1f}%")
    print(f"\n  escrito: {args.out}")

    if sin_mapeo > 0.05 * total:
        print("\n  [!] Mas del 5% sin mapeo: verifica que el GFF y el FAA")
        print("      vengan del MISMO accession de ensamblado.")

    print("\n  Siguiente paso - verificar con BUSCO que D sea bajo:")
    print(f"      busco -i {args.out} -l insecta_odb10 -m proteins -c 8 \\")
    print(f"            -o busco_{args.out.stem}")


if __name__ == "__main__":
    main()
