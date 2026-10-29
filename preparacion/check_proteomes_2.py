#!/usr/bin/env python3
"""
check_proteomes.py  (v2)

Audita que ensamblados hay en NCBI para armar el set de proteomas de un
analisis de ortologia (OrthoFinder) con Dalbulus maidis como blanco.

Cambios respecto de v1
----------------------
1. BUG ARREGLADO: 'NaN or 0.0' devuelve NaN en Python (NaN es truthy), asi
   que los ensamblados anotados SIN BUSCO publicado recibian puntaje NaN y
   quedaban por debajo de los no anotados. Por eso Drosophila y Laodelphax
   aparecian como "sin anotacion" teniendo versiones anotadas.
2. El bonus por 'reference genome' baja de 10 a 3: hacia ganar al ensamblado
   designado como referencia aunque fuera viejo y peor (caso Diaphorina).
3. TRIAGE: ahora lista TODOS los ensamblados anotados por especie, no solo
   el mejor, para poder elegir a mano.
4. Seccion aparte para ensamblados SIN anotar de la especie blanco, con sus
   metricas, para evaluar si conviene anotarlos.
5. Diagnostico del problema de duplicados: separa el D esperable por
   isoformas del D por duplicacion haplotipica.

Uso:
    pip install requests pandas
    python check_proteomes.py

Salidas:
    proteomas_ncbi_completo.csv
    proteomas_ncbi_anotados.csv
    descargar_proteomas.sh
    (y un triage completo por pantalla)
"""

from __future__ import annotations

import sys
import time
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# CONFIGURACION
# ---------------------------------------------------------------------------

API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
NCBI_API_KEY: Optional[str] = None      # opcional: sube el limite a 10 req/s
PAUSA = 0.4                              # segundos entre consultas

# La especie blanco se trata aparte: nos interesan sus ensamblados aunque
# NO esten anotados, para decidir si vale la pena anotarlos.
ESPECIE_BLANCO = "Dalbulus maidis"

ESPECIES = [
    # taxon, rol, grupo, anillo (0 = fuera del gradiente hemiptero)
    {"taxon": "Tribolium castaneum",        "rol": "fenotipo (iBeetle)",      "grupo": "Coleoptera",     "anillo": 0},
    {"taxon": "Drosophila melanogaster",    "rol": "anotacion / GO",          "grupo": "Diptera",        "anillo": 0},
    {"taxon": "Apis mellifera",             "rol": "calibrador no-blanco",    "grupo": "Hymenoptera",    "anillo": 0},

    {"taxon": "Macrosteles quadrilineatus", "rol": "anillo 1 Deltocephalinae","grupo": "Cicadellidae",   "anillo": 1},
    {"taxon": "Homalodisca vitripennis",    "rol": "anillo 1 Cicadellinae",   "grupo": "Cicadellidae",   "anillo": 1},
    {"taxon": "Empoasca onukii",            "rol": "anillo 1 Typhlocybinae",  "grupo": "Cicadellidae",   "anillo": 1},
    {"taxon": "Cicadella viridis",          "rol": "anillo 1 Cicadellinae",   "grupo": "Cicadellidae",   "anillo": 1},
    {"taxon": "Nephotettix cincticeps",     "rol": "anillo 1 Deltocephalinae","grupo": "Cicadellidae",   "anillo": 1},
    {"taxon": "Psammotettix striatus",      "rol": "anillo 1 Deltocephalinae","grupo": "Cicadellidae",   "anillo": 1},
    {"taxon": "Recilia dorsalis",           "rol": "anillo 1 Deltocephalinae","grupo": "Cicadellidae",   "anillo": 1},

    {"taxon": "Nilaparvata lugens",         "rol": "anillo 2 Delphacidae",    "grupo": "Fulgoromorpha",  "anillo": 2},
    {"taxon": "Laodelphax striatellus",     "rol": "anillo 2 Delphacidae",    "grupo": "Fulgoromorpha",  "anillo": 2},
    {"taxon": "Sogatella furcifera",        "rol": "anillo 2 Delphacidae",    "grupo": "Fulgoromorpha",  "anillo": 2},

    {"taxon": "Diaphorina citri",           "rol": "anillo 3 psilido",        "grupo": "Sternorrhyncha", "anillo": 3},
    {"taxon": "Acyrthosiphon pisum",        "rol": "anillo 3 afido",          "grupo": "Sternorrhyncha", "anillo": 3},
    {"taxon": "Bemisia tabaci",             "rol": "anillo 3 mosca blanca",   "grupo": "Sternorrhyncha", "anillo": 3},

    {"taxon": "Daphnia pulex",              "rol": "outgroup",                "grupo": "Crustacea",      "anillo": 0},

    {"taxon": ESPECIE_BLANCO,               "rol": "BLANCO",                  "grupo": "Cicadellidae",   "anillo": 1},
]

UMBRALES = {
    "busco_completo_min": 90.0,
    "busco_duplicado_max": 5.0,
    "busco_fragmentado_max": 5.0,
    "genes_min": 8000,
}

# Especies donde conviene NO usar NCBI
FUENTE_ALTERNATIVA = {
    "Drosophila melanogaster": "FlyBase (mejor anotacion y fuente de los GO)",
}


# ---------------------------------------------------------------------------
# UTILIDADES
# ---------------------------------------------------------------------------

def _num(x, default=0.0) -> float:
    """NaN/None -> default.  NO usar 'or': en Python NaN es truthy."""
    if x is None:
        return default
    try:
        if pd.isna(x):
            return default
    except (TypeError, ValueError):
        pass
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def _get(d: Any, *claves, default=None):
    for k in claves:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
        if d is None:
            return default
    return d


def _headers() -> Dict[str, str]:
    h = {"Accept": "application/json"}
    if NCBI_API_KEY:
        h["api-key"] = NCBI_API_KEY
    return h


# ---------------------------------------------------------------------------
# CONSULTA
# ---------------------------------------------------------------------------

def fetch_assemblies(taxon: str) -> List[Dict[str, Any]]:
    url = f"{API_BASE}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report"
    params: Dict[str, Any] = {"page_size": 100, "filters.assembly_version": "current"}
    reportes: List[Dict[str, Any]] = []
    token = None

    while True:
        if token:
            params["page_token"] = token
        try:
            r = requests.get(url, params=params, headers=_headers(), timeout=60)
        except requests.RequestException as e:
            print(f"  error de red: {e}", file=sys.stderr)
            return reportes

        if r.status_code == 404:
            return reportes
        if not r.ok:
            print(f"  HTTP {r.status_code}: {r.text[:200]}", file=sys.stderr)
            return reportes

        data = r.json()
        reportes.extend(data.get("reports") or [])
        token = data.get("next_page_token")
        if not token:
            break
        time.sleep(PAUSA)

    return reportes


def parse_report(rep: Dict[str, Any], meta: Dict[str, Any]) -> Dict[str, Any]:
    ann = rep.get("annotation_info") or {}
    stats = rep.get("assembly_stats") or {}
    info = rep.get("assembly_info") or {}
    busco = ann.get("busco") or {}

    def pct(x):
        return round(x * 100, 2) if isinstance(x, (int, float)) else None

    def intify(x):
        try:
            return int(x)
        except (TypeError, ValueError):
            return None

    total_len = intify(stats.get("total_sequence_length"))

    return {
        "taxon_consultado": meta["taxon"],
        "rol": meta["rol"],
        "grupo": meta["grupo"],
        "anillo": meta["anillo"],
        "organismo": _get(rep, "organism", "organism_name"),
        "accession": rep.get("accession"),
        "nombre_ensamblado": info.get("assembly_name"),
        "nivel": info.get("assembly_level"),
        "categoria_refseq": info.get("refseq_category"),
        "fecha": info.get("release_date"),
        "anotado": bool(ann),
        "proveedor_anotacion": ann.get("provider"),
        "release_anotacion": ann.get("release_date"),
        "genes_codificantes": _get(ann, "stats", "gene_counts", "protein_coding"),
        "genes_totales": _get(ann, "stats", "gene_counts", "total"),
        "tamano_mb": round(total_len / 1e6, 1) if total_len else None,
        "contig_n50": intify(stats.get("contig_n50")),
        "scaffold_n50": intify(stats.get("scaffold_n50")),
        "n_contigs": intify(stats.get("number_of_contigs")),
        "n_scaffolds": intify(stats.get("number_of_scaffolds")),
        "busco_lineage": busco.get("busco_lineage"),
        "busco_completo_pct": pct(busco.get("complete")),
        "busco_unico_pct": pct(busco.get("single_copy")),
        "busco_duplicado_pct": pct(busco.get("duplicated")),
        "busco_fragmentado_pct": pct(busco.get("fragmented")),
        "busco_ausente_pct": pct(busco.get("missing")),
    }


# ---------------------------------------------------------------------------
# EVALUACION
# ---------------------------------------------------------------------------

def evaluar(f: pd.Series) -> str:
    if not f["anotado"]:
        return "SIN ANOTACION"

    motivos: List[str] = []
    genes = f["genes_codificantes"]
    if genes is not None and not pd.isna(genes) and genes < UMBRALES["genes_min"]:
        motivos.append(f"pocos genes ({int(genes)})")

    c, d, fr = f["busco_completo_pct"], f["busco_duplicado_pct"], f["busco_fragmentado_pct"]

    if c is None or pd.isna(c):
        return "REVISAR: sin BUSCO publicado -> correrlo local"

    if c < UMBRALES["busco_completo_min"]:
        motivos.append(f"completo {c}%")
    if not pd.isna(_num(d, float("nan"))) and _num(d) > UMBRALES["busco_duplicado_max"]:
        motivos.append(f"DUPLICADOS {d}%")
    if not pd.isna(_num(fr, float("nan"))) and _num(fr) > UMBRALES["busco_fragmentado_max"]:
        motivos.append(f"fragmentado {fr}%")

    return "ACEPTA" if not motivos else "REVISAR: " + "; ".join(motivos)


def puntaje(f: pd.Series) -> float:
    """Ordena candidatos dentro de una especie. Clave: _num(), no 'or'."""
    if not f["anotado"]:
        return -100.0
    c = _num(f["busco_completo_pct"], 50.0)      # sin BUSCO -> neutro
    d = _num(f["busco_duplicado_pct"], 10.0)
    n50 = _num(f["contig_n50"], 0.0)
    bonus = 3.0 if f.get("categoria_refseq") else 0.0   # era 10: sesgaba a refs viejas
    return c - 2.0 * d + bonus + min(n50 / 1e7, 5.0)


def diagnostico_duplicados(f: pd.Series) -> str:
    """Hipotesis sobre el origen de un BUSCO duplicado alto."""
    d = _num(f["busco_duplicado_pct"], 0.0)
    if d <= UMBRALES["busco_duplicado_max"]:
        return ""
    nivel = (f.get("nivel") or "").lower()
    genes = _num(f["genes_codificantes"], 0.0)

    pistas = []
    if genes > 18000:
        pistas.append("muchos genes: posible exceso de isoformas o haplotipos")
    if "scaffold" in nivel or "contig" in nivel:
        pistas.append("nivel scaffold/contig: compatible con haplotipos sin colapsar")
    if "chromosome" in nivel:
        pistas.append("nivel cromosoma: mas probable que sea isoformas")

    base = f"D={d}%. " + ("; ".join(pistas) if pistas else "")
    return base + " -> TEST: filtrar isoformas y correr BUSCO -m proteins"


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    filas: List[Dict[str, Any]] = []
    print(f"Consultando NCBI para {len(ESPECIES)} especies...\n")

    for meta in ESPECIES:
        print(f"  {meta['taxon']:<34}", end="", flush=True)
        reps = fetch_assemblies(meta["taxon"])
        if not reps:
            print("sin ensamblados")
        else:
            n_anot = sum(1 for r in reps if r.get("annotation_info"))
            print(f"{len(reps):>4} ensamblados ({n_anot} anotados)")
            for r in reps:
                filas.append(parse_report(r, meta))
        time.sleep(PAUSA)

    if not filas:
        sys.exit("\nNo se obtuvo nada. Revisa conexion y API_BASE.")

    df = pd.DataFrame(filas)
    df["veredicto"] = df.apply(evaluar, axis=1)
    df["_puntaje"] = df.apply(puntaje, axis=1)
    df["diagnostico"] = df.apply(diagnostico_duplicados, axis=1)
    df.to_csv("proteomas_ncbi_completo.csv", index=False)

    anot = df[df["anotado"]].copy()
    anot = anot.sort_values(["anillo", "taxon_consultado", "_puntaje"],
                            ascending=[True, True, False])
    anot.to_csv("proteomas_ncbi_anotados.csv", index=False)

    # ------------------------------------------------------------------
    # TRIAGE: todos los anotados, por especie
    # ------------------------------------------------------------------
    cols = ["accession", "nivel", "categoria_refseq", "fecha",
            "genes_codificantes", "tamano_mb", "contig_n50",
            "busco_completo_pct", "busco_duplicado_pct",
            "busco_fragmentado_pct", "veredicto"]

    print("\n" + "=" * 120)
    print("TRIAGE - TODOS LOS ENSAMBLADOS ANOTADOS, POR ESPECIE")
    print("=" * 120)

    for anillo in sorted(anot["anillo"].unique()):
        sub = anot[anot["anillo"] == anillo]
        etiqueta = {0: "REFERENCIA / OUTGROUP", 1: "ANILLO 1 - Cicadellidae",
                    2: "ANILLO 2 - Fulgoromorpha",
                    3: "ANILLO 3 - Sternorrhyncha"}.get(anillo, f"anillo {anillo}")
        print(f"\n\n### {etiqueta} ###")
        for sp, g in sub.groupby("taxon_consultado", sort=False):
            alt = FUENTE_ALTERNATIVA.get(sp)
            nota = f"   [usar en su lugar: {alt}]" if alt else ""
            print(f"\n-- {sp}{nota}")
            with pd.option_context("display.width", 220, "display.max_colwidth", 38):
                print(g[cols].to_string(index=False))

    # ------------------------------------------------------------------
    # ESPECIE BLANCO: ensamblados sin anotar
    # ------------------------------------------------------------------
    blanco = df[df["taxon_consultado"] == ESPECIE_BLANCO]
    print("\n\n" + "=" * 120)
    print(f"ENSAMBLADOS DE {ESPECIE_BLANCO.upper()} (anotados o no)")
    print("=" * 120)
    if blanco.empty:
        print("  ninguno")
    else:
        cb = ["accession", "nivel", "fecha", "anotado", "tamano_mb",
              "contig_n50", "scaffold_n50", "n_contigs", "n_scaffolds",
              "busco_completo_pct"]
        with pd.option_context("display.width", 220):
            print(blanco.sort_values("fecha", ascending=False)[cb].to_string(index=False))
        mejor_b = blanco.sort_values(
            ["anotado", "scaffold_n50"], ascending=[False, False]).iloc[0]
        if not mejor_b["anotado"]:
            print("\n  Ninguno esta anotado. Si el nivel es 'Chromosome' y el N50 es alto,")
            print("  conviene ANOTARLO con BRAKER3/GALBA usando el transcriptoma como")
            print("  evidencia, en lugar de traducir el transcriptoma con TransDecoder:")
            print("  da estructura genica real y un conteo de paralogos confiable.")

    # ------------------------------------------------------------------
    # DIAGNOSTICO DE DUPLICADOS
    # ------------------------------------------------------------------
    problemas = anot[anot["diagnostico"] != ""].groupby(
        "taxon_consultado", as_index=False).head(1)
    if not problemas.empty:
        print("\n\n" + "=" * 120)
        print("DUPLICADOS ALTOS - DIAGNOSTICO")
        print("=" * 120)
        for _, r in problemas.iterrows():
            print(f"\n  {r['taxon_consultado']}  ({r['accession']}, {r['nivel']})")
            print(f"    {r['diagnostico']}")
        print("\n  Dos causas posibles, con consecuencias opuestas:")
        print("    A) isoformas multiples   -> SE ARREGLA con una_isoforma_por_gen.py")
        print("    B) haplotipos sin colapsar -> NO se arregla; evaluar purge_dups")
        print("  El test decisivo es empirico:")
        print("    python una_isoforma_por_gen.py --faa X.faa --gff X.gff \\")
        print("           --out proteomas/Especie.faa --etiqueta Esp")
        print("    busco -i proteomas/Especie.faa -l insecta_odb10 -m proteins -c 8 \\")
        print("          -o busco_Especie_filtrado")

    # ------------------------------------------------------------------
    # SET PROPUESTO + SCRIPT DE DESCARGA
    # ------------------------------------------------------------------
    elegidos = (anot[~anot["taxon_consultado"].isin(FUENTE_ALTERNATIVA)]
                .groupby("taxon_consultado", as_index=False).head(1)
                .sort_values(["anillo", "taxon_consultado"]))

    print("\n\n" + "=" * 120)
    print("SET PROPUESTO (mejor anotado por especie)")
    print("=" * 120)
    cs = ["taxon_consultado", "anillo", "rol", "accession", "nivel",
          "genes_codificantes", "busco_completo_pct", "busco_duplicado_pct", "veredicto"]
    with pd.option_context("display.width", 220, "display.max_colwidth", 34):
        print(elegidos[cs].to_string(index=False))

    n_anillo1 = int((elegidos["anillo"] == 1).sum())
    if n_anillo1 < 2:
        print("\n  [!] Menos de 2 cicadelidos en el anillo 1. Con uno solo no se puede")
        print("      distinguir una duplicacion propia de D. maidis de un artefacto")
        print("      de ensamblado: el filtro 1:1 queda sin respaldo.")

    with open("descargar_proteomas.sh", "w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("# Generado por check_proteomes.py v2\n")
        fh.write("# Requiere: conda install -c conda-forge ncbi-datasets-cli\n")
        fh.write("set -euo pipefail\nmkdir -p proteomas_crudos proteomas\n\n")
        for _, r in elegidos.iterrows():
            slug = str(r["taxon_consultado"]).replace(" ", "_")
            etq = "".join(w[0] for w in str(r["taxon_consultado"]).split()[:1]) + \
                  str(r["taxon_consultado"]).split()[-1][:3]
            fh.write(f"# {r['taxon_consultado']} [{r['rol']}] "
                     f"BUSCO C={r['busco_completo_pct']} D={r['busco_duplicado_pct']}\n")
            fh.write(f"datasets download genome accession {r['accession']} \\\n"
                     f"    --include protein,gff3 \\\n"
                     f"    --filename proteomas_crudos/{slug}.zip\n")
            fh.write(f"unzip -o proteomas_crudos/{slug}.zip -d proteomas_crudos/{slug}\n")
            fh.write(f"python una_isoforma_por_gen.py \\\n"
                     f"    --faa proteomas_crudos/{slug}/ncbi_dataset/data/*/protein.faa \\\n"
                     f"    --gff proteomas_crudos/{slug}/ncbi_dataset/data/*/genomic.gff \\\n"
                     f"    --out proteomas/{slug}.faa --etiqueta {etq}\n\n")
        for sp, fuente in FUENTE_ALTERNATIVA.items():
            fh.write(f"# {sp}: NO desde NCBI -> {fuente}\n")
            fh.write(f"#   descargar manualmente y correr una_isoforma_por_gen.py\n\n")

    print("\n\nArchivos escritos:")
    print("  proteomas_ncbi_completo.csv")
    print("  proteomas_ncbi_anotados.csv")
    print("  descargar_proteomas.sh")

    sin_nada = [m["taxon"] for m in ESPECIES
                if m["taxon"] not in set(df["taxon_consultado"])]
    sin_anot = sorted(set(df[~df["anotado"]]["taxon_consultado"]) -
                      set(anot["taxon_consultado"]))
    if sin_nada or sin_anot:
        print("\nNO UTILIZABLES PARA ORTHOFINDER:")
        for t in sin_nada:
            print(f"  - {t}  (sin ensamblados)")
        for t in sin_anot:
            marca = "  <- es la especie blanco, se anota aparte" if t == ESPECIE_BLANCO else ""
            print(f"  - {t}  (ensamblados sin anotar){marca}")


if __name__ == "__main__":
    main()
