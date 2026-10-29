#!/usr/bin/env python3
"""
check_proteomes.py

Busca en NCBI que ensamblados hay disponibles para una lista de especies,
si estan anotados, y con que calidad (BUSCO, N50, nivel de ensamblado).

Pensado para armar el set de proteomas de un analisis de ortologia
(OrthoFinder) con Dalbulus maidis como especie blanco.

Uso:
    pip install requests pandas
    python check_proteomes.py

Salidas:
    proteomas_ncbi_completo.csv   - todos los ensamblados encontrados
    proteomas_ncbi_mejores.csv    - el mejor candidato por especie
    descargar_proteomas.sh        - comandos de descarga para los elegidos

Nota: la API de NCBI cambia de version cada tanto. Si algo falla con un
error 404, revisa API_BASE abajo contra la documentacion vigente en
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
"""

import json
import sys
import time
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# CONFIGURACION
# ---------------------------------------------------------------------------

API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

# Especies a evaluar. El campo "rol" documenta por que esta cada una:
# cada especie del set tiene que tener un trabajo asignado.
ESPECIES = [
    # --- Donantes de fenotipo y anotacion ---
    {"taxon": "Tribolium castaneum",      "rol": "fenotipo (iBeetle)",        "grupo": "Coleoptera"},
    {"taxon": "Drosophila melanogaster",  "rol": "anotacion (FlyBase/GO)",    "grupo": "Diptera"},
    {"taxon": "Apis mellifera",           "rol": "calibrador no-blanco",      "grupo": "Hymenoptera"},

    # --- Anillo 1: Cicadellidae (misma familia que D. maidis) ---
    {"taxon": "Macrosteles quadrilineatus", "rol": "anillo 1 - Deltocephalinae", "grupo": "Cicadellidae"},
    {"taxon": "Homalodisca vitripennis",    "rol": "anillo 1 - Cicadellinae",    "grupo": "Cicadellidae"},
    {"taxon": "Graminella nigrifrons",      "rol": "anillo 1 - verificar",       "grupo": "Cicadellidae"},
    {"taxon": "Nephotettix cincticeps",     "rol": "anillo 1 - alternativa",     "grupo": "Cicadellidae"},

    # --- Anillo 2: Auchenorrhyncha no-Cicadellidae ---
    {"taxon": "Nilaparvata lugens",       "rol": "anillo 2 - Delphacidae",    "grupo": "Fulgoromorpha"},
    {"taxon": "Laodelphax striatellus",   "rol": "anillo 2 - alternativa",    "grupo": "Fulgoromorpha"},
    {"taxon": "Peregrinus maidis",        "rol": "anillo 2 - verificar",      "grupo": "Fulgoromorpha"},

    # --- Anillo 3: Sternorrhyncha ---
    {"taxon": "Diaphorina citri",         "rol": "anillo 3 - psilido",        "grupo": "Sternorrhyncha"},
    {"taxon": "Acyrthosiphon pisum",      "rol": "anillo 3 - afido",          "grupo": "Sternorrhyncha"},
    {"taxon": "Bemisia tabaci",           "rol": "anillo 3 - mosca blanca",   "grupo": "Sternorrhyncha"},

    # --- Outgroup (raiz del arbol) ---
    {"taxon": "Daphnia pulex",            "rol": "outgroup",                  "grupo": "Crustacea"},
    {"taxon": "Tetranychus urticae",      "rol": "outgroup alternativo",      "grupo": "Chelicerata"},

    # --- La especie blanco (probablemente sin proteoma; se verifica igual) ---
    {"taxon": "Dalbulus maidis",          "rol": "BLANCO",                    "grupo": "Cicadellidae"},
]

# Umbrales de aceptacion para un proteoma de referencia.
# El de duplicados es el que mas importa: un D alto en una especie de
# referencia mete redundancia en TODOS los ortogrupos, no solo en los suyos.
UMBRALES = {
    "busco_completo_min": 90.0,   # %
    "busco_duplicado_max": 5.0,   # %
    "busco_fragmentado_max": 5.0, # %
    "genes_min": 8000,            # genes codificantes; muy por debajo = anotacion pobre
}

PAUSA_ENTRE_CONSULTAS = 0.4  # segundos; NCBI limita ~3 req/s sin API key

# Si tenes una API key de NCBI (opcional, sube el limite a 10 req/s):
NCBI_API_KEY: Optional[str] = None


# ---------------------------------------------------------------------------
# CONSULTA A LA API
# ---------------------------------------------------------------------------

def _headers() -> Dict[str, str]:
    h = {"Accept": "application/json"}
    if NCBI_API_KEY:
        h["api-key"] = NCBI_API_KEY
    return h


def fetch_assemblies(taxon: str, solo_anotados: bool = False) -> List[Dict[str, Any]]:
    """Devuelve todos los reportes de ensamblado de NCBI para un taxon.

    No filtra por 'reference_only' a proposito: para muchas especies no-modelo
    el ensamblado util no esta marcado como referencia, y filtrarlo de entrada
    hace que parezca que no hay nada disponible.
    """
    url = f"{API_BASE}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report"
    params: Dict[str, Any] = {
        "page_size": 100,
        "filters.assembly_version": "current",
    }
    if solo_anotados:
        params["filters.has_annotation"] = "true"

    reportes: List[Dict[str, Any]] = []
    token = None

    while True:
        if token:
            params["page_token"] = token
        try:
            r = requests.get(url, params=params, headers=_headers(), timeout=60)
        except requests.RequestException as e:
            print(f"  [red] error de red para '{taxon}': {e}", file=sys.stderr)
            return reportes

        if r.status_code == 404:
            # Taxon no encontrado o sin ensamblados
            return reportes
        if not r.ok:
            print(f"  [!] HTTP {r.status_code} para '{taxon}': {r.text[:200]}",
                  file=sys.stderr)
            return reportes

        data = r.json()
        reportes.extend(data.get("reports", []) or [])
        token = data.get("next_page_token")
        if not token:
            break
        time.sleep(PAUSA_ENTRE_CONSULTAS)

    return reportes


# ---------------------------------------------------------------------------
# PARSEO
# ---------------------------------------------------------------------------

def _get(d: Any, *claves, default=None):
    """Acceso anidado tolerante a claves faltantes."""
    for k in claves:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
        if d is None:
            return default
    return d


def parse_report(rep: Dict[str, Any], meta: Dict[str, str]) -> Dict[str, Any]:
    """Aplana un reporte de NCBI a una fila de tabla."""
    ann = rep.get("annotation_info") or {}
    stats = rep.get("assembly_stats") or {}
    info = rep.get("assembly_info") or {}

    # BUSCO: NCBI lo calcula y lo publica dentro de annotation_info para
    # los ensamblados anotados. Los valores vienen como fraccion (0-1).
    busco = ann.get("busco") or {}

    def pct(x):
        return round(x * 100, 2) if isinstance(x, (int, float)) else None

    total_len = stats.get("total_sequence_length")
    try:
        total_len = int(total_len)
    except (TypeError, ValueError):
        total_len = None

    return {
        "taxon_consultado": meta["taxon"],
        "rol": meta["rol"],
        "grupo": meta["grupo"],
        "organismo": _get(rep, "organism", "organism_name"),
        "tax_id": _get(rep, "organism", "tax_id"),
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
        "contig_n50": stats.get("contig_n50"),
        "scaffold_n50": stats.get("scaffold_n50"),
        "n_contigs": stats.get("number_of_contigs"),
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

def evaluar(fila: pd.Series) -> str:
    """Devuelve ACEPTA / REVISAR / RECHAZA y el motivo."""
    motivos = []

    if not fila["anotado"]:
        return "RECHAZA: sin anotacion (no hay proteoma)"

    genes = fila["genes_codificantes"]
    if pd.notna(genes) and genes < UMBRALES["genes_min"]:
        motivos.append(f"pocos genes ({int(genes)})")

    c = fila["busco_completo_pct"]
    d = fila["busco_duplicado_pct"]
    f = fila["busco_fragmentado_pct"]

    if pd.isna(c):
        motivos.append("sin BUSCO publicado -> correrlo localmente")
        return "REVISAR: " + "; ".join(motivos)

    if c < UMBRALES["busco_completo_min"]:
        motivos.append(f"completo bajo ({c}%)")
    if pd.notna(d) and d > UMBRALES["busco_duplicado_max"]:
        motivos.append(f"DUPLICADOS alto ({d}%)")
    if pd.notna(f) and f > UMBRALES["busco_fragmentado_max"]:
        motivos.append(f"fragmentado alto ({f}%)")

    if not motivos:
        return "ACEPTA"
    return "REVISAR: " + "; ".join(motivos)


def puntaje(fila: pd.Series) -> float:
    """Puntaje para ordenar candidatos dentro de una misma especie.
    Prioriza: anotado > BUSCO completo alto > duplicados bajo > N50 alto.
    """
    if not fila["anotado"]:
        return -1.0
    c = fila["busco_completo_pct"] or 0.0
    d = fila["busco_duplicado_pct"] or 0.0
    n50 = fila["contig_n50"] or 0
    bonus_ref = 10.0 if fila["categoria_refseq"] else 0.0
    return c - 2.0 * d + bonus_ref + min(n50 / 1e7, 5.0)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    filas: List[Dict[str, Any]] = []

    print(f"Consultando NCBI para {len(ESPECIES)} especies...\n")

    for meta in ESPECIES:
        taxon = meta["taxon"]
        print(f"  {taxon:<32}", end="", flush=True)
        reps = fetch_assemblies(taxon)
        if not reps:
            print("sin ensamblados")
            filas.append({
                "taxon_consultado": taxon, "rol": meta["rol"], "grupo": meta["grupo"],
                "organismo": None, "accession": None, "anotado": False,
            })
        else:
            anotados = sum(1 for r in reps if r.get("annotation_info"))
            print(f"{len(reps):>3} ensamblados ({anotados} anotados)")
            for r in reps:
                filas.append(parse_report(r, meta))
        time.sleep(PAUSA_ENTRE_CONSULTAS)

    df = pd.DataFrame(filas)
    if df.empty:
        print("\nNo se obtuvo nada. Revisa conexion y API_BASE.")
        return

    # Rellenar columnas que puedan faltar si ninguna especie las trajo
    for col in ["busco_completo_pct", "busco_duplicado_pct", "busco_fragmentado_pct",
                "genes_codificantes", "contig_n50", "categoria_refseq"]:
        if col not in df.columns:
            df[col] = pd.NA

    df["veredicto"] = df.apply(evaluar, axis=1)
    df["_puntaje"] = df.apply(puntaje, axis=1)

    df = df.sort_values(["taxon_consultado", "_puntaje"], ascending=[True, False])
    df.to_csv("proteomas_ncbi_completo.csv", index=False)

    # Mejor candidato por especie
    mejores = df.groupby("taxon_consultado", as_index=False).head(1).copy()
    mejores = mejores.sort_values("grupo")
    mejores.to_csv("proteomas_ncbi_mejores.csv", index=False)

    # ---- Resumen en pantalla ----
    cols = ["taxon_consultado", "grupo", "accession", "nivel", "anotado",
            "genes_codificantes", "busco_completo_pct", "busco_duplicado_pct",
            "veredicto"]
    print("\n" + "=" * 110)
    print("MEJOR CANDIDATO POR ESPECIE")
    print("=" * 110)
    with pd.option_context("display.max_colwidth", 42, "display.width", 200):
        print(mejores[cols].to_string(index=False))

    # ---- Script de descarga para los que pasan ----
    ok = mejores[mejores["veredicto"].str.startswith(("ACEPTA", "REVISAR"))
                 & mejores["anotado"]]

    with open("descargar_proteomas.sh", "w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("# Generado por check_proteomes.py\n")
        fh.write("# Requiere: conda install -c conda-forge ncbi-datasets-cli\n")
        fh.write("set -euo pipefail\nmkdir -p proteomas_crudos\n\n")
        for _, r in ok.iterrows():
            slug = str(r["taxon_consultado"]).replace(" ", "_")
            fh.write(f"# {r['taxon_consultado']}  [{r['rol']}]  "
                     f"BUSCO C={r['busco_completo_pct']} D={r['busco_duplicado_pct']}\n")
            fh.write(f"datasets download genome accession {r['accession']} \\\n")
            fh.write(f"    --include protein,gff3 \\\n")
            fh.write(f"    --filename proteomas_crudos/{slug}.zip\n")
            fh.write(f"unzip -o proteomas_crudos/{slug}.zip "
                     f"-d proteomas_crudos/{slug}\n\n")

    print("\nArchivos escritos:")
    print("  proteomas_ncbi_completo.csv")
    print("  proteomas_ncbi_mejores.csv")
    print("  descargar_proteomas.sh")

    # ---- Avisos ----
    faltan = mejores[~mejores["anotado"]]["taxon_consultado"].tolist()
    if faltan:
        print("\nSIN PROTEOMA ANOTADO (no sirven para OrthoFinder):")
        for t in faltan:
            print(f"  - {t}")
        print("  -> buscar reemplazo del mismo anillo filogenetico")

    sin_busco = mejores[mejores["anotado"] & mejores["busco_completo_pct"].isna()]
    if not sin_busco.empty:
        print("\nANOTADOS PERO SIN BUSCO PUBLICADO (correrlo local):")
        for t in sin_busco["taxon_consultado"]:
            print(f"  - {t}")
        print("  busco -i <proteoma>.faa -l insecta_odb10 -m proteins -c 8 -o busco_<sp>")


if __name__ == "__main__":
    main()
