"""
gtex.py
=======
GTEx bulk-tissue median gene expression (GTEx v8) via the GTEx Portal API.

For each variant's gene we report the median TPM of the top-3 expressing
tissues, and always also include skeletal muscle, cultured fibroblasts, and
whole blood (appended if they are not already in the top 3).

API: https://gtexportal.org/api/v2
  1. /reference/gene?geneId=<symbol>            -> versioned gencodeId
  2. /expression/medianGeneExpression?gencodeId=<id>&datasetId=gtex_v8
                                                 -> median TPM per tissue
"""

GTEX_API = "https://gtexportal.org/api/v2"
DATASET = "gtex_v8"

# Tissues to always surface, with friendly display names.
ALWAYS_INCLUDE = [
    ("Muscle_Skeletal", "Muscle (skeletal)"),
    ("Cells_Cultured_fibroblasts", "Fibroblasts (cultured)"),
    ("Whole_Blood", "Whole blood"),
]

TOP_N = 3


def _pretty(tissue_id):
    return tissue_id.replace("_", " ")


def _resolve_gencode_id(client, gene_symbol, gene_id):
    """Look up the GTEx versioned gencodeId for a gene."""
    for query in (gene_symbol, gene_id):
        if not query:
            continue
        data = client.get_json(f"{GTEX_API}/reference/gene",
                               params={"geneId": query},
                               cache_key=f"gtex_gene:{query}")
        rows = data.get("data") if isinstance(data, dict) else None
        if rows:
            # Prefer an exact (case-insensitive) symbol match if present.
            if gene_symbol:
                for r in rows:
                    if (r.get("geneSymbol", "").upper()
                            == gene_symbol.upper()):
                        return r.get("gencodeId")
            return rows[0].get("gencodeId")
    return None


def query_gtex(client, gene_symbol, gene_id=None, enabled=True):
    """
    Return median-TPM expression for the top tissues plus muscle/fibroblast/
    blood for `gene_symbol`.
    """
    result = {"found": False, "tissues": [], "max_median": 0.0,
              "error": "", "url": ""}
    if not enabled:
        result["error"] = "skipped"
        return result
    if not gene_symbol and not gene_id:
        result["error"] = "no gene"
        return result

    result["url"] = f"https://gtexportal.org/home/gene/{gene_symbol or gene_id}"

    gencode = _resolve_gencode_id(client, gene_symbol, gene_id)
    if not gencode:
        result["error"] = "gene not found in GTEx"
        return result

    data = client.get_json(f"{GTEX_API}/expression/medianGeneExpression",
                           params={"gencodeId": gencode, "datasetId": DATASET},
                           cache_key=f"gtex_med:{gencode}")
    if not isinstance(data, dict) or "_error" in data:
        result["error"] = (data.get("_error", "GTEx query failed")
                           if isinstance(data, dict) else "GTEx query failed")
        return result

    rows = data.get("data") or []
    if not rows:
        result["error"] = "no GTEx expression data"
        return result

    by_id = {}
    for r in rows:
        tid = r.get("tissueSiteDetailId")
        med = r.get("median")
        if tid is None or med is None:
            continue
        by_id[tid] = float(med)

    if not by_id:
        result["error"] = "no GTEx expression data"
        return result

    ranked = sorted(by_id.items(), key=lambda kv: kv[1], reverse=True)
    rank_of = {tid: i + 1 for i, (tid, _) in enumerate(ranked)}

    tissues = []
    chosen = set()
    for tid, med in ranked[:TOP_N]:
        tissues.append({"name": _pretty(tid), "median": med,
                        "top": True, "rank": rank_of[tid]})
        chosen.add(tid)

    # Append the always-include tissues that didn't make the top N.
    for tid, label in ALWAYS_INCLUDE:
        if tid in chosen:
            # Update the display name to the friendly label.
            for t in tissues:
                if t["name"] == _pretty(tid):
                    t["name"] = label
            continue
        med = by_id.get(tid)
        tissues.append({
            "name": label,
            "median": med if med is not None else 0.0,
            "top": False,
            "rank": rank_of.get(tid),
        })

    result["found"] = True
    result["tissues"] = tissues
    result["max_median"] = max(t["median"] for t in tissues) or 1.0
    return result
