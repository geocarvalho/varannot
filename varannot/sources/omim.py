"""
omim.py
=======
OMIM lookup using the locally downloaded OMIM data files (no API key needed).

OMIM grants two separate things: an API key, and the bulk download files.
If you have the download files (mim2gene.txt, mimTitles.txt, genemap2.txt,
morbidmap.txt), this module reads `genemap2.txt` directly — fully offline,
no key, no rate limits.

genemap2.txt provides, per gene:
  - the gene MIM number and approved symbol
  - the associated phenotypes, each with its phenotype MIM number,
    mapping key, and inheritance

The Phenotypes field is parsed following OMIM's own recommended logic.

(An API-based fallback, query_omim_api, is kept for users who instead have a
licensed API key.)
"""

import os
import re

OMIM_ENTRY_URL = "https://www.omim.org/entry"

# Mapping-key meanings (the "(n)" after each phenotype MIM number in genemap2)
MAPPING_KEY_DESC = {
    "1": "mapped by wildtype gene",
    "2": "disorder mapped",
    "3": "molecular basis known",
    "4": "contiguous gene del/dup syndrome",
}


# ---------------------------------------------------------------------------
# Loading + indexing genemap2.txt
# ---------------------------------------------------------------------------
def load_genemap2(path):
    """
    Parse genemap2.txt into an index keyed by gene symbol (uppercased).

    Returns a dict: { "SYMBOL": record, ... } where record has the gene MIM,
    approved symbol, and a list of parsed phenotypes. Both the approved symbol
    and each alias in the "Gene Symbols" column are indexed, so lookups are
    robust to symbol differences between VEP and OMIM.
    """
    if not path or not os.path.exists(path):
        return None

    index = {}
    header_cols = None

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                # The last comment line before data is the column header.
                # It begins with "# Chromosome" — capture it.
                stripped = line.lstrip("#").strip()
                if stripped.lower().startswith("chromosome"):
                    header_cols = [c.strip() for c in stripped.split("\t")]
                continue

            if header_cols is None:
                # No header seen; fall back to known column order.
                header_cols = _DEFAULT_COLUMNS

            fields = line.split("\t")
            row = dict(zip(header_cols, fields))

            mim = _col(row, "MIM Number")
            approved = _col(row, "Approved Symbol") or _col(row, "Approved Gene Symbol")
            gene_symbols = _col(row, "Gene Symbols") or _col(row, "Gene/Locus And Other Related Symbols")
            gene_name = _col(row, "Gene Name")
            phenotypes_raw = _col(row, "Phenotypes")

            if not mim:
                continue

            record = {
                "gene_mim": str(mim),
                "approved_symbol": approved,
                "gene_name": gene_name,
                "phenotypes": _parse_phenotypes(phenotypes_raw),
            }

            # Index by approved symbol and every alias
            keys = set()
            if approved:
                keys.add(approved.upper())
            for alias in (gene_symbols or "").split(","):
                alias = alias.strip()
                if alias:
                    keys.add(alias.upper())
            for k in keys:
                # Prefer a record that actually has phenotypes if there is a clash
                if k not in index or (record["phenotypes"] and not index[k]["phenotypes"]):
                    index[k] = record

    return index


# Column order used only if the header comment is missing (defensive default).
_DEFAULT_COLUMNS = [
    "Chromosome", "Genomic Position Start", "Genomic Position End",
    "Cyto Location", "Computed Cyto Location", "MIM Number",
    "Gene Symbols", "Gene Name", "Approved Symbol",
    "Entrez Gene ID", "Ensembl Gene ID", "Comments", "Phenotypes",
    "Mouse Gene Symbol/ID",
]


def _col(row, name):
    """Case/space-tolerant column getter."""
    if name in row:
        return (row[name] or "").strip()
    # try a loose match
    target = name.lower().replace(" ", "")
    for k, v in row.items():
        if k.lower().replace(" ", "") == target:
            return (v or "").strip()
    return ""


def _parse_phenotypes(phenotypes_raw):
    """
    Parse the genemap2 'Phenotypes' field into structured entries.

    Format (phenotypes separated by ';'):
        Phenotype name, 600000 (3), Autosomal recessive
    Some phenotypes lack a MIM number:
        Phenotype name (2), Autosomal dominant

    Follows OMIM's own recommended regex approach.
    """
    out = []
    if not phenotypes_raw:
        return out

    for chunk in phenotypes_raw.split(";"):
        pheno = chunk.strip()
        if not pheno:
            continue

        name = pheno
        mim = ""
        mapping_key = ""
        inheritance = ""

        # With a phenotype MIM number: "name, 600000 (3), Inheritance"
        m = re.match(r"^(.*),\s*(\d{6})\s*\((\d)\)(?:,\s*(.*))?$", pheno)
        if m:
            name = m.group(1).strip()
            mim = m.group(2)
            mapping_key = m.group(3)
            inheritance = (m.group(4) or "").strip()
        else:
            # Without a MIM number: "name (3), Inheritance"
            m2 = re.match(r"^(.*)\((\d)\)(?:,\s*(.*))?$", pheno)
            if m2:
                name = m2.group(1).strip().rstrip(",").strip()
                mapping_key = m2.group(2)
                inheritance = (m2.group(3) or "").strip()

        # Strip a leading brace/bracket/question mark marker OMIM uses for
        # non-disease / provisional phenotypes, but keep it readable.
        out.append({
            "phenotype": name,
            "mim": mim,
            "inheritance": inheritance,
            "mapping_key": mapping_key,
        })

    return out


# ---------------------------------------------------------------------------
# Query (local files)
# ---------------------------------------------------------------------------
def query_omim_local(index, gene_symbol):
    """Look up a gene in the loaded genemap2 index."""
    result = {"found": False, "gene_mim": "", "gene_url": "", "phenotypes": []}
    if index is None:
        result["error"] = "genemap2.txt not loaded"
        return result
    if not gene_symbol:
        return result

    record = index.get(gene_symbol.upper())
    if not record:
        return result

    result["found"] = True
    result["gene_mim"] = record["gene_mim"]
    result["gene_url"] = f"{OMIM_ENTRY_URL}/{record['gene_mim']}"
    result["phenotypes"] = record["phenotypes"]
    return result


# ---------------------------------------------------------------------------
# Query (API fallback — only if the user has a licensed key instead of files)
# ---------------------------------------------------------------------------
OMIM_API = "https://api.omim.org/api"


def query_omim_api(client, gene_symbol, api_key):
    """Look up OMIM via the official API (requires a licensed key)."""
    result = {"found": False, "gene_mim": "", "gene_url": "", "phenotypes": []}
    if not api_key:
        result["error"] = "no OMIM API key provided"
        return result
    if not gene_symbol:
        return result

    params = {
        "search": f"+{gene_symbol}",
        "filter": "prefix:gene",
        "include": "geneMap",
        "format": "json",
        "apiKey": api_key,
    }
    data = client.get_json(f"{OMIM_API}/entry/search", params=params,
                           cache_key=f"omim:{gene_symbol}")
    if "_error" in data:
        result["error"] = data["_error"]
        return result

    search_response = (data.get("omim") or {}).get("searchResponse", {})
    entries = search_response.get("entryList", [])
    if not entries:
        return result

    entry = entries[0].get("entry", {})
    mim = entry.get("mimNumber", "")
    result["found"] = True
    result["gene_mim"] = str(mim)
    result["gene_url"] = f"{OMIM_ENTRY_URL}/{mim}"

    gene_map = entry.get("geneMap", {})
    for pheno_wrap in gene_map.get("phenotypeMapList", []):
        pheno = pheno_wrap.get("phenotypeMap", {})
        result["phenotypes"].append({
            "phenotype": pheno.get("phenotype", ""),
            "mim": str(pheno.get("phenotypeMimNumber", "")),
            "inheritance": pheno.get("phenotypeInheritance", "") or "",
            "mapping_key": pheno.get("phenotypeMappingKey", ""),
        })
    return result
