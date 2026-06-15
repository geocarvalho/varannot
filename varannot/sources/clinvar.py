"""
clinvar.py
==========
ClinVar lookup via NCBI E-utilities (esearch + esummary).

Provides:
  - Whether the variant is in ClinVar
  - Clinical significance / classification
  - Number of submissions/entries
  - Review status (star rating)
  - A ClinVar link

Endpoint: https://eutils.ncbi.nlm.nih.gov/entrez/eutils
"""

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation"


def query_clinvar(client, chrom, pos, ref, alt, api_key=None):
    """
    Search ClinVar for a variant by genomic position.

    Uses an esearch on the position/alleles then esummary on the top hit.
    """
    chrom = str(chrom).replace("chr", "")
    result = {
        "found": False,
        "clinical_significance": "",
        "review_status": "",
        "n_submissions": 0,
        "variation_id": "",
        "url": "",
        "title": "",
    }

    # Search term: by genomic location and base change.
    # ClinVar indexes variants; we search the position and filter by alleles in title.
    term = f"{chrom}[chr] AND {pos}[chrpos] AND single_gene[prop]"
    params = {
        "db": "clinvar",
        "term": f"{chrom}[chr] AND {pos}[chrpos37] OR {pos}[chrpos38]",
        "retmode": "json",
        "retmax": 20,
    }
    if api_key:
        params["api_key"] = api_key

    search = client.get_json(f"{EUTILS}/esearch.fcgi", params=params,
                             cache_key=f"clinvar_search:{chrom}:{pos}:{ref}:{alt}")
    if "_error" in search:
        result["error"] = search["_error"]
        return result

    idlist = (search.get("esearchresult") or {}).get("idlist", [])
    if not idlist:
        return result

    # Pull summaries for candidate IDs and match the allele change
    sum_params = {
        "db": "clinvar",
        "id": ",".join(idlist),
        "retmode": "json",
    }
    if api_key:
        sum_params["api_key"] = api_key

    summary = client.get_json(f"{EUTILS}/esummary.fcgi", params=sum_params,
                              cache_key=f"clinvar_summary:{chrom}:{pos}:{ref}:{alt}")
    if "_error" in summary:
        result["error"] = summary["_error"]
        return result

    docs = (summary.get("result") or {})
    uids = docs.get("uids", [])

    best = None
    for uid in uids:
        rec = docs.get(uid, {})
        title = rec.get("title", "")
        # Prefer a record whose title mentions the exact substitution
        if f"{ref}>{alt}" in title or f"{ref}/{alt}" in title:
            best = rec
            break
    if best is None and uids:
        best = docs.get(uids[0], {})

    if not best:
        return result

    result["found"] = True
    result["variation_id"] = best.get("uid", "")
    result["title"] = best.get("title", "")
    result["url"] = f"{CLINVAR_URL}/{best.get('uid', '')}"

    germline = best.get("germline_classification", {}) or {}
    result["clinical_significance"] = germline.get("description", "") or \
        (best.get("clinical_significance", {}) or {}).get("description", "")
    result["review_status"] = germline.get("review_status", "")

    # Count submissions if available
    trait_set = germline.get("trait_set", [])
    result["n_submissions"] = best.get("supporting_submissions", {}).get("scv", []).__len__() \
        if isinstance(best.get("supporting_submissions", {}), dict) else 0

    return result
