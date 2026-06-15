"""
revel.py
========
REVEL score (and CADD as a bonus) via MyVariant.info, which serves dbNSFP
fields. No API key required.

REVEL is an ensemble missense pathogenicity predictor (0-1, higher = more
likely pathogenic). It is only defined for missense SNVs.

API: https://myvariant.info/v1/variant/<hgvs>?assembly=hg38&fields=dbnsfp.revel,dbnsfp.cadd

The variant id uses HGVS genomic notation, e.g. chr3:g.177026399C>T.
"""

MYVARIANT_API = "https://myvariant.info/v1/variant"


def _hgvs_g(chrom, pos, ref, alt):
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return f"{chrom}:g.{pos}{ref}>{alt}"


def _first(value):
    """dbNSFP fields are sometimes scalars, sometimes lists (one per transcript)."""
    if isinstance(value, list):
        nums = [v for v in value if isinstance(v, (int, float))]
        return max(nums) if nums else (value[0] if value else None)
    return value


def query_revel(client, chrom, pos, ref, alt, assembly="hg38"):
    """
    Fetch REVEL (and CADD phred) for a single SNV from MyVariant.info / dbNSFP.

    Returns a dict with the REVEL score and a CADD phred score when available.
    """
    result = {
        "found": False,
        "revel": None,
        "cadd_phred": None,
        "url": "",
        "error": "",
    }

    hgvs = _hgvs_g(chrom, pos, ref, alt)
    params = {
        "assembly": assembly,
        "fields": "dbnsfp.revel,dbnsfp.cadd,cadd.phred",
    }
    # URL-encode happens in requests; cache key keeps it readable.
    url = f"{MYVARIANT_API}/{hgvs}"
    cache_key = f"myvariant:{assembly}:{hgvs}"
    result["url"] = f"https://myvariant.info/v1/variant/{hgvs}?assembly={assembly}"

    data = client.get_json(url, params=params, cache_key=cache_key)

    if "_error" in data:
        # 404 just means the variant isn't in dbNSFP (e.g. non-missense)
        if data.get("_status") == 404:
            return result
        result["error"] = data["_error"]
        return result

    dbnsfp = data.get("dbnsfp") or {}

    # REVEL: dbnsfp.revel.score (newer) or dbnsfp.revel (older flat)
    revel = dbnsfp.get("revel")
    if isinstance(revel, dict):
        result["revel"] = _first(revel.get("score"))
    else:
        result["revel"] = _first(revel)

    # CADD phred: prefer top-level cadd.phred, fall back to dbnsfp.cadd.phred
    cadd = data.get("cadd") or {}
    cadd_phred = cadd.get("phred")
    if cadd_phred is None:
        dbnsfp_cadd = dbnsfp.get("cadd") or {}
        cadd_phred = dbnsfp_cadd.get("phred")
    result["cadd_phred"] = _first(cadd_phred)

    if result["revel"] is not None or result["cadd_phred"] is not None:
        result["found"] = True

    return result


def interpret(revel):
    """Common REVEL interpretation bands (no official cutoff; these are typical)."""
    if revel is None:
        return ""
    if revel >= 0.7:
        return "likely pathogenic"
    if revel >= 0.5:
        return "possibly pathogenic"
    if revel >= 0.25:
        return "uncertain"
    return "likely benign"
