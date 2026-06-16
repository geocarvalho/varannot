"""
autoacmg.py
===========
ACMG/AMP variant classification via a **self-hosted** AutoACMG instance
(bihealth/auto-acmg).

AutoACMG does **not** offer a public hosted server — you must run your own:
    https://auto-acmg.readthedocs.io/   (Docker; needs SeqRepo + a REEV proxy)

Once it's running, it exposes:
    GET {base}/api/v1/predict/seqvar?variant_name=chr1:228282272:G:A&genome_release=GRCh38

The response gives a per-criterion ACMG prediction (PVS1, PS1 ... BP7), each
with a ``prediction`` (applicable / not_applicable / not_automated / ...) and a
``strength``. AutoACMG itself does not emit a single Pathogenic/Benign verdict,
so we collect the **applicable** criteria and combine them into an overall
classification using the ACMG/AMP combining rules (Richards et al., 2015).

This source is opt-in: pass ``enabled=True`` and a ``base_url`` (or rely on the
default ``http://localhost:8080``). When the instance is unreachable the report
just shows an informative note instead of failing.
"""

DEFAULT_URL = "http://localhost:8080"

# AutoACMG predictions can be slow (they fan out to several backends per
# variant), so allow a generous timeout.
REQUEST_TIMEOUT = 180

# Map the strength enum AutoACMG returns onto the ACMG evidence buckets used by
# the combining rules. Using the *returned* strength (not the criterion name)
# respects AutoACMG's strength adjustments (e.g. a downgraded PVS1).
_STRENGTH_BUCKET = {
    "pathogenic_very_strong": "PVS",
    "pathogenic_strong": "PS",
    "pathogenic_moderate": "PM",
    "pathogenic_supporting": "PP",
    "benign_stand_alone": "BA",
    "benign_strong": "BS",
    "benign_supporting": "BP",
}

# Human-readable strength labels for the report.
_BUCKET_LABEL = {
    "PVS": "very strong",
    "PS": "strong",
    "PM": "moderate",
    "PP": "supporting",
    "BA": "stand-alone",
    "BS": "strong",
    "BP": "supporting",
}

# Sort order so pathogenic-strong criteria come first in the displayed list.
_BUCKET_ORDER = {"PVS": 0, "PS": 1, "PM": 2, "PP": 3, "BA": 4, "BS": 5, "BP": 6}


def _variant_name(chrom, pos, ref, alt):
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return f"{chrom}:{pos}:{ref}:{alt}"


def classify(criteria):
    """Combine applicable ACMG criteria into an overall classification.

    ``criteria`` is the list of applicable-criteria dicts (each with a
    ``bucket`` of PVS/PS/PM/PP/BA/BS/BP). Returns ``(label, abbrev)`` following
    the ACMG/AMP combining rules (Richards et al., 2015).
    """
    counts = {"PVS": 0, "PS": 0, "PM": 0, "PP": 0, "BA": 0, "BS": 0, "BP": 0}
    for c in criteria:
        b = c.get("bucket")
        if b in counts:
            counts[b] += 1
    pvs, ps, pm, pp = counts["PVS"], counts["PS"], counts["PM"], counts["PP"]
    ba, bs, bp = counts["BA"], counts["BS"], counts["BP"]

    pathogenic = (
        (pvs >= 1 and (ps >= 1 or pm >= 2 or (pm == 1 and pp >= 1) or pp >= 2))
        or ps >= 2
        or (ps == 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)))
    )
    likely_pathogenic = (
        (pvs >= 1 and pm == 1)
        or (ps == 1 and 1 <= pm <= 2)
        or (ps == 1 and pp >= 2)
        or pm >= 3
        or (pm == 2 and pp >= 2)
        or (pm == 1 and pp >= 4)
    )
    benign = ba >= 1 or bs >= 2
    likely_benign = (bs == 1 and bp >= 1) or bp >= 2

    path_side = pathogenic or likely_pathogenic
    benign_side = benign or likely_benign

    # Both the pathogenic and benign rules fire -> conflicting -> VUS.
    if path_side and benign_side:
        return "Uncertain significance (conflicting)", "VUS"
    if pathogenic:
        return "Pathogenic", "P"
    if likely_pathogenic:
        return "Likely pathogenic", "LP"
    if benign:
        return "Benign", "B"
    if likely_benign:
        return "Likely benign", "LB"
    return "Uncertain significance", "VUS"


def query_autoacmg(client, chrom, pos, ref, alt, base_url=None, enabled=False,
                   genome_release="GRCh38"):
    """Fetch and classify ACMG criteria for one variant from AutoACMG.

    Returns a dict with the combined ``classification``/``abbrev`` plus the list
    of applicable ``criteria``. Safe to call when disabled or when the instance
    is unreachable — the result simply carries ``found=False`` and a note.
    """
    result = {
        "enabled": enabled,
        "found": False,
        "error": "",
        "url": "",
        "criteria": [],
        "classification": "",
        "abbrev": "",
        "gene_symbol": "",
        "phgvs": "",
        "transcript_id": "",
    }
    if not enabled:
        return result

    base = (base_url or DEFAULT_URL).rstrip("/")
    name = _variant_name(chrom, pos, ref, alt)
    url = f"{base}/api/v1/predict/seqvar"
    params = {"variant_name": name, "genome_release": genome_release}
    result["url"] = f"{url}?variant_name={name}&genome_release={genome_release}"
    cache_key = f"autoacmg:{base}:{genome_release}:{name}"

    data = client.get_json(url, params=params, cache_key=cache_key,
                           timeout=REQUEST_TIMEOUT)

    if "_error" in data:
        result["error"] = data["_error"]
        return result

    pred = (data or {}).get("prediction") or {}
    pdata = pred.get("data") or {}
    result["gene_symbol"] = pdata.get("gene_symbol", "") or ""
    result["phgvs"] = pdata.get("pHGVS", "") or ""
    result["transcript_id"] = pdata.get("transcript_id", "") or ""

    criteria = pred.get("criteria") or {}
    applicable = []
    for key, crit in criteria.items():
        if not isinstance(crit, dict):
            continue
        if crit.get("prediction") != "applicable":
            continue
        strength = crit.get("strength", "") or ""
        bucket = _STRENGTH_BUCKET.get(strength, "")
        applicable.append({
            "name": crit.get("name") or key.upper(),
            "strength": strength,
            "bucket": bucket,
            "strength_label": _BUCKET_LABEL.get(bucket, ""),
            "summary": (crit.get("summary") or "").strip(),
        })

    applicable.sort(key=lambda c: _BUCKET_ORDER.get(c["bucket"], 99))
    result["criteria"] = applicable
    result["found"] = True

    label, abbrev = classify(applicable)
    result["classification"] = label
    result["abbrev"] = abbrev
    return result
