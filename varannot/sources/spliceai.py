"""
spliceai.py
===========
SpliceAI splice-disruption scores via the Broad Institute SpliceAI-lookup API.

Returns the four SpliceAI delta scores and the overall maximum (delta score):
  - DS_AG  acceptor gain
  - DS_AL  acceptor loss
  - DS_DG  donor gain
  - DS_DL  donor loss
plus the corresponding base-pair offsets (DP_*).

API: https://spliceai-38-xwkwwwxdwq-uc.a.run.app  (GRCh38 Cloud Run service)

IMPORTANT: This is the public *interactive* endpoint, rate limited to only a
few requests per user per minute. For batches it is queried with an extra
delay (see SPLICEAI_MIN_INTERVAL). For heavy/automated use, run your own
instance (see github.com/broadinstitute/SpliceAI-lookup) and point
--spliceai-url at it.
"""

import time

# GRCh38 SpliceAI Cloud Run service used by spliceailookup.broadinstitute.org
SPLICEAI_URL_38 = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app"
SPLICEAI_LOOKUP_WEB = "https://spliceailookup.broadinstitute.org"

# Public interactive endpoint is heavily rate limited; be polite.
SPLICEAI_MIN_INTERVAL = 6.0  # seconds between SpliceAI calls

# The Cloud Run service can take well over the default 30s to cold-start and
# run the model, so SpliceAI requests get a generous timeout of their own.
SPLICEAI_TIMEOUT = 180

_last_call = {"t": 0.0}


def _throttle():
    elapsed = time.time() - _last_call["t"]
    if elapsed < SPLICEAI_MIN_INTERVAL:
        time.sleep(SPLICEAI_MIN_INTERVAL - elapsed)
    _last_call["t"] = time.time()


def _max_delta(scores):
    """Largest of the four delta scores (the standard SpliceAI 'DS' max)."""
    vals = [
        scores.get("DS_AG"), scores.get("DS_AL"),
        scores.get("DS_DG"), scores.get("DS_DL"),
    ]
    vals = [v for v in vals if isinstance(v, (int, float))]
    return max(vals) if vals else None


def query_spliceai(client, chrom, pos, ref, alt, base_url=None, enabled=True,
                   distance=50):
    """
    Query SpliceAI for a single SNV/indel (GRCh38).

    Returns a dict with per-transcript scores and the overall max delta.
    If `enabled` is False, returns a skipped marker without any network call.
    """
    result = {
        "enabled": enabled,
        "found": False,
        "max_delta": None,
        "transcripts": [],   # list of {symbol, DS_AG, DS_AL, DS_DG, DS_DL, DP_*}
        "web_url": "",
        "error": "",
    }
    if not enabled:
        result["error"] = "skipped (use --spliceai to enable)"
        return result

    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    variant = f"{chrom}-{pos}-{ref}-{alt}"
    # The Broad service expects the request on the "/spliceai/" path, not the
    # bare host. Accept either form for a custom --spliceai-url.
    base = (base_url or SPLICEAI_URL_38).rstrip("/")
    if not base.endswith("/spliceai"):
        base = base + "/spliceai/"
    else:
        base = base + "/"
    result["web_url"] = f"{SPLICEAI_LOOKUP_WEB}/#variant={variant}&hg=38"

    params = {
        "hg": "38",
        "variant": variant,
        "distance": str(distance),
        "mask": "0",
    }
    cache_key = f"spliceai38:{variant}:d{distance}"

    # Only throttle when we will actually hit the network (cache miss).
    cache_path = client._cache_path(cache_key)
    import os
    if not os.path.exists(cache_path):
        _throttle()

    # This is a Cloud Run service that can cold-start slowly (model loading),
    # so give it a much longer timeout than the default API calls.
    data = client.get_json(base, params=params, cache_key=cache_key,
                           timeout=SPLICEAI_TIMEOUT)

    if "_error" in data:
        result["error"] = data["_error"]
        return result

    # The API returns {"scores": [...]} where each item is a transcript-level
    # score string or dict. Newer versions return structured fields.
    scores_list = data.get("scores") or data.get("variants") or []
    if isinstance(scores_list, dict):
        scores_list = [scores_list]

    parsed = []
    for item in scores_list:
        # item may be a raw string like
        # "TRANSCRIPT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"
        if isinstance(item, str):
            fields = item.split("|")
            if len(fields) >= 6:
                try:
                    entry = {
                        "symbol": fields[1] if len(fields) > 1 else "",
                        "DS_AG": float(fields[2]),
                        "DS_AL": float(fields[3]),
                        "DS_DG": float(fields[4]),
                        "DS_DL": float(fields[5]),
                    }
                    if len(fields) >= 10:
                        entry.update({
                            "DP_AG": _safe_int(fields[6]),
                            "DP_AL": _safe_int(fields[7]),
                            "DP_DG": _safe_int(fields[8]),
                            "DP_DL": _safe_int(fields[9]),
                        })
                    parsed.append(entry)
                except (ValueError, IndexError):
                    continue
        elif isinstance(item, dict):
            entry = {
                "symbol": (item.get("g_name") or item.get("SYMBOL")
                           or item.get("symbol") or item.get("t_name", "")),
                "DS_AG": _safe_float(item.get("DS_AG")),
                "DS_AL": _safe_float(item.get("DS_AL")),
                "DS_DG": _safe_float(item.get("DS_DG")),
                "DS_DL": _safe_float(item.get("DS_DL")),
                "DP_AG": _safe_int(item.get("DP_AG")),
                "DP_AL": _safe_int(item.get("DP_AL")),
                "DP_DG": _safe_int(item.get("DP_DG")),
                "DP_DL": _safe_int(item.get("DP_DL")),
            }
            parsed.append(entry)

    if not parsed:
        # Some responses wrap the message under "error" or "message"
        msg = data.get("error") or data.get("message") or "no SpliceAI scores returned"
        result["error"] = str(msg)
        return result

    result["found"] = True
    result["transcripts"] = parsed
    # Overall max across all transcripts
    overall = [_max_delta(t) for t in parsed]
    overall = [v for v in overall if v is not None]
    result["max_delta"] = max(overall) if overall else None
    return result


def _safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _safe_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def interpret(max_delta):
    """Standard SpliceAI interpretation thresholds."""
    if max_delta is None:
        return ""
    if max_delta >= 0.8:
        return "high precision"
    if max_delta >= 0.5:
        return "likely splice-altering"
    if max_delta >= 0.2:
        return "possible / high recall"
    return "low"
