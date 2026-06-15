"""
gene_model.py
=============
Build a small "where is the variant in the gene" diagram for each variant.

It fetches the transcript's exon/intron structure from the Ensembl REST API
(GRCh38) and renders an inline SVG gene model: exons drawn to genomic scale
(coding regions thick, UTRs thin), introns as a connecting line, and the
variant marked with a pointer. For intronic / near-exon variants it also draws
a dashed connector to the nearest exon boundary and labels the distance.

The SVG is inline (no image files, no extra dependencies) and uses the report's
CSS custom properties so it matches the rest of the theme.
"""

ENSEMBL_SERVER = "https://rest.ensembl.org"


# ---------------------------------------------------------------------------
# Data fetch
# ---------------------------------------------------------------------------
def fetch_transcript_model(client, ensembl_tx_id):
    """Look up exon coordinates + CDS bounds for an Ensembl transcript."""
    if not ensembl_tx_id:
        return {"found": False, "error": "no Ensembl transcript id"}

    tx = ensembl_tx_id.split(".")[0]
    url = f"{ENSEMBL_SERVER}/lookup/id/{tx}"
    data = client.get_json(url, params={"expand": "1"},
                           cache_key=f"ensembl_lookup:{tx}")
    if not isinstance(data, dict) or "_error" in data:
        return {"found": False,
                "error": data.get("_error", "lookup failed")
                if isinstance(data, dict) else "lookup failed"}

    exons = data.get("Exon") or []
    if not exons:
        return {"found": False, "error": "no exon structure returned"}

    transl = data.get("Translation") or {}
    return {
        "found": True,
        "error": "",
        "strand": data.get("strand", 1),
        "exons": sorted(
            ({"start": int(e["start"]), "end": int(e["end"])} for e in exons),
            key=lambda e: e["start"],
        ),
        "cds_start": transl.get("start"),
        "cds_end": transl.get("end"),
        "n_exons": len(exons),
    }


# ---------------------------------------------------------------------------
# Locate the variant relative to the exon structure
# ---------------------------------------------------------------------------
def locate_variant(model, pos):
    """
    Work out whether `pos` sits in an exon and, if not, the distance to the
    nearest exon boundary. Exons are numbered in transcript (5'->3') order.
    """
    pos = int(pos)
    exons = model["exons"]
    strand = model.get("strand", 1) or 1
    # transcript order: 5'->3'
    ordered = exons if strand >= 0 else list(reversed(exons))

    for i, e in enumerate(ordered):
        if e["start"] <= pos <= e["end"]:
            return {
                "in_exon": True,
                "exon_index": i + 1,
                "n_exons": len(ordered),
                "distance": 0,
                "nearest_boundary": pos,
            }

    best = None
    for i, e in enumerate(ordered):
        for boundary in (e["start"], e["end"]):
            d = abs(pos - boundary)
            if best is None or d < best["distance"]:
                best = {
                    "in_exon": False,
                    "exon_index": i + 1,
                    "n_exons": len(ordered),
                    "distance": d,
                    "nearest_boundary": boundary,
                }
    return best or {"in_exon": False, "exon_index": 0, "n_exons": 0,
                    "distance": None, "nearest_boundary": pos}


def describe_location(loc):
    """Human-readable one-liner for the diagram caption."""
    if not loc or loc.get("n_exons") == 0:
        return "transcript structure unavailable"
    if loc["in_exon"]:
        return f"In exon {loc['exon_index']} of {loc['n_exons']}"
    d = loc["distance"]
    dist = f"{d:,} bp" if d is not None else "unknown distance"
    return f"Intronic — {dist} from nearest exon (exon {loc['exon_index']})"


# ---------------------------------------------------------------------------
# SVG rendering
# ---------------------------------------------------------------------------
def render_svg(model, pos, loc, width=560):
    """Render the gene model as an inline SVG string."""
    pos = int(pos)
    exons = model["exons"]
    strand = model.get("strand", 1) or 1
    cds_start = model.get("cds_start")
    cds_end = model.get("cds_end")

    ml, mr = 16, 16
    plot_w = width - ml - mr
    height = 78
    center_y = 46
    utr_h, cds_h = 9, 17

    g_min = min(e["start"] for e in exons)
    g_max = max(e["end"] for e in exons)
    lo = min(g_min, pos)
    hi = max(g_max, pos)
    span = (hi - lo) or 1

    def gx(g):
        frac = (g - lo) / span
        if strand < 0:
            frac = 1 - frac
        return ml + frac * plot_w

    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="100%" '
        f'preserveAspectRatio="xMidYMid meet" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'font-family="var(--sans, sans-serif)" role="img" '
        f'aria-label="gene model showing variant position">'
    ]

    # intron line across the exonic span
    x_a, x_b = gx(g_min), gx(g_max)
    parts.append(
        f'<line x1="{min(x_a, x_b):.1f}" y1="{center_y}" '
        f'x2="{max(x_a, x_b):.1f}" y2="{center_y}" '
        f'stroke="var(--gap, #b9b6ad)" stroke-width="1.5"/>'
    )

    var_exon = loc.get("exon_index") if loc and loc.get("in_exon") else None
    # exons numbered in transcript order; map back to genomic-sorted index
    ordered = exons if strand >= 0 else list(reversed(exons))

    def rect(x1, x2, h, fill, stroke=None):
        left = min(x1, x2)
        w = max(abs(x2 - x1), 2.0)
        y = center_y - h / 2
        s = (f'<rect x="{left:.1f}" y="{y:.1f}" width="{w:.1f}" '
             f'height="{h}" rx="1.5" fill="{fill}"')
        if stroke:
            s += f' stroke="{stroke}" stroke-width="1.5"'
        return s + '/>'

    for idx, e in enumerate(ordered, start=1):
        x1, x2 = gx(e["start"]), gx(e["end"])
        highlight = "var(--warn, #b5532a)" if idx == var_exon else None
        # UTR-height base box for the whole exon
        parts.append(rect(x1, x2, utr_h, "var(--accent-soft, #e6f0ef)",
                          stroke=highlight))
        # coding portion (taller, darker)
        if cds_start is not None and cds_end is not None:
            cs = max(e["start"], cds_start)
            ce = min(e["end"], cds_end)
            if cs <= ce:
                parts.append(rect(gx(cs), gx(ce), cds_h,
                                  "var(--accent, #2b6f6a)", stroke=highlight))

    # dashed connector to nearest exon boundary for intronic variants
    if loc and not loc.get("in_exon") and loc.get("nearest_boundary") is not None:
        bx = gx(loc["nearest_boundary"])
        vx = gx(pos)
        parts.append(
            f'<line x1="{vx:.1f}" y1="{center_y}" x2="{bx:.1f}" y2="{center_y}" '
            f'stroke="var(--warn, #b5532a)" stroke-width="1.2" '
            f'stroke-dasharray="3 2"/>'
        )

    # variant marker: stem + downward triangle
    vx = gx(pos)
    parts.append(
        f'<line x1="{vx:.1f}" y1="18" x2="{vx:.1f}" y2="{center_y + cds_h/2:.1f}" '
        f'stroke="var(--warn, #b5532a)" stroke-width="1.5"/>'
    )
    parts.append(
        f'<path d="M{vx - 5:.1f},8 L{vx + 5:.1f},8 L{vx:.1f},18 Z" '
        f'fill="var(--warn, #b5532a)"/>'
    )

    # 5' / 3' orientation labels
    parts.append(
        f'<text x="{ml}" y="{height - 6}" font-size="9" '
        f'fill="var(--ink-soft, #5a6275)">5\u2032</text>'
    )
    parts.append(
        f'<text x="{width - mr}" y="{height - 6}" font-size="9" '
        f'text-anchor="end" fill="var(--ink-soft, #5a6275)">3\u2032</text>'
    )

    parts.append('</svg>')
    return "".join(parts)


# ---------------------------------------------------------------------------
# Top-level entry point used by annotate.py
# ---------------------------------------------------------------------------
def build_gene_model(client, vep_parsed, pos, enabled=True):
    """Fetch structure, locate the variant, and render the SVG diagram."""
    result = {"found": False, "svg": "", "caption": "", "error": "",
              "location": None}
    if not enabled:
        result["error"] = "skipped"
        return result

    ens_tx = vep_parsed.get("ensembl_transcript_id")
    model = fetch_transcript_model(client, ens_tx)
    if not model.get("found"):
        result["error"] = model.get("error", "no transcript model")
        return result

    loc = locate_variant(model, pos)
    result["found"] = True
    result["location"] = loc
    result["caption"] = describe_location(loc)
    result["svg"] = render_svg(model, pos, loc)
    return result
