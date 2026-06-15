"""
protein_domains.py
==================
Build a "where in the protein" diagram for each (coding) variant.

It fetches Pfam domain annotations and the protein length for the variant's
Ensembl protein (ENSP) from the Ensembl REST API, then renders an inline SVG:
the full-length protein as a bar, Pfam domains as labelled colored boxes, and
the variant's residue marked with a pointer. For truncating variants
(stop-gain / frameshift) the lost C-terminal portion is shaded.

Like the gene model, the SVG is inline (no image files, no extra deps) and uses
the report's CSS variables, with hardcoded color fallbacks for standalone use.
"""

ENSEMBL_SERVER = "https://rest.ensembl.org"

# Distinct, muted colors for domains (cycled by unique domain name).
_PALETTE = [
    "#2b6f6a", "#8a5a2b", "#3f5d9e", "#7a8a2b",
    "#8a2b6a", "#2b8a8a", "#9e6a3f", "#5a3f9e",
]

_TRUNCATING = ("stop_gained", "frameshift_variant", "start_lost",
               "stop_lost", "transcript_ablation")


def fetch_protein_length(client, ensp):
    tx = ensp.split(".")[0]
    data = client.get_json(f"{ENSEMBL_SERVER}/lookup/id/{tx}",
                           cache_key=f"ens_translation:{tx}")
    if isinstance(data, dict) and "_error" not in data:
        return data.get("length")
    return None


def fetch_pfam_domains(client, ensp):
    tx = ensp.split(".")[0]
    data = client.get_json(f"{ENSEMBL_SERVER}/overlap/translation/{tx}",
                           params={"type": "Pfam"},
                           cache_key=f"ens_pfam:{tx}")
    feats = []
    if isinstance(data, list):
        for f in data:
            try:
                start, end = int(f["start"]), int(f["end"])
            except (KeyError, TypeError, ValueError):
                continue
            name = f.get("description") or f.get("id") or "domain"
            feats.append({"start": start, "end": end, "name": name,
                          "id": f.get("id", "")})
    feats.sort(key=lambda d: d["start"])
    return feats


def _assign_colors(domains):
    """Map each unique domain name to a stable color."""
    colors = {}
    for d in domains:
        if d["name"] not in colors:
            colors[d["name"]] = _PALETTE[len(colors) % len(_PALETTE)]
        d["color"] = colors[d["name"]]
    return colors


def _locate(domains, aa_pos):
    """Find the domain containing aa_pos, else the nearest one."""
    for d in domains:
        if d["start"] <= aa_pos <= d["end"]:
            return {"in_domain": True, "name": d["name"], "distance": 0}
    best = None
    for d in domains:
        dist = min(abs(aa_pos - d["start"]), abs(aa_pos - d["end"]))
        if best is None or dist < best["distance"]:
            best = {"in_domain": False, "name": d["name"], "distance": dist}
    return best


def _render_svg(length, domains, aa_pos, truncating, width=560):
    ml, mr = 16, 16
    plot_w = width - ml - mr
    height = 64
    bar_y = 30
    bar_h = 12
    dom_h = 18

    def ax(aa):
        frac = (aa - 1) / max(length - 1, 1)
        return ml + frac * plot_w

    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="100%" '
        f'preserveAspectRatio="xMidYMid meet" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'font-family="var(--sans, sans-serif)" role="img" '
        f'aria-label="protein domain diagram showing variant position">'
    ]

    # backbone bar
    parts.append(
        f'<rect x="{ml}" y="{bar_y}" width="{plot_w:.1f}" height="{bar_h}" '
        f'rx="3" fill="var(--line, #e3e0d8)"/>'
    )

    # truncated (lost) C-terminal region for stop/frameshift variants
    if truncating and aa_pos and aa_pos < length:
        x0 = ax(aa_pos)
        parts.append(
            f'<rect x="{x0:.1f}" y="{bar_y}" width="{(ml + plot_w) - x0:.1f}" '
            f'height="{bar_h}" rx="3" fill="var(--warn, #b5532a)" '
            f'opacity="0.18"/>'
        )

    # domains
    cy = bar_y + bar_h / 2
    for d in domains:
        x1, x2 = ax(d["start"]), ax(d["end"])
        w = max(x2 - x1, 2.0)
        parts.append(
            f'<rect x="{x1:.1f}" y="{cy - dom_h/2:.1f}" width="{w:.1f}" '
            f'height="{dom_h}" rx="3" fill="{d["color"]}" opacity="0.92">'
            f'<title>{_esc(d["name"])} ({d["start"]}\u2013{d["end"]})</title>'
            f'</rect>'
        )

    # variant marker
    if aa_pos:
        vx = ax(aa_pos)
        parts.append(
            f'<line x1="{vx:.1f}" y1="6" x2="{vx:.1f}" y2="{cy + dom_h/2 + 2:.1f}" '
            f'stroke="var(--warn, #b5532a)" stroke-width="1.5"/>'
        )
        parts.append(
            f'<path d="M{vx - 5:.1f},2 L{vx + 5:.1f},2 L{vx:.1f},11 Z" '
            f'fill="var(--warn, #b5532a)"/>'
        )

    # N / C terminus labels
    parts.append(
        f'<text x="{ml}" y="{height - 6}" font-size="9" '
        f'fill="var(--ink-soft, #5a6275)">N</text>'
    )
    parts.append(
        f'<text x="{width - mr}" y="{height - 6}" font-size="9" '
        f'text-anchor="end" fill="var(--ink-soft, #5a6275)">'
        f'C ({length} aa)</text>'
    )

    parts.append('</svg>')
    return "".join(parts)


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def _caption(aa_pos, length, loc, truncating):
    # Plain text only: the template renders this with autoescaping, so HTML
    # entities like &middot; would be shown literally. Use the Unicode middot.
    dot = " \u00b7 "
    pos = f"aa {aa_pos:,} of {length:,}"
    if truncating and aa_pos:
        lost = length - aa_pos
        pct = 100 * lost / length if length else 0
        trunc = f"{dot}truncates ~{lost:,} aa ({pct:.0f}% of protein)"
    else:
        trunc = ""
    if loc is None:
        return f"{pos}{trunc}{dot}no annotated Pfam domains"
    if loc["in_domain"]:
        return f"{pos}{trunc}{dot}in {loc['name']}"
    return (f"{pos}{trunc}{dot}not in a Pfam domain "
            f"(nearest: {loc['name']}, {loc['distance']:,} aa away)")


def build_protein_diagram(client, vep_parsed, enabled=True):
    """Fetch domains + length, locate the residue, and render the SVG."""
    result = {"found": False, "svg": "", "caption": "", "error": "",
              "domains": [], "length": None}
    if not enabled:
        result["error"] = "skipped"
        return result

    aa_pos = vep_parsed.get("protein_start")
    ensp = vep_parsed.get("ensembl_protein_id")
    if not aa_pos:
        result["error"] = "non-coding variant (no protein residue)"
        return result
    if not ensp:
        result["error"] = "no Ensembl protein id"
        return result
    aa_pos = int(aa_pos)

    length = fetch_protein_length(client, ensp)
    if not length:
        result["error"] = "protein length unavailable"
        return result

    domains = fetch_pfam_domains(client, ensp)
    _assign_colors(domains)
    loc = _locate(domains, aa_pos) if domains else None

    terms = set(vep_parsed.get("consequence_terms") or [])
    if vep_parsed.get("most_severe_consequence"):
        terms.add(vep_parsed["most_severe_consequence"])
    truncating = bool(terms & set(_TRUNCATING))

    result["found"] = True
    result["length"] = length
    result["caption"] = _caption(aa_pos, length, loc, truncating)
    result["svg"] = _render_svg(length, domains, aa_pos, truncating)
    # de-duplicated legend (unique domain names, in N->C order of first occurrence)
    seen = {}
    for d in domains:
        seen.setdefault(d["name"], d["color"])
    result["domains"] = [{"name": n, "color": c} for n, c in seen.items()]
    return result
