"""
alphagenome.py
==============
Optional AlphaGenome (Google DeepMind) variant-effect predictions.

When an API key is supplied (``--alphagenome-key`` / ``ALPHAGENOME_API_KEY``),
this runs ``predict_variant`` for each variant and renders, per requested
output type, a REF-vs-ALT figure centered on the variant. Figures are encoded
inline as base64 PNG data URIs so the report stays a single self-contained file.

The AlphaGenome client (``pip install alphagenome``) and matplotlib are imported
lazily, so the rest of varannot works without them installed.

Docs: https://www.alphagenomedocs.com  |  https://github.com/google-deepmind/alphagenome
"""

import base64
import hashlib
import io
import os

# All supported output types (the public OutputType enum names).
ALL_OUTPUT_TYPES = [
    "ATAC", "CAGE", "DNASE", "RNA_SEQ", "CHIP_HISTONE", "CHIP_TF",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
    "CONTACT_MAPS", "PROCAP",
]

# OutputType name -> attribute on the prediction Output object.
_ATTR = {
    "ATAC": "atac",
    "CAGE": "cage",
    "DNASE": "dnase",
    "RNA_SEQ": "rna_seq",
    "CHIP_HISTONE": "chip_histone",
    "CHIP_TF": "chip_tf",
    "SPLICE_SITES": "splice_sites",
    "SPLICE_SITE_USAGE": "splice_site_usage",
    "SPLICE_JUNCTIONS": "splice_junctions",
    "CONTACT_MAPS": "contact_maps",
    "PROCAP": "procap",
}

DEFAULT_SEQUENCE_LENGTH = "SEQUENCE_LENGTH_100KB"
DEFAULT_ZOOM = 2 ** 15  # ~32 kb window around the variant for 1D tracks
DEFAULT_TOP_TRACKS = 1   # plot only the N most variant-relevant tracks/tissues


def is_available():
    """True if the alphagenome package is importable."""
    try:
        import alphagenome  # noqa: F401
        return True
    except Exception:
        return False


# Bump when the plotting/track-selection logic changes so old PNGs are redrawn.
CACHE_VERSION = "v2"


def _cache_path(cache_dir, variant_str, out_name, length_name, zoom, ontology,
                top_n):
    raw = (f"{CACHE_VERSION}|{variant_str}|{out_name}|{length_name}|{zoom}|"
           f"{','.join(ontology or [])}|top{top_n}")
    h = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:32]
    return os.path.join(cache_dir, "alphagenome", f"{h}.png")


def _png_to_data_uri(png_bytes):
    return "data:image/png;base64," + base64.b64encode(png_bytes).decode("ascii")


def _top_track_indices(ref, alt, n):
    """
    Rank tracks by how strongly the variant changes them and return the indices
    of the top ``n`` (the tracks/tissues most relevant to this variant).

    When both REF and ALT are present we rank by the total absolute REF->ALT
    difference per track; otherwise we fall back to the strongest overall signal.
    Tracks with no signal at all (e.g. empty splice-junction tracks) are dropped,
    so they never fill up the plot -- unless *every* track is empty, in which
    case we keep the strongest few so the figure isn't blank.
    """
    import numpy as np
    if ref is not None and alt is not None and ref.num_tracks == alt.num_tracks:
        diff = np.abs(np.asarray(alt.values) - np.asarray(ref.values))
        score = diff.reshape(-1, diff.shape[-1]).sum(axis=0)
        # An informative track must carry some signal in REF or ALT (not just be
        # identical zeros); rank those, but fall back to the raw signal if the
        # variant changes nothing.
        signal = (np.abs(np.asarray(ref.values)).reshape(-1, ref.num_tracks)
                  .sum(axis=0)
                  + np.abs(np.asarray(alt.values)).reshape(-1, alt.num_tracks)
                  .sum(axis=0))
    else:
        base = ref if ref is not None else alt
        score = np.abs(np.asarray(base.values)).reshape(-1, base.num_tracks).sum(axis=0)
        signal = score

    if score.shape[0] == 0:
        return []

    # Prefer tracks that the variant actually perturbs; otherwise the ones with
    # any signal; otherwise (all empty) anything.
    rank = score
    candidates = np.where(score > 0)[0]
    if candidates.size == 0:
        rank = signal
        candidates = np.where(signal > 0)[0]
    if candidates.size == 0:
        candidates = np.arange(score.shape[0])

    order = candidates[np.argsort(rank[candidates])[::-1]]
    n = max(1, min(n, order.shape[0]))
    idx = order[:n]
    # Keep the original track order among the chosen ones for stable plotting.
    return sorted(int(i) for i in idx)


def _filter_tracks(data, idx):
    """Return ``data`` keeping only the track indices in ``idx``."""
    import numpy as np
    mask = np.zeros(data.num_tracks, dtype=bool)
    mask[idx] = True
    return data.filter_tracks(mask)


def _figure_for_output(plot_components, outputs, variant, out_name, zoom,
                       top_n=DEFAULT_TOP_TRACKS):
    """Build a matplotlib Figure for one output type.

    Returns ``(figure, error, track_name)``. AlphaGenome returns hundreds of
    tracks when no ontology filter is set (more than its plotter allows), so we
    keep only the ``top_n`` tracks where this variant has the largest effect.
    """
    attr = _ATTR[out_name]
    ref = getattr(outputs.reference, attr, None)
    alt = getattr(outputs.alternate, attr, None)
    if ref is None and alt is None:
        return None, "no data returned for this output", ""

    # Keep only the most variant-relevant tracks (drops empty/zero-signal tracks,
    # e.g. the many blank splice-junction tracks). Works for TrackData and
    # JunctionData alike via filter_tracks.
    idx = _top_track_indices(ref, alt, top_n)
    if not idx:
        return None, "no tracks with signal near the variant", ""
    if ref is not None:
        ref = _filter_tracks(ref, idx)
    if alt is not None:
        alt = _filter_tracks(alt, idx)

    annotations = [plot_components.VariantAnnotation([variant], alpha=0.8)]
    base = ref if ref is not None else alt  # at least one is non-None here

    track_name = ""
    try:
        names = list(base.names)
        if names:
            track_name = str(names[0])
    except Exception:
        track_name = ""

    if out_name == "CONTACT_MAPS":
        comps = []
        if ref is not None:
            comps.append(plot_components.ContactMaps(tdata=ref))
        if alt is not None:
            comps.append(plot_components.ContactMaps(tdata=alt))
        return (plot_components.plot(comps, interval=base.interval,
                                     annotations=annotations), None, track_name)

    if out_name == "SPLICE_JUNCTIONS":
        comps = []
        if ref is not None:
            comps.append(plot_components.Sashimi(ref))
        if alt is not None:
            comps.append(plot_components.Sashimi(alt))
        return (plot_components.plot(comps, interval=base.interval.resize(zoom),
                                     annotations=annotations), None, track_name)

    tdata = {}
    if ref is not None:
        tdata["REF"] = ref
    if alt is not None:
        tdata["ALT"] = alt
    comp = plot_components.OverlaidTracks(
        tdata=tdata, colors={"REF": "dimgrey", "ALT": "red"})
    return (plot_components.plot([comp], interval=base.interval.resize(zoom),
                                 annotations=annotations), None, track_name)


def run_for_variants(variants_meta, api_key, cache_dir=".varannot_cache",
                     output_types=None, ontology_terms=None,
                     sequence_length=DEFAULT_SEQUENCE_LENGTH, zoom=DEFAULT_ZOOM,
                     top_n=DEFAULT_TOP_TRACKS):
    """
    Run AlphaGenome for each variant and build per-output plot data URIs.

    Parameters
    ----------
    variants_meta : list of dicts with keys chrom, pos, ref, alt, label, gene
    api_key       : AlphaGenome API key
    output_types  : list of OutputType names (defaults to all)
    ontology_terms: optional list of ontology CURIEs (e.g. ['UBERON:0001157'])
    top_n         : per output, plot only the N most variant-relevant tracks

    Returns a dict consumed by the template:
      {enabled, error,
       variants:[{label, gene, id, plots:[{output,img,error,track}]}]}
    """
    result = {"enabled": True, "error": "", "variants": []}
    output_types = output_types or list(ALL_OUTPUT_TYPES)

    if not api_key:
        result["enabled"] = False
        result["error"] = "no API key"
        return result

    # Point matplotlib at a writable cache dir (avoids warnings / speeds import).
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(cache_dir, "mpl"))
    try:
        os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    except OSError:
        pass

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from alphagenome.data import genome
        from alphagenome.models import dna_client
        from alphagenome.visualization import plot_components
    except Exception as exc:
        result["enabled"] = False
        result["error"] = (f"AlphaGenome/matplotlib not installed ({exc}); "
                           f"pip install alphagenome")
        return result

    os.makedirs(os.path.join(cache_dir, "alphagenome"), exist_ok=True)

    try:
        seq_len = getattr(dna_client, sequence_length)
    except AttributeError:
        seq_len = getattr(dna_client, DEFAULT_SEQUENCE_LENGTH)

    try:
        model = dna_client.create(api_key)
    except Exception as exc:
        result["enabled"] = False
        result["error"] = f"could not create AlphaGenome client: {exc}"
        return result

    for i, vm in enumerate(variants_meta):
        chrom = str(vm["chrom"])
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        variant_str = f"{chrom}-{vm['pos']}-{vm['ref']}-{vm['alt']}"

        entry = {"label": vm.get("label", variant_str),
                 "gene": vm.get("gene", ""),
                 "id": f"ag-{i}", "plots": [], "error": ""}

        # Which outputs still need computing (not already cached as PNG)?
        cache_paths = {
            o: _cache_path(cache_dir, variant_str, o, sequence_length, zoom,
                           ontology_terms, top_n)
            for o in output_types
        }
        missing = [o for o in output_types if not os.path.exists(cache_paths[o])]

        outputs = None
        if missing:
            try:
                variant = genome.Variant(
                    chromosome=chrom, position=int(vm["pos"]),
                    reference_bases=vm["ref"], alternate_bases=vm["alt"])
                interval = variant.reference_interval.resize(seq_len)
                requested = {dna_client.OutputType[o] for o in missing}
                outputs = model.predict_variant(
                    interval=interval, variant=variant,
                    requested_outputs=requested,
                    ontology_terms=ontology_terms)
            except Exception as exc:
                entry["error"] = f"prediction failed: {exc}"

        for o in output_types:
            path = cache_paths[o]
            track_path = path + ".track"
            plot = {"output": o, "img": "", "error": "", "track": ""}

            if os.path.exists(path):
                with open(path, "rb") as fh:
                    plot["img"] = _png_to_data_uri(fh.read())
                if os.path.exists(track_path):
                    with open(track_path, "r", encoding="utf-8") as fh:
                        plot["track"] = fh.read().strip()
                entry["plots"].append(plot)
                continue

            if outputs is None:
                plot["error"] = entry["error"] or "not computed"
                entry["plots"].append(plot)
                continue

            try:
                fig, err, track_name = _figure_for_output(
                    plot_components, outputs, variant, o, zoom, top_n=top_n)
                if fig is None:
                    plot["error"] = err or "no figure"
                else:
                    buf = io.BytesIO()
                    fig.savefig(buf, format="png", dpi=110,
                                bbox_inches="tight")
                    plt.close(fig)
                    png = buf.getvalue()
                    with open(path, "wb") as fh:
                        fh.write(png)
                    plot["img"] = _png_to_data_uri(png)
                    plot["track"] = track_name or ""
                    if track_name:
                        with open(track_path, "w", encoding="utf-8") as fh:
                            fh.write(track_name)
            except Exception as exc:
                plot["error"] = f"plot failed: {exc}"

            entry["plots"].append(plot)

        result["variants"].append(entry)

    return result
