#!/usr/bin/env python3
"""
varannot — variant annotation framework
========================================
Given a list of variants (chr,pos,ref,alt) produce a single combined HTML
report with HGVS, conservation, gnomAD counts, AlphaMissense, ClinVar,
gene constraint, and OMIM.

Usage:
    python -m varannot.annotate \\
        --input variants.csv \\
        --output report.html \\
        --omim-key YOUR_KEY

    # OMIM key can also come from the OMIM_API_KEY environment variable.

Input format (one variant per line, header optional):
    chr3,177026399,C,T
    chrX,32389644,G,A

Requires: requests, jinja2
"""

import argparse
import datetime
import os
import sys

from jinja2 import Environment, FileSystemLoader, select_autoescape

from .http_client import CachedSession
from .sources import vep as vep_src
from .sources import gnomad as gnomad_src
from .sources import clinvar as clinvar_src
from .sources import omim as omim_src
from .sources import conservation as cons_src
from .sources import spliceai as spliceai_src
from .sources import revel as revel_src
from .sources import gene_model as gene_model_src
from .sources import protein_domains as protein_src
from .sources import gtex as gtex_src
from .sources import alphagenome as alphagenome_src


TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")


def parse_variants(path):
    """Read a CSV/TSV of variants. Accepts comma or tab; skips header if present."""
    variants = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            sep = "," if "," in line else "\t"
            parts = [p.strip() for p in line.split(sep)]
            if len(parts) < 4:
                print(f"  WARNING: skipping malformed line: {line}")
                continue
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            # Skip a header row like "chr,pos,ref,alt"
            if not pos.isdigit():
                continue
            variants.append({"chrom": chrom, "pos": int(pos), "ref": ref.upper(), "alt": alt.upper()})
    return variants


def annotate_one(client, var, omim_key, ncbi_key=None,
                 spliceai_enabled=False, spliceai_url=None,
                 cons46_enabled=False, exonaa_path=None, liftover=None,
                 omim_index=None):
    """Run all sources for a single variant and return a merged record."""
    chrom, pos, ref, alt = var["chrom"], var["pos"], var["ref"], var["alt"]
    print(f"  {chrom}:{pos} {ref}>{alt} ...", end="", flush=True)

    vep_raw = vep_src.query_vep(client, chrom, pos, ref, alt)
    vep = vep_src.parse_vep(vep_raw)
    gene = vep.get("gene_symbol", "")

    gnomad = gnomad_src.query_variant(client, chrom, pos, ref, alt)
    constraint = gnomad_src.query_constraint(client, gene) if gene else {"found": False}
    clinvar = clinvar_src.query_clinvar(client, chrom, pos, ref, alt, api_key=ncbi_key)
    if gene:
        if omim_index is not None:
            omim = omim_src.query_omim_local(omim_index, gene)
        elif omim_key:
            omim = omim_src.query_omim_api(client, gene, omim_key)
        else:
            omim = {"found": False}
    else:
        omim = {"found": False}

    conservation = cons_src.build_protein_track(
        client, vep, enabled=cons46_enabled,
        exonaa_path=exonaa_path, liftover=liftover,
        chrom=chrom, pos=pos,
    )
    conservation_score = cons_src.conservation_score_from_vep(vep)

    spliceai = spliceai_src.query_spliceai(
        client, chrom, pos, ref, alt,
        base_url=spliceai_url, enabled=spliceai_enabled,
    )
    spliceai["interpretation"] = spliceai_src.interpret(spliceai.get("max_delta"))

    revel = revel_src.query_revel(client, chrom, pos, ref, alt)
    revel["interpretation"] = revel_src.interpret(revel.get("revel"))

    gene_model = gene_model_src.build_gene_model(client, vep, pos)
    protein_diagram = protein_src.build_protein_diagram(client, vep)
    gtex = gtex_src.query_gtex(client, gene, gene_id=vep.get("gene_id"))

    print(" done")
    return {
        "input": var,
        "vep": vep,
        "gnomad": gnomad,
        "constraint": constraint,
        "clinvar": clinvar,
        "omim": omim,
        "conservation": conservation,
        "conservation_score": conservation_score,
        "spliceai": spliceai,
        "revel": revel,
        "gene_model": gene_model,
        "protein_diagram": protein_diagram,
        "gtex": gtex,
    }


def render_report(records, output_path, alphagenome=None):
    """Render the combined HTML report.

    When ``alphagenome`` (the dict from ``alphagenome_src.run_for_variants``) is
    supplied, the report is rendered with two tabs: Annotations + AlphaGenome.
    Otherwise it's the plain single-view annotation report.
    """
    env = Environment(
        loader=FileSystemLoader(TEMPLATE_DIR),
        autoescape=select_autoescape(["html", "xml", "j2"]),
    )
    template = env.get_template("report.html.j2")
    html = template.render(
        variants=records,
        n_variants=len(records),
        generated=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        alphagenome=alphagenome,
    )
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def main():
    parser = argparse.ArgumentParser(description="Variant annotation -> combined HTML report")
    parser.add_argument("--input", required=True, help="CSV/TSV of variants (chr,pos,ref,alt)")
    parser.add_argument("--output", default="variant_report.html", help="Output HTML file")
    parser.add_argument("--omim-genemap2", default=os.environ.get("OMIM_GENEMAP2"),
                        help="Path to OMIM genemap2.txt download file (preferred; no key needed)")
    parser.add_argument("--omim-key", default=os.environ.get("OMIM_API_KEY"),
                        help="OMIM API key (alternative to --omim-genemap2; or set OMIM_API_KEY)")
    parser.add_argument("--ncbi-key", default=os.environ.get("NCBI_API_KEY"),
                        help="Optional NCBI E-utilities API key (raises rate limit)")
    parser.add_argument("--cache-dir", default=".varannot_cache",
                        help="Directory for cached API responses")
    parser.add_argument("--min-interval", type=float, default=0.4,
                        help="Minimum seconds between live API calls (politeness)")
    parser.add_argument("--spliceai", action="store_true",
                        help="Enable SpliceAI scores (public endpoint is rate limited to ~few/min)")
    parser.add_argument("--spliceai-url", default=None,
                        help="Custom SpliceAI API base URL (e.g. your own local instance)")
    parser.add_argument("--conservation46", action="store_true",
                        help="Build the UCSC multiz46way amino-acid track (downloads ~355MB once, hg38->hg19)")
    parser.add_argument("--exonaa-dir", default=".varannot_cache",
                        help="Where to store/find the downloaded 46-way exonAA file")
    parser.add_argument("--max-variants", type=int, default=50,
                        help="Safety cap on number of variants (default 50)")
    parser.add_argument("--alphagenome-key",
                        default=os.environ.get("ALPHAGENOME_API_KEY"),
                        help="AlphaGenome API key. When given, the report gains a "
                             "second 'AlphaGenome' tab (combined report). Or set "
                             "ALPHAGENOME_API_KEY. Requires `pip install alphagenome`.")
    parser.add_argument("--alphagenome-ontology", default=None,
                        help="Comma-separated ontology CURIEs to restrict "
                             "AlphaGenome tracks. Default: all.")
    parser.add_argument("--alphagenome-length", default="SEQUENCE_LENGTH_100KB",
                        help="AlphaGenome input window "
                             "(SEQUENCE_LENGTH_16KB|100KB|500KB|1MB)")
    parser.add_argument("--alphagenome-outputs", default=None,
                        help="Comma-separated AlphaGenome output types "
                             "(default: all). E.g. 'RNA_SEQ,DNASE,ATAC'.")
    parser.add_argument("--alphagenome-top-tracks", type=int,
                        default=alphagenome_src.DEFAULT_TOP_TRACKS,
                        help="Per output, plot only the N tracks/tissues where the "
                             "variant has the largest effect (default 1).")
    args = parser.parse_args()

    # OMIM: prefer the local genemap2.txt download; fall back to API key.
    omim_index = None
    if args.omim_genemap2:
        print(f"Loading OMIM genemap2.txt from {args.omim_genemap2} ...")
        omim_index = omim_src.load_genemap2(args.omim_genemap2)
        if omim_index is None:
            print("  WARNING: could not read genemap2.txt; OMIM section will be skipped.\n")
        else:
            print(f"  Loaded {len(omim_index)} gene symbols.\n")
    elif not args.omim_key:
        print("NOTE: no OMIM source provided — OMIM section will be skipped.")
        print("      Pass --omim-genemap2 /path/to/genemap2.txt (preferred) or --omim-key.\n")

    variants = parse_variants(args.input)
    if not variants:
        print("ERROR: no valid variants found in input.")
        sys.exit(1)

    if len(variants) > args.max_variants:
        print(f"ERROR: {len(variants)} variants exceeds --max-variants={args.max_variants}.")
        print("       This tool is designed for small batches; raise the cap if intended.")
        sys.exit(1)

    if args.spliceai and not args.spliceai_url:
        print("NOTE: SpliceAI uses the public Broad endpoint (rate limited to a few/min).")
        print(f"      With a ~6s delay per variant, {len(variants)} variants will take "
              f"about {len(variants) * 6 // 60 + 1} min.")
        print("      For speed, run your own instance and pass --spliceai-url.\n")

    print(f"Annotating {len(variants)} variant(s)...")
    client = CachedSession(cache_dir=args.cache_dir, min_interval=args.min_interval)

    # One-time 46-way setup: download the exonAA alignment. The alignment is
    # indexed by RefSeq accession + protein residue, so no genomic liftover is
    # needed (the residue is found directly from the MANE RefSeq id).
    exonaa_path = None
    liftover = None
    if args.conservation46:
        print("\nSetting up UCSC multiz46way conservation track...")
        print("  (residues are matched by MANE RefSeq id + protein position)")
        try:
            exonaa_path = cons_src.ensure_exonaa(args.exonaa_dir)
        except Exception as exc:
            print(f"  WARNING: 46-way setup failed ({exc}); falling back to phyloP score only.")
            exonaa_path = None
        print("")

    records = []
    for var in variants:
        try:
            records.append(annotate_one(
                client, var, args.omim_key, args.ncbi_key,
                spliceai_enabled=args.spliceai, spliceai_url=args.spliceai_url,
                cons46_enabled=args.conservation46,
                exonaa_path=exonaa_path, liftover=liftover,
                omim_index=omim_index,
            ))
        except Exception as exc:  # keep going on per-variant failures
            print(f" ERROR: {exc}")
            records.append({
                "input": var,
                "vep": {"error": str(exc)},
                "gnomad": {"found": False, "browser_url": "#", "error": str(exc)},
                "constraint": {"found": False},
                "clinvar": {"found": False, "error": str(exc)},
                "omim": {"found": False, "error": str(exc)},
                "conservation": {"available": False, "note": str(exc),
                                 "species": [], "track": "", "human_aa": "",
                                 "n_species": 0, "assembly": "hg19/46way"},
                "conservation_score": None,
                "spliceai": {"enabled": args.spliceai, "found": False,
                             "error": str(exc), "transcripts": [], "max_delta": None,
                             "web_url": "", "interpretation": ""},
                "revel": {"found": False, "error": str(exc), "revel": None,
                          "cadd_phred": None, "url": "", "interpretation": ""},
                "gene_model": {"found": False, "svg": "", "caption": "",
                               "error": str(exc), "location": None},
                "protein_diagram": {"found": False, "svg": "", "caption": "",
                                    "error": str(exc), "domains": [],
                                    "length": None},
                "gtex": {"found": False, "tissues": [], "max_median": 0.0,
                         "error": str(exc), "url": ""},
            })

    # Optional: combined report with an AlphaGenome tab (only when a key is given).
    alphagenome = None
    if args.alphagenome_key:
        if not alphagenome_src.is_available():
            print("\nNOTE: --alphagenome-key given but the 'alphagenome' package "
                  "is not installed; skipping the AlphaGenome tab.")
            print("      Install it with: pip install alphagenome matplotlib")
        else:
            print("\nRunning AlphaGenome predictions for the combined report "
                  "(this can take a while per variant)...")
            ontology = ([t.strip() for t in args.alphagenome_ontology.split(",")
                         if t.strip()] if args.alphagenome_ontology else None)
            outputs = ([o.strip().upper() for o in args.alphagenome_outputs.split(",")
                        if o.strip()] if args.alphagenome_outputs else None)
            variants_meta = [{
                "chrom": rec["input"]["chrom"], "pos": rec["input"]["pos"],
                "ref": rec["input"]["ref"], "alt": rec["input"]["alt"],
                "label": (f"{rec['input']['chrom']}:{rec['input']['pos']} "
                          f"{rec['input']['ref']}>{rec['input']['alt']}"),
                "gene": rec.get("vep", {}).get("gene_symbol", ""),
            } for rec in records]
            alphagenome = alphagenome_src.run_for_variants(
                variants_meta, args.alphagenome_key, cache_dir=args.cache_dir,
                output_types=outputs, ontology_terms=ontology,
                sequence_length=args.alphagenome_length,
                top_n=args.alphagenome_top_tracks)
            if alphagenome.get("error"):
                print(f"  AlphaGenome note: {alphagenome['error']}")
            else:
                print(f"  AlphaGenome plots generated for "
                      f"{len(alphagenome['variants'])} variant(s).")

    render_report(records, args.output, alphagenome=alphagenome)
    print(f"\nDone! Report written to: {args.output}")
    print(f"Cache stored in: {args.cache_dir}/ (delete to force fresh queries)")


if __name__ == "__main__":
    main()
