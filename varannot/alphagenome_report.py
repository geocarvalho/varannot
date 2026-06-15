#!/usr/bin/env python3
"""
varannot AlphaGenome report
===========================
Standalone command that generates a *separate* HTML report with AlphaGenome
(Google DeepMind) variant-effect plots. This is kept apart from the main
annotation report because the predictions are slow and require an API key.

Usage:
    python -m varannot.alphagenome_report \\
        --input variants.csv \\
        --output alphagenome_report.html \\
        --alphagenome-key YOUR_KEY

    # Key can also come from the ALPHAGENOME_API_KEY environment variable.
    # Requires: pip install alphagenome

Output types (default: all):
    ATAC, CAGE, DNASE, RNA_SEQ, CHIP_HISTONE, CHIP_TF, SPLICE_SITES,
    SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS, PROCAP
"""

import argparse
import datetime
import os
import sys

from jinja2 import Environment, FileSystemLoader, select_autoescape

from .annotate import parse_variants
from .sources import alphagenome as alphagenome_src

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")


def render_report(alphagenome, output_path):
    env = Environment(
        loader=FileSystemLoader(TEMPLATE_DIR),
        autoescape=select_autoescape(["html", "xml", "j2"]),
    )
    template = env.get_template("alphagenome_report.html.j2")
    html = template.render(
        alphagenome=alphagenome,
        n_variants=len(alphagenome.get("variants", [])),
        generated=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    )
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def main():
    parser = argparse.ArgumentParser(
        description="AlphaGenome variant-effect plots -> standalone HTML report")
    parser.add_argument("--input", required=True,
                        help="CSV/TSV of variants (chr,pos,ref,alt)")
    parser.add_argument("--output", default="alphagenome_report.html",
                        help="Output HTML file")
    parser.add_argument("--cache-dir", default=".varannot_cache",
                        help="Directory for cached AlphaGenome plot PNGs")
    parser.add_argument("--alphagenome-key",
                        default=os.environ.get("ALPHAGENOME_API_KEY"),
                        help="AlphaGenome API key (or set ALPHAGENOME_API_KEY). "
                             "Requires `pip install alphagenome`.")
    parser.add_argument("--alphagenome-ontology", default=None,
                        help="Comma-separated ontology CURIEs to restrict tracks "
                             "(e.g. 'UBERON:0001157,UBERON:0002048'). Default: all.")
    parser.add_argument("--alphagenome-length", default="SEQUENCE_LENGTH_100KB",
                        help="Input window (SEQUENCE_LENGTH_16KB|100KB|500KB|1MB)")
    parser.add_argument("--alphagenome-outputs", default=None,
                        help="Comma-separated output types (default: all). "
                             "E.g. 'RNA_SEQ,DNASE,ATAC'.")
    parser.add_argument("--alphagenome-top-tracks", type=int,
                        default=alphagenome_src.DEFAULT_TOP_TRACKS,
                        help="Per output, plot only the N tracks/tissues where "
                             "the variant has the largest effect (default 1). "
                             "Avoids AlphaGenome's 'too many tracks' limit when "
                             "no --alphagenome-ontology filter is given.")
    parser.add_argument("--max-variants", type=int, default=50,
                        help="Safety cap on number of variants (default 50)")
    args = parser.parse_args()

    if not args.alphagenome_key:
        print("ERROR: an AlphaGenome API key is required.")
        print("       Pass --alphagenome-key or set ALPHAGENOME_API_KEY.")
        sys.exit(1)

    if not alphagenome_src.is_available():
        print("ERROR: the 'alphagenome' package is not installed.")
        print("       Install it with: pip install alphagenome")
        sys.exit(1)

    variants = parse_variants(args.input)
    if not variants:
        print("ERROR: no valid variants found in input.")
        sys.exit(1)
    if len(variants) > args.max_variants:
        print(f"ERROR: {len(variants)} variants exceeds "
              f"--max-variants={args.max_variants}.")
        sys.exit(1)

    ontology = ([t.strip() for t in args.alphagenome_ontology.split(",")
                 if t.strip()] if args.alphagenome_ontology else None)
    outputs = ([o.strip().upper() for o in args.alphagenome_outputs.split(",")
                if o.strip()] if args.alphagenome_outputs else None)

    variants_meta = [{
        "chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
        "label": f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}",
        "gene": "",
    } for v in variants]

    print(f"Running AlphaGenome for {len(variants_meta)} variant(s) "
          f"(this can take a while per variant)...")
    alphagenome = alphagenome_src.run_for_variants(
        variants_meta, args.alphagenome_key, cache_dir=args.cache_dir,
        output_types=outputs, ontology_terms=ontology,
        sequence_length=args.alphagenome_length,
        top_n=args.alphagenome_top_tracks)

    if not alphagenome.get("enabled"):
        print(f"ERROR: AlphaGenome did not run ({alphagenome.get('error')}).")
        sys.exit(1)

    render_report(alphagenome, args.output)
    print(f"\nDone! AlphaGenome report written to: {args.output}")
    print(f"Plot cache stored in: {args.cache_dir}/alphagenome/")


if __name__ == "__main__":
    main()
