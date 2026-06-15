"""
gnomad.py
=========
gnomAD v4 GraphQL client.

Provides:
  - Per-variant allele counts, homozygote/heterozygote counts (exome + genome)
  - Gene constraint metrics: pLI, missense Z, pRec, LOEUF
  - A browser link for the variant

Endpoint: https://gnomad.broadinstitute.org/api
"""

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_BROWSER = "https://gnomad.broadinstitute.org/variant"

VARIANT_QUERY = """
query VariantInfo($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    exome {
      ac
      an
      af
      ac_hom
    }
    genome {
      ac
      an
      af
      ac_hom
    }
  }
}
"""

CONSTRAINT_QUERY = """
query GeneConstraint($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    gene_id
    symbol
    gnomad_constraint {
      pli
      mis_z
      syn_z
      lof_z
      oe_lof
      oe_lof_lower
      oe_lof_upper
      oe_mis
      exp_lof
      obs_lof
    }
  }
}
"""


def _variant_id(chrom, pos, ref, alt):
    chrom = str(chrom).replace("chr", "")
    return f"{chrom}-{pos}-{ref}-{alt}"


def query_variant(client, chrom, pos, ref, alt, dataset="gnomad_r4"):
    """Query gnomAD for variant-level allele/homozygote counts."""
    vid = _variant_id(chrom, pos, ref, alt)
    payload = {
        "query": VARIANT_QUERY,
        "variables": {"variantId": vid, "dataset": dataset},
    }
    cache_key = f"gnomad_var:{dataset}:{vid}"
    data = client.post_json(GNOMAD_API, payload, cache_key=cache_key)

    result = {
        "variant_id": vid,
        "browser_url": f"{GNOMAD_BROWSER}/{vid}?dataset={dataset}",
        "found": False,
        "exome": None,
        "genome": None,
        "total_ac": 0,
        "total_an": 0,
        "total_hom": 0,
        "total_het": 0,
        "af": None,
    }

    if "_error" in data:
        result["error"] = data["_error"]
        return result

    variant = (data.get("data") or {}).get("variant")
    if not variant:
        # Not observed in gnomAD
        return result

    result["found"] = True
    exome = variant.get("exome")
    genome = variant.get("genome")
    result["exome"] = exome
    result["genome"] = genome

    ac = an = hom = 0
    for block in (exome, genome):
        if block:
            ac += block.get("ac") or 0
            an += block.get("an") or 0
            hom += block.get("ac_hom") or 0

    result["total_ac"] = ac
    result["total_an"] = an
    result["total_hom"] = hom
    # heterozygous carriers = allele count - 2 * homozygotes
    result["total_het"] = max(ac - 2 * hom, 0)
    result["af"] = (ac / an) if an else None
    return result


def query_constraint(client, gene_symbol, dataset="gnomad_r4"):
    """Query gnomAD for gene-level constraint metrics."""
    if not gene_symbol:
        return {"found": False}
    payload = {
        "query": CONSTRAINT_QUERY,
        "variables": {"geneSymbol": gene_symbol, "referenceGenome": "GRCh38"},
    }
    cache_key = f"gnomad_constraint:{gene_symbol}"
    data = client.post_json(GNOMAD_API, payload, cache_key=cache_key)

    result = {"found": False, "gene_symbol": gene_symbol}
    if "_error" in data:
        result["error"] = data["_error"]
        return result

    gene = (data.get("data") or {}).get("gene")
    if not gene:
        return result
    constraint = gene.get("gnomad_constraint")
    if not constraint:
        return result

    result["found"] = True
    result["gene_id"] = gene.get("gene_id")
    if result["gene_id"]:
        result["browser_url"] = (f"https://gnomad.broadinstitute.org/gene/"
                                 f"{result['gene_id']}?dataset={dataset}")
    else:
        result["browser_url"] = (f"https://gnomad.broadinstitute.org/gene/"
                                 f"{gene_symbol}?dataset={dataset}")
    result["pli"] = constraint.get("pli")
    result["mis_z"] = constraint.get("mis_z")
    result["syn_z"] = constraint.get("syn_z")
    result["lof_z"] = constraint.get("lof_z")
    result["oe_lof"] = constraint.get("oe_lof")
    result["oe_lof_upper"] = constraint.get("oe_lof_upper")  # this is LOEUF
    result["oe_mis"] = constraint.get("oe_mis")
    # Note: gnomAD v4 constraint no longer exposes pRec directly; LOEUF preferred
    return result
