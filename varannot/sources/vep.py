"""
vep.py
======
Ensembl VEP REST client (GRCh38).

Provides:
  - HGVS genomic/coding/protein nomenclature
  - Most severe consequence and gene/transcript info
  - AlphaMissense score and class (via VEP plugin output)
  - Conservation (PhyloP / GERP) when available

Endpoint: https://rest.ensembl.org/vep/human/region
"""

VEP_SERVER = "https://rest.ensembl.org"


def query_vep(client, chrom, pos, ref, alt):
    """
    Query VEP for a single SNV/indel.

    Input coordinates are 1-based VCF style. Returns the raw VEP record
    (first element of the response list) or a dict with _error.
    """
    chrom = str(chrom).replace("chr", "")
    # VEP region input format: "chr pos . ref alt"
    variant_str = f"{chrom} {pos} . {ref} {alt} . . ."

    url = f"{VEP_SERVER}/vep/human/region"
    payload = {
        "variants": [variant_str],
        "hgvs": True,
        "AlphaMissense": True,
        "Conservation": True,
        "canonical": True,
        "mane": True,
        "numbers": True,
        "domains": True,
        # "merged" reports BOTH Ensembl and RefSeq transcripts, so we can show
        # RefSeq (NM_/NP_) HGVS while keeping Ensembl-only annotations like
        # AlphaMissense.
        "merged": True,
    }
    # Cache key is versioned ("m" = merged transcript set) so older cached
    # Ensembl-only responses are not reused after this change.
    cache_key = f"vep:m:{chrom}:{pos}:{ref}:{alt}"
    data = client.post_json(url, payload, cache_key=cache_key)

    if isinstance(data, dict) and "_error" in data:
        return data
    if isinstance(data, list) and data:
        return data[0]
    return {"_error": "empty VEP response"}


def _accession_from_hgvs(hgvs):
    """Pull the accession from an HGVS string like 'ENSP00000413251.3:p.Val...'."""
    if not hgvs or ":" not in hgvs:
        return ""
    return hgvs.split(":", 1)[0]


def _is_refseq(tc):
    """True if this transcript consequence is a RefSeq transcript."""
    src = (tc.get("source") or "").lower()
    tid = tc.get("transcript_id") or ""
    return src == "refseq" or tid.startswith(("NM_", "NR_", "XM_", "XR_"))


def _is_mane(tc):
    return bool(tc.get("mane_select") or tc.get("mane_plus_clinical"))


def _pick_display_transcript(tcs):
    """
    Pick the transcript whose HGVS we display, preferring RefSeq + MANE.

    Order: RefSeq MANE -> RefSeq canonical -> any MANE -> any canonical -> first.
    """
    if not tcs:
        return None
    preferences = [
        lambda tc: _is_refseq(tc) and _is_mane(tc),
        lambda tc: _is_refseq(tc) and tc.get("canonical"),
        lambda tc: _is_refseq(tc),
        _is_mane,
        lambda tc: tc.get("canonical"),
    ]
    for pred in preferences:
        for tc in tcs:
            if pred(tc):
                return tc
    return tcs[0]


def _pick_info_transcript(tcs):
    """
    Pick the Ensembl transcript used for gene/protein/AlphaMissense fields.

    AlphaMissense and some other annotations are only attached to Ensembl
    transcripts, so this prefers Ensembl MANE -> Ensembl canonical -> first
    Ensembl -> first.
    """
    if not tcs:
        return None
    ens = [tc for tc in tcs if not _is_refseq(tc)]
    for tc in ens:
        if _is_mane(tc):
            return tc
    for tc in ens:
        if tc.get("canonical"):
            return tc
    if ens:
        return ens[0]
    return tcs[0]


def _mane_refseq_accession(tcs):
    """Return the MANE Select RefSeq mRNA accession (e.g. 'NM_024665.7')."""
    for tc in tcs:
        if _is_refseq(tc) and _is_mane(tc):
            return tc.get("transcript_id", "")
    for tc in tcs:
        ms = tc.get("mane_select") or ""
        if ms.startswith(("NM_", "NR_")):
            return ms
    return ""


# Backwards-compatible alias.
def _pick_transcript(vep_record):
    return _pick_info_transcript(vep_record.get("transcript_consequences", []))


def _empty_parsed():
    """Return a parsed-VEP dict with every key present and empty."""
    return {
        "most_severe_consequence": "",
        "hgvsg": "",
        "hgvsc": "",
        "hgvsp": "",
        "gene_symbol": "",
        "gene_id": "",
        "transcript_id": "",
        "ensembl_transcript_id": "",
        "ensembl_protein_id": "",
        "strand": None,
        "mane_select": "",
        "consequence_terms": [],
        "alphamissense_score": None,
        "alphamissense_class": "",
        "conservation": None,
        "exon": "",
        "intron": "",
        "error": "",
        "protein_id": "",
        "protein_start": None,
        "protein_end": None,
        "amino_acids": "",
        "aa_ref": "",
    }


def parse_vep(vep_record):
    """Extract the fields we care about from a VEP record."""
    out = _empty_parsed()
    if not vep_record or "_error" in vep_record:
        out["error"] = vep_record.get("_error", "no data") if vep_record else "no data"
        return out

    out["most_severe_consequence"] = vep_record.get("most_severe_consequence", "")
    out["conservation"] = vep_record.get("conservation", None)

    tcs = vep_record.get("transcript_consequences", [])
    # `info` (Ensembl) carries gene/protein/AlphaMissense; `disp` (RefSeq when
    # available) carries the HGVS we want to display.
    info = _pick_info_transcript(tcs)
    disp = _pick_display_transcript(tcs)

    if info:
        out["gene_symbol"] = info.get("gene_symbol", "")
        out["gene_id"] = info.get("gene_id", "")
        # Ensembl transcript id (kept even when we display RefSeq HGVS) so the
        # gene-model diagram can look up exon coordinates from Ensembl.
        out["ensembl_transcript_id"] = info.get("transcript_id", "")
        # Ensembl protein id (ENSP) for the protein-domain diagram lookup. VEP
        # does not return a dedicated protein_id field, so fall back to the
        # accession embedded in the Ensembl hgvsp (e.g. "ENSP00000413251.3").
        out["ensembl_protein_id"] = (info.get("protein_id")
                                     or _accession_from_hgvs(info.get("hgvsp")))
        out["strand"] = info.get("strand")
        out["consequence_terms"] = info.get("consequence_terms", [])
        out["exon"] = info.get("exon", "")
        out["intron"] = info.get("intron", "")
        out["protein_start"] = info.get("protein_start")
        out["protein_end"] = info.get("protein_end")
        out["amino_acids"] = info.get("amino_acids", "")
        # amino_acids is like "T/A" (ref/alt); expose the reference residue
        aa = out["amino_acids"]
        out["aa_ref"] = aa.split("/")[0] if aa else ""
        # AlphaMissense fields (naming can vary by VEP version)
        am_score = info.get("alphamissense", {})
        if isinstance(am_score, dict):
            out["alphamissense_score"] = am_score.get("am_pathogenicity")
            out["alphamissense_class"] = am_score.get("am_class", "")
        else:
            out["alphamissense_score"] = info.get("am_pathogenicity")
            out["alphamissense_class"] = info.get("am_class", "")

    # Display transcript + HGVS: prefer RefSeq, fall back to Ensembl.
    if disp:
        out["transcript_id"] = disp.get("transcript_id", "")
        out["hgvsc"] = disp.get("hgvsc", "")
        out["hgvsp"] = disp.get("hgvsp", "")
        out["protein_id"] = disp.get("protein_id", "")
    # If the RefSeq display transcript lacked a coding HGVS, fall back to the
    # Ensembl one as a last resort so the field is never empty unnecessarily.
    if info and not out["hgvsc"]:
        out["transcript_id"] = info.get("transcript_id", out["transcript_id"])
        out["hgvsc"] = info.get("hgvsc", "")
        out["hgvsp"] = out["hgvsp"] or info.get("hgvsp", "")
        out["protein_id"] = out["protein_id"] or info.get("protein_id", "")

    # mane_select always holds the RefSeq NM accession (used by the MANE badge
    # and the 46-way conservation lookup).
    out["mane_select"] = _mane_refseq_accession(tcs)

    # hgvsg often appears at the top level or on the first transcript
    out["hgvsg"] = (vep_record.get("hgvsg", "")
                    or (disp.get("hgvsg", "") if disp else "")
                    or (info.get("hgvsg", "") if info else ""))

    return out
