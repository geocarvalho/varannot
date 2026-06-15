"""
conservation.py
===============
Amino-acid (protein) conservation track from the UCSC **multiz46way** alignment,
matching the per-species residue display used by GeneBe.

The output is a string like:

    TxTTTTTxxTTTTTTTTTTTTTTTTxTTTTTxVIITTTTIIVVIx

one character per non-human species (45 of them), in UCSC's fixed 46-way order,
giving the amino acid that species has at the codon aligned to the variant's
residue. Gaps / missing are shown as the gap character.

IMPORTANT ABOUT ASSEMBLIES
--------------------------
The multiz46way alignment exists ONLY for hg19. UCSC never produced a 46-way
for hg38. So for an hg38 input variant we lift the coordinate over to hg19
first (via pyliftover), then look up the residue in the hg19 exonAA alignment.

DATA SOURCE
-----------
UCSC pre-computed per-species protein alignments for CDS regions:
  https://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz46way/alignments/refGene.exonAA.fa.gz

This module downloads that file once (≈355 MB) into a cache directory and
builds a lightweight position index, so per-variant lookups are fast.

Because the file is large, the download is OFF by default. Enable the 46-way
track with --conservation46 (see annotate.py). Without it, the cell falls back
to the numeric phyloP/GERP score from VEP.
"""

import gzip
import os

# UCSC 46-way species in their canonical order (human first). Sourced from the
# multiz46way alignments README. The track we emit excludes human (the
# reference), giving 45 characters.
SPECIES_46WAY = [
    "hg19",      # Human (reference, excluded from emitted track)
    "panTro2",   # Chimp
    "gorGor1",   # Gorilla
    "ponAbe2",   # Orangutan
    "rheMac2",   # Rhesus
    "papHam1",   # Baboon
    "calJac1",   # Marmoset
    "tarSyr1",   # Tarsier
    "micMur1",   # Mouse lemur
    "otoGar1",   # Bushbaby
    "tupBel1",   # Tree shrew
    "mm9",       # Mouse
    "rn4",       # Rat
    "dipOrd1",   # Kangaroo rat
    "cavPor3",   # Guinea Pig
    "speTri1",   # Squirrel
    "oryCun2",   # Rabbit
    "ochPri2",   # Pika
    "vicPac1",   # Alpaca
    "turTru1",   # Dolphin
    "bosTau4",   # Cow
    "equCab2",   # Horse
    "felCat3",   # Cat
    "canFam2",   # Dog
    "myoLuc1",   # Microbat
    "pteVam1",   # Megabat
    "eriEur1",   # Hedgehog
    "sorAra1",   # Shrew
    "loxAfr3",   # Elephant
    "proCap1",   # Rock hyrax
    "echTel1",   # Tenrec
    "dasNov2",   # Armadillo
    "choHof1",   # Sloth
    "macEug1",   # Wallaby
    "monDom5",   # Opossum
    "ornAna1",   # Platypus
    "galGal3",   # Chicken
    "taeGut1",   # Zebra finch
    "anoCar1",   # Lizard
    "xenTro2",   # X. tropicalis
    "tetNig2",   # Tetraodon
    "fr2",       # Fugu
    "gasAcu1",   # Stickleback
    "oryLat2",   # Medaka
    "danRer6",   # Zebrafish
    "petMar1",   # Lamprey
]

EXONAA_URL = (
    "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/"
    "multiz46way/alignments/refGene.exonAA.fa.gz"
)


# ---------------------------------------------------------------------------
# Numeric conservation score (kept from VEP, e.g. GERP)
# ---------------------------------------------------------------------------
def conservation_score_from_vep(vep_parsed):
    """Pull a GERP/PhyloP-style numeric score if VEP returned one."""
    cons = vep_parsed.get("conservation")
    if isinstance(cons, (int, float)):
        return float(cons)
    return None


# ---------------------------------------------------------------------------
# 46-way amino-acid track
# ---------------------------------------------------------------------------
def _empty(note=""):
    return {
        "available": False, "track": "", "human_aa": "",
        "species": [], "n_species": 0, "note": note, "assembly": "hg19/46way",
    }


def build_protein_track(client, vep_parsed, enabled=False,
                        exonaa_path=None, liftover=None,
                        chrom=None, pos=None):
    """
    Build the 46-way per-species amino-acid track for a variant.

    Parameters
    ----------
    vep_parsed : dict  - needs mane_select (RefSeq NM id), protein_start, aa_ref
    enabled    : bool  - if False, returns a 'not run' marker (no big download)
    exonaa_path: str   - path to the downloaded refGene.exonAA.fa.gz

    The UCSC refGene.exonAA alignment is keyed by RefSeq accession and indexed
    by protein residue number, so the residue is located purely from the
    transcript's RefSeq id + protein position (no genomic coordinate or
    hg38->hg19 liftover is needed). `liftover`, `chrom`, and `pos` are accepted
    for backwards compatibility but are unused.

    Returns the standard conservation dict.
    """
    if not enabled:
        return _empty("46-way track not run (enable with --conservation46)")

    protein_pos = vep_parsed.get("protein_start")
    aa_ref = vep_parsed.get("aa_ref") or ""
    # The exonAA file uses RefSeq NM accessions; VEP's MANE Select field gives
    # exactly that (e.g. "NM_024665.7"). Strip the version for matching.
    refseq = (vep_parsed.get("mane_select") or "").split(".")[0]

    if not protein_pos:
        return _empty("no protein residue (non-coding variant)")
    if not refseq.startswith(("NM_", "NR_")):
        return _empty("no RefSeq (MANE) transcript for 46-way lookup")
    if not exonaa_path or not os.path.exists(exonaa_path):
        return _empty("46-way exonAA file not available; run the downloader")

    # Look up the residue column in the exonAA alignment for this transcript.
    block = _find_exon_block(exonaa_path, refseq, int(protein_pos))
    if block is None:
        return _empty(f"residue not found in 46-way exonAA alignment "
                      f"({refseq} p.{protein_pos})")

    human_seq = block.get("hg19", "")
    col = block["column"]
    if col is None or col >= len(human_seq):
        return _empty("residue column outside aligned block")

    human_aa = human_seq[col]
    track_chars = []
    species_names = []
    for sp in SPECIES_46WAY[1:]:  # skip human
        seq = block.get(sp, "")
        species_names.append(sp)
        if not seq or col >= len(seq):
            track_chars.append("-")
            continue
        aa = seq[col]
        track_chars.append(aa if aa not in (".", " ") else "-")

    out = _empty()
    out["available"] = True
    out["track"] = "".join(track_chars)
    out["human_aa"] = human_aa if human_aa not in ("-", ".", " ") else (aa_ref or "")
    out["species"] = species_names
    out["n_species"] = len(species_names)
    out["note"] = ""
    return out


def _norm_chrom(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else "chr" + chrom


def _parse_record_name(rec):
    """
    Parse a UCSC exonAA record name into (refseq, species, exon_index).

    The record name looks like ``NM_024665_hg19_1_14``:
      <refseq>_<species>_<exonIndex>_<totalExons>
    where the RefSeq accession itself contains an underscore (NM_024665), so we
    parse from the right: the last two tokens are the (1-based) exon index and
    the exon count; the token before them is the species; everything before
    that is the RefSeq accession.

    Returns (refseq, species, exon_index) or (None, None, None) if it doesn't
    fit the expected shape.
    """
    parts = rec.split("_")
    if len(parts) < 5:
        return None, None, None
    try:
        exon_index = int(parts[-2])
    except ValueError:
        return None, None, None
    species = parts[-3]
    refseq = "_".join(parts[:-3])
    return refseq, species, exon_index


def _find_exon_block(exonaa_path, refseq, protein_pos):
    """
    Scan the refGene exonAA FASTA for the records of `refseq`, find the exon
    that contains residue `protein_pos`, and compute the alignment column for
    that residue.

    The UCSC refGene.exonAA FASTA stores one record per (transcript, species,
    exon). Each header looks like:
      >NM_024665_hg19_1_14 19 0 1 chr3:176782708-176782765-
    followed by the amino-acid sequence for that exon (gaps as '-'). Records
    for one exon appear consecutively across species (hg19 first), and within
    an exon every species' sequence is aligned to the same length.

    We collect this transcript's per-exon, per-species sequences, walk the
    exons in order accumulating ungapped human residues, and stop at the exon
    spanning `protein_pos`.
    """
    target = refseq.split(".")[0]

    # exon_index -> {species: aligned_seq}
    exons = {}
    saw_target = False

    with gzip.open(exonaa_path, "rt") as fh:
        cur_species = None
        cur_exon = None
        in_target = False
        seq_parts = []

        def flush():
            if in_target and cur_exon is not None and cur_species is not None:
                exons.setdefault(cur_exon, {})[cur_species] = "".join(seq_parts)

        for line in fh:
            if line.startswith(">"):
                flush()
                rec = line[1:].split()[0]
                name, species, exon_index = _parse_record_name(rec)
                in_target = (name == target)
                if in_target:
                    saw_target = True
                elif saw_target:
                    # Records are grouped by transcript; once we've passed the
                    # target transcript's block we can stop scanning the file.
                    break
                cur_species = species
                cur_exon = exon_index
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        else:
            flush()  # only runs if the loop wasn't broken (target at EOF)

    if not exons:
        return None

    residues_before = 0
    for idx in sorted(exons):
        block = exons[idx]
        human = block.get("hg19", "")
        exon_res = len(human.replace("-", ""))
        if residues_before < protein_pos <= residues_before + exon_res:
            want = protein_pos - residues_before
            col = _ungapped_to_col(human, want)
            result = dict(block)
            result["column"] = col
            return result
        residues_before += exon_res

    return None


def _ungapped_to_col(aligned_seq, residue_n_1based):
    """Map the Nth ungapped residue (1-based) to its column in a gapped seq."""
    count = 0
    for col, ch in enumerate(aligned_seq):
        if ch != "-":
            count += 1
            if count == residue_n_1based:
                return col
    return None


# ---------------------------------------------------------------------------
# Downloader / setup helper
# ---------------------------------------------------------------------------
def ensure_exonaa(cache_dir, url=EXONAA_URL):
    """
    Ensure the refGene.exonAA.fa.gz file is present in cache_dir; download if
    not. Returns the local path. Intended to be called once on the user's
    machine. Requires `requests`.
    """
    import requests

    os.makedirs(cache_dir, exist_ok=True)
    dest = os.path.join(cache_dir, "refGene.exonAA.fa.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        return dest

    print(f"Downloading 46-way exonAA alignment (~355 MB) to {dest} ...")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = 0
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if total:
                    pct = 100 * done // total
                    print(f"\r  {pct:3d}%  ({done >> 20} / {total >> 20} MB)",
                          end="", flush=True)
    print("\n  done.")
    return dest


def make_liftover():
    """Create an hg38->hg19 LiftOver (downloads the chain file on first use)."""
    try:
        from pyliftover import LiftOver
        return LiftOver("hg38", "hg19")
    except Exception as exc:  # pragma: no cover
        print(f"  WARNING: liftover unavailable ({exc}); 46-way needs hg19 coords.")
        return None


# Backwards-compatible no-op kept so older imports don't break.
def build_match_track(client, chrom, pos, window=0):
    return _empty("nucleotide track deprecated; use the 46-way protein track")
