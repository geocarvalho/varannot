# varannot

A small command-line framework that annotates a list of genomic variants
(GRCh38) and produces a **single combined HTML report**. Built for small
batches (up to 50 variants) in rare-disease analysis.

Each variant gets:

- **HGVS nomenclature** (g. / c. / p.) via Ensembl VEP, preferring **RefSeq**
  transcripts (MANE Select first, then RefSeq canonical) and falling back to
  Ensembl only as a last resort
- **Gene-model diagram**: an inline SVG showing where the variant lands in the
  exon/intron structure (or how far it is from the nearest exon)
- **Protein-domain diagram**: an inline SVG placing the variant on the protein
  relative to Pfam domains (with the truncated region shaded for stop-gain /
  frameshift variants)
- **Conservation**: a GERP/PhyloP score plus the **UCSC multiz46way per-species
  amino-acid track** (`TxTTTTTxx…`) — the residue each of 45 species has at the
  codon aligned to the variant, in UCSC's fixed 46-way order, matching the
  GeneBe display. Opt-in (`--conservation46`); see notes below.
- **gnomAD v4**: homozygote and heterozygote counts, allele frequency,
  and a direct browser link
- **AlphaMissense**: pathogenicity score and class (for missense variants)
- **REVEL**: ensemble missense pathogenicity score (plus CADD phred as a bonus),
  via MyVariant.info / dbNSFP
- **SpliceAI**: splice-disruption delta scores (acceptor/donor gain/loss) and
  the overall max delta, via the Broad SpliceAI-lookup API (opt-in)
- **ClinVar**: classification, review status, submission count, link
- **Gene constraint** (gnomAD): pLI, missense Z, LOEUF, pLoF Z
- **GTEx**: median gene expression (TPM) across the top tissues, always
  including muscle, fibroblast, and whole blood
- **OMIM**: gene MIM number and associated phenotypes with inheritance

**AlphaGenome** (Google DeepMind) variant-effect plots are available either as a
**separate report** (its own command) or folded into the main report as a second
**tab** — see [AlphaGenome report](#alphagenome-report). They're opt-in because
the predictions are slower and need an API key.

Everything runs from your command line and calls public APIs. Responses are
cached on disk, so re-running on the same variants is instant and offline.

> **Note on "offline":** the underlying datasets (gnomAD, ClinVar, OMIM,
> conservation) are far too large to bundle locally, so the first run needs
> network access to the public APIs. After that, the on-disk cache lets you
> re-render reports without hitting the network again.

## Install

```bash
pip install -r requirements.txt
```

The AlphaGenome report additionally needs `alphagenome` and `matplotlib`
(optional; only required for the separate AlphaGenome command):

```bash
pip install alphagenome matplotlib
```

## Data sources & keys

| Source        | Endpoint                              | Key needed?            |
|---------------|---------------------------------------|------------------------|
| Ensembl VEP   | rest.ensembl.org                      | No                     |
| Gene model / Pfam | rest.ensembl.org (lookup / overlap) | No                   |
| gnomAD v4     | gnomad.broadinstitute.org/api         | No                     |
| AlphaMissense | via Ensembl VEP                       | No                     |
| REVEL + CADD  | myvariant.info (dbNSFP)               | No                     |
| SpliceAI      | spliceailookup-api (Broad, Cloud Run) | No (rate limited)      |
| ClinVar       | eutils.ncbi.nlm.nih.gov               | Optional (NCBI key)    |
| GTEx v8       | gtexportal.org/api/v2                  | No                     |
| OMIM          | local genemap2.txt (download files)   | No (or licensed key)   |
| Conservation  | UCSC multiz46way exonAA + VEP          | No (355 MB download)   |
| AlphaGenome   | AlphaGenome API (DeepMind)             | Yes (API key)          |

OMIM offers two kinds of access — **download files** or an **API key** — and they
are granted separately. This tool prefers the **download files**: if you have
`genemap2.txt`, point `--omim-genemap2` at it. No key, no rate limits, offline.

```bash
python -m varannot.annotate --input test/example_variants.csv --output report.html \
    --omim-genemap2 /path/to/genemap2.txt
```

If you instead have a licensed API key, you can still use `--omim-key` (or the
`OMIM_API_KEY` env var). Request either at <https://omim.org/api>.
An NCBI key (optional, raises ClinVar rate limits) is at
<https://www.ncbi.nlm.nih.gov/account/>.

## Usage

Input is a CSV/TSV with `chr,pos,ref,alt` (a header row is optional):

```
chr3,177026399,C,T
chrX,32389644,G,A
```

Run:

```bash
# Using the OMIM download file (preferred)
python -m varannot.annotate \
    --input test/example_variants.csv \
    --output report.html \
    --omim-genemap2 /path/to/genemap2.txt

# Or via environment variables
export OMIM_GENEMAP2=/path/to/genemap2.txt
export NCBI_API_KEY=YOUR_NCBI_KEY      # optional
python -m varannot.annotate --input test/example_variants.csv --output report.html
```

Open `report.html` in any browser.

## Options

| Flag             | Default            | Description                                  |
|------------------|--------------------|----------------------------------------------|
| `--input`        | (required)         | CSV/TSV of variants                          |
| `--output`       | variant_report.html| Output HTML file                             |
| `--omim-genemap2`| `$OMIM_GENEMAP2`   | Path to OMIM genemap2.txt (preferred)        |
| `--omim-key`     | `$OMIM_API_KEY`    | OMIM API key (alternative to genemap2)       |
| `--ncbi-key`     | `$NCBI_API_KEY`    | Optional NCBI E-utilities key                |
| `--cache-dir`    | `.varannot_cache`  | Where cached API responses live              |
| `--min-interval` | 0.4                | Seconds between live API calls (politeness)  |
| `--spliceai`     | off                | Enable SpliceAI scores (opt-in; rate limited)|
| `--spliceai-url` | (Broad endpoint)   | Custom SpliceAI API base URL (own instance)  |
| `--conservation46`| off               | Enable the UCSC 46-way amino-acid track      |
| `--exonaa-dir`   | `.varannot_cache`  | Where the UCSC exonAA file is downloaded      |
| `--max-variants` | 50                 | Safety cap on batch size                     |

### SpliceAI notes

SpliceAI is **opt-in** because the public Broad endpoint only allows a few
requests per minute. When enabled, the tool waits ~6 s between SpliceAI calls,
so a 50-variant batch takes about 5 minutes for that section alone (everything
else still runs at full speed). Scores are cached like all other sources.

```bash
# Enable SpliceAI using the public endpoint
python -m varannot.annotate --input test/example_variants.csv --output report.html \
    --omim-genemap2 genemap2.txt --spliceai

# Faster: point at your own SpliceAI-lookup instance (no rate limit)
python -m varannot.annotate --input test/example_variants.csv --output report.html \
    --omim-genemap2 genemap2.txt --spliceai --spliceai-url http://localhost:8080
```

The report shows the **max delta score** (0–1) with the standard interpretation
bands (≥0.2 possible, ≥0.5 likely, ≥0.8 high-precision) plus per-transcript
acceptor/donor gain/loss values for any transcript scoring ≥0.2.

To force fresh queries, delete the cache directory:

```bash
rm -rf .varannot_cache
```

## AlphaGenome report

AlphaGenome (Google DeepMind) variant-effect predictions are opt-in (they're
slower and need an API key). You can get them two ways:

1. **Separate report** (default, keeps the main report fast) — its own command
   and HTML file.
2. **Combined report** — add `--alphagenome-key` to the normal `annotate`
   command and the report gains a second **AlphaGenome tab**.

Requirements for both: `pip install alphagenome matplotlib` and an AlphaGenome
API key (pass `--alphagenome-key` or set `ALPHAGENOME_API_KEY`).

```bash
# 1) Separate AlphaGenome-only report
python -m varannot.alphagenome_report \
    --input test/example_variants.csv \
    --output alphagenome_report.html \
    --alphagenome-key YOUR_KEY

# 2) Combined report: annotations + an AlphaGenome tab in one file
python -m varannot.annotate \
    --input test/example_variants.csv \
    --output report.html \
    --omim-genemap2 genemap2.txt \
    --alphagenome-key YOUR_KEY
```

The same `--alphagenome-*` options (`--alphagenome-outputs`,
`--alphagenome-ontology`, `--alphagenome-length`, `--alphagenome-top-tracks`)
work on both commands.

The report has a variant selector and a grid of plots per variant. By default it
requests all available output types; restrict them for speed:

```bash
python -m varannot.alphagenome_report --input test/example_variants.csv \
    --alphagenome-key YOUR_KEY \
    --alphagenome-outputs RNA_SEQ,DNASE,ATAC,SPLICE_JUNCTIONS
```

| Flag                     | Default                  | Description                                   |
|--------------------------|--------------------------|-----------------------------------------------|
| `--input`                | (required)               | CSV/TSV of variants                           |
| `--output`               | alphagenome_report.html  | Output HTML file                              |
| `--alphagenome-key`      | `$ALPHAGENOME_API_KEY`   | AlphaGenome API key (required)                |
| `--alphagenome-outputs`  | all                      | Comma-separated output types                  |
| `--alphagenome-ontology` | all                      | Comma-separated ontology CURIEs to restrict tracks |
| `--alphagenome-top-tracks` | 1                      | Plot only the N most variant-relevant tracks per output |
| `--alphagenome-length`   | SEQUENCE_LENGTH_100KB    | Input window (16KB / 100KB / 500KB / 1MB)     |
| `--cache-dir`            | `.varannot_cache`        | Where generated plot PNGs are cached          |
| `--max-variants`         | 50                       | Safety cap on batch size                      |

Available output types: `ATAC`, `CAGE`, `DNASE`, `RNA_SEQ`, `CHIP_HISTONE`,
`CHIP_TF`, `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`,
`CONTACT_MAPS`, `PROCAP`.

Without an ontology filter, AlphaGenome returns hundreds of tracks per output
(more than its plotter allows). By default the report therefore shows only the
**single track/tissue where the variant has the largest predicted effect** for
each output, and labels the plot with that track's name. Raise
`--alphagenome-top-tracks` to show more (e.g. `--alphagenome-top-tracks 5`), or
restrict tissues explicitly with `--alphagenome-ontology`.

## Project layout

```
varannot/
├── annotate.py            # main CLI orchestrator + HTML rendering
├── alphagenome_report.py  # separate CLI for the AlphaGenome report
├── http_client.py         # cached, rate-limited HTTP session
├── sources/
│   ├── vep.py             # HGVS (RefSeq/MANE), consequence, AlphaMissense
│   ├── gene_model.py      # exon/intron SVG diagram (Ensembl lookup)
│   ├── protein_domains.py # protein/Pfam SVG diagram (Ensembl overlap)
│   ├── gnomad.py          # variant counts + gene constraint (GraphQL)
│   ├── clinvar.py         # NCBI E-utilities lookup
│   ├── omim.py            # OMIM API (licensed key)
│   ├── revel.py           # REVEL + CADD (MyVariant.info / dbNSFP)
│   ├── spliceai.py        # SpliceAI delta scores (Broad lookup API)
│   ├── gtex.py            # GTEx v8 median gene expression
│   ├── conservation.py    # UCSC multiz species match track
│   └── alphagenome.py     # AlphaGenome predictions + matplotlib plots
└── templates/
    ├── report.html.j2              # main report (gains tabs when combined)
    ├── alphagenome_report.html.j2  # standalone AlphaGenome report
    ├── _ag_styles.html.j2          # shared AlphaGenome CSS (grid + lightbox)
    ├── _ag_panel.html.j2           # shared AlphaGenome plot grid + selector
    └── _ag_lightbox.html.j2        # shared click-to-zoom lightbox + JS
```

## Extending

Each source module exposes a single `query_*` function that takes the shared
`CachedSession` and returns a plain dict. To add a new annotation (e.g. SpliceAI,
REVEL), add a module under `sources/`, call it from `annotate_one()` in
`annotate.py`, and add a cell to `templates/report.html.j2`.

## Caveats

- ClinVar matching is by genomic position + allele; for complex indels confirm
  the matched record's title in the report.
### Conservation (46-way) notes

The amino-acid conservation track reproduces the UCSC **multiz46way** display
that GeneBe shows: for the variant's codon, the amino acid each of 45 vertebrate
species carries, in UCSC's fixed 46-way species order.

Two things to know:

1. **It's opt-in** (`--conservation46`) because it downloads UCSC's
   `refGene.exonAA.fa.gz` (~355 MB) once into `--exonaa-dir`. After that,
   lookups are local and fast.
2. **No coordinate liftover is needed.** The exonAA file is indexed by RefSeq
   accession (e.g. `NM_024665`) and protein residue, so the track is looked up
   directly from VEP's MANE/RefSeq transcript and the variant's protein
   position — no hg38→hg19 conversion involved.

```bash
python -m varannot.annotate --input test/example_variants.csv --output report.html \
    --omim-genemap2 genemap2.txt --conservation46 --exonaa-dir ~/varannot_data
```

Only coding (missense) variants get a track; for others the cell shows the
numeric phyloP/GERP score instead. The track should match GeneBe's character
output on standard codons; rare edge cases (alignment gaps at the exact codon,
unusual transcripts) may differ slightly.
- gnomAD v4 constraint reports LOEUF (preferred) rather than the legacy pRec.
- Annotations are for research use and always require expert review.
