#!/usr/bin/env python3
"""
varannot web UI
===============
A small Flask app that wraps the annotation pipeline: paste variants (one per
line, ``chr,pos,ref,alt``) and optionally paste API/license keys, then get the
same HTML report the CLI produces.

Keys are handled **in memory only** for the duration of the request:
  * the form posts over POST (keys never land in the URL or server access logs),
  * key fields are password inputs and are never echoed back into the page,
  * keys are never written to disk or logged. The on-disk cache only stores
    public API responses (never the keys themselves).

Run it:
    python -m varannot.web                 # http://127.0.0.1:8000
    python -m varannot.web --port 8080 --host 0.0.0.0

(Port 8000 by default: on macOS port 5000 is used by the AirPlay Receiver and
returns HTTP 403, so we avoid it.)

Requires: flask (pip install flask). AlphaGenome predictions additionally need
`pip install alphagenome matplotlib` and a key.
"""

import argparse
import os
import threading
import time
import uuid

from flask import (Flask, Response, jsonify, redirect, render_template,
                   request, url_for)

from . import annotate
from .http_client import CachedSession
from .sources import alphagenome as ag_src
from .sources import conservation as cons_src
from .sources import omim as omim_src

app = Flask(__name__)

MAX_VARIANTS = 50

# In-memory job registry for background runs. Keyed by a random id. Note: this
# lives in the process memory, so run the dev server (or a single worker) — keys
# are kept here only while a job runs and are dropped as soon as it finishes.
_JOBS = {}
_JOBS_LOCK = threading.Lock()
_JOB_TTL = 2 * 60 * 60  # forget finished/abandoned jobs after 2 hours


def _job_update(job_id, **fields):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is not None:
            job.update(fields)


def _job_public(job):
    """The subset of a job that's safe to expose (no html, no keys)."""
    return {
        "state": job.get("state"),
        "message": job.get("message", ""),
        "current": job.get("current", 0),
        "total": job.get("total", 0),
        "error": job.get("error", ""),
    }


def _purge_jobs():
    now = time.time()
    with _JOBS_LOCK:
        stale = [k for k, j in _JOBS.items()
                 if now - j.get("created", now) > _JOB_TTL]
        for k in stale:
            _JOBS.pop(k, None)

# OMIM is read from the bundled download file (no key prompt). Override with the
# OMIM_GENEMAP2 env var if you keep it elsewhere.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENEMAP2_PATH = os.environ.get("OMIM_GENEMAP2",
                               os.path.join(REPO_ROOT, "data", "genemap2.txt"))

_OMIM_INDEX = None
_OMIM_LOADED = False


def _omim_index():
    """Load (once) the OMIM genemap2 index from data/genemap2.txt, if present."""
    global _OMIM_INDEX, _OMIM_LOADED
    if not _OMIM_LOADED:
        _OMIM_LOADED = True
        if os.path.exists(GENEMAP2_PATH):
            try:
                _OMIM_INDEX = omim_src.load_genemap2(GENEMAP2_PATH)
            except Exception:
                _OMIM_INDEX = None
    return _OMIM_INDEX

# Example placeholder shown in the textarea.
EXAMPLE = "chr3,177026399,C,T\nchrX,32389644,G,A\nchr17,43094692,C,T"


def _form_context(**extra):
    ctx = {
        "ag_installed": ag_src.is_available(),
        "max_variants": MAX_VARIANTS,
        "example": EXAMPLE,
        "examples": [ln for ln in EXAMPLE.splitlines() if ln.strip()],
        "error": None,
        "variants_text": "",
    }
    ctx.update(extra)
    return ctx


@app.get("/")
def index():
    return render_template("form.html.j2", **_form_context())


@app.post("/report")
def report():
    """Validate the form, start a background job, and show the running page."""
    variants_text = request.form.get("variants", "")
    variants = annotate.parse_variants_text(variants_text)

    if not variants:
        return render_template("form.html.j2", **_form_context(
            error="No valid variants found. Use one 'chr,pos,ref,alt' per line.",
            variants_text=variants_text)), 400

    if len(variants) > MAX_VARIANTS:
        return render_template("form.html.j2", **_form_context(
            error=f"{len(variants)} variants exceeds the limit of {MAX_VARIANTS}.",
            variants_text=variants_text)), 400

    # AlphaGenome key: held only inside the job dict while it runs, then dropped.
    ag_key = (request.form.get("alphagenome_key") or "").strip() or None

    try:
        top_n = int(request.form.get("alphagenome_top_tracks") or
                    ag_src.DEFAULT_TOP_TRACKS)
    except ValueError:
        top_n = ag_src.DEFAULT_TOP_TRACKS

    opts = {
        "spliceai": request.form.get("spliceai") == "on",
        "cons46": request.form.get("conservation46") == "on",
        "run_ag": request.form.get("run_alphagenome") == "on" and bool(ag_key),
        "ag_key": ag_key,
        "outputs": _split_csv(request.form.get("alphagenome_outputs")),
        "ontology": _split_csv(request.form.get("alphagenome_ontology")),
        "top_n": top_n,
    }

    _purge_jobs()
    job_id = uuid.uuid4().hex
    with _JOBS_LOCK:
        _JOBS[job_id] = {
            "state": "running", "message": "Starting...",
            "current": 0, "total": len(variants), "error": "",
            "html": None, "created": time.time(),
        }

    thread = threading.Thread(target=_run_job, args=(job_id, variants, opts),
                              daemon=True)
    thread.start()

    return redirect(url_for("running", job_id=job_id))


def _run_job(job_id, variants, opts):
    """Background worker: annotate all variants, optionally run AlphaGenome."""
    try:
        client = CachedSession()
        omim_index = _omim_index()
        ncbi_key = os.environ.get("NCBI_API_KEY")
        total = len(variants)

        exonaa_path = None
        if opts["cons46"]:
            _job_update(job_id, message="Preparing 46-way conservation "
                                        "(first run downloads ~355 MB)...")
            try:
                exonaa_path = cons_src.ensure_exonaa(client.cache_dir)
            except Exception:
                exonaa_path = None  # fall back to phyloP score only

        records = []
        for i, var in enumerate(variants, 1):
            _job_update(job_id, current=i, total=total,
                        message=(f"Annotating {var['chrom']}:{var['pos']} "
                                 f"{var['ref']}>{var['alt']} ({i}/{total})"))
            try:
                records.append(annotate.annotate_one(
                    client, var, None, ncbi_key,
                    spliceai_enabled=opts["spliceai"],
                    cons46_enabled=opts["cons46"], exonaa_path=exonaa_path,
                    omim_index=omim_index,
                ))
            except Exception as exc:  # keep going on per-variant failures
                records.append(annotate._error_record(
                    var, exc, spliceai_enabled=opts["spliceai"]))

        alphagenome = None
        if opts["run_ag"] and ag_src.is_available():
            _job_update(job_id, current=0, total=total,
                        message="Running AlphaGenome predictions "
                                "(this can take a while per variant)...")

            def _ag_progress(done, tot, label):
                _job_update(job_id, current=done, total=tot,
                            message=f"AlphaGenome: {label} ({done}/{tot})")

            meta = [{
                "chrom": r["input"]["chrom"], "pos": r["input"]["pos"],
                "ref": r["input"]["ref"], "alt": r["input"]["alt"],
                "label": (f"{r['input']['chrom']}:{r['input']['pos']} "
                          f"{r['input']['ref']}>{r['input']['alt']}"),
                "gene": r.get("vep", {}).get("gene_symbol", ""),
            } for r in records]
            alphagenome = ag_src.run_for_variants(
                meta, opts["ag_key"], cache_dir=client.cache_dir,
                output_types=[o.upper() for o in opts["outputs"]] or None,
                ontology_terms=opts["ontology"] or None, top_n=opts["top_n"],
                progress=_ag_progress)

        _job_update(job_id, message="Rendering report...")
        html = annotate.render_report_str(records, alphagenome=alphagenome,
                                          web=True)
        _job_update(job_id, state="done", message="Done", html=html)
    except Exception as exc:
        _job_update(job_id, state="error", error=str(exc))
    finally:
        # Drop the key as soon as the work is finished.
        opts.pop("ag_key", None)


@app.get("/running/<job_id>")
def running(job_id):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
    if job is None:
        return redirect(url_for("index"))
    if job["state"] == "done":
        return redirect(url_for("result", job_id=job_id))
    return render_template("running.html.j2", job_id=job_id,
                           initial=_job_public(job))


@app.get("/progress/<job_id>")
def progress(job_id):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
    if job is None:
        return jsonify({"state": "missing"}), 404
    return jsonify(_job_public(job))


@app.get("/result/<job_id>")
def result(job_id):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
    if job is None:
        return redirect(url_for("index"))
    if job["state"] != "done":
        return redirect(url_for("running", job_id=job_id))
    html = job.get("html") or ""
    # Hand the report over once, then free the memory it held.
    with _JOBS_LOCK:
        _JOBS.pop(job_id, None)
    return Response(html, mimetype="text/html")


def _split_csv(value):
    return [v.strip() for v in (value or "").split(",") if v.strip()]


def main():
    parser = argparse.ArgumentParser(description="varannot web UI")
    parser.add_argument("--host", default="127.0.0.1",
                        help="Interface to bind (default 127.0.0.1, localhost only)")
    parser.add_argument("--port", type=int, default=8000,
                        help="Port (default 8000; avoid 5000, used by macOS AirPlay)")
    args = parser.parse_args()
    # debug=False on purpose: avoids the reloader/debugger echoing request data.
    # threaded=True so the background job runs while the browser polls progress.
    app.run(host=args.host, port=args.port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
