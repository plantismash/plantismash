# License: GNU Affero General Public License v3 or later

from typing import Any, Dict, List, Tuple
from pyquery import PyQuery as pq
import logging
import os
import csv
import io
import urllib.parse
import gzip
import posixpath

from Bio.SeqFeature import CompoundLocation
from antismash import utils

# Reuse your dataclasses/types from the detector
from .tfbs_detection import TFBSFinderResults, TFBSHit

# Validation/annotation helpers (TAIR & connectTF)
from .validation import load_tair_symbols, load_validated_map, evidence_summary, norm_agi

CONNECTF_URL = "https://connectf.org/"
GLOBAL_CSV_NAME = "tfbs_hits_all_bgcs.csv"


def _safe_filename(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in str(s))


def _get_output_dir(options) -> str:
    """
    Return the plantiSMASH run's results directory.
    Prefer the absolute path set by run_antismash.py:
      - options.full_outputfolder_path
      - else abspath(options.outputfoldername)
    """
    val = getattr(options, "full_outputfolder_path", None)
    if isinstance(val, str) and val.strip():
        return val

    val = getattr(options, "outputfoldername", None)
    if isinstance(val, str) and val.strip():
        return os.path.abspath(val)

    # Fallbacks
    candidates = [
        "outputfolder", "output_dir", "outdir", "output", "output_folder",
        "output_path", "results_dir", "results_path", "result_dir", "work_dir",
    ]
    for attr in candidates:
        v = getattr(options, attr, None)
        if isinstance(v, str) and v.strip():
            return os.path.abspath(v)

    general = getattr(options, "general", None)
    if general:
        for attr in ["full_outputfolder_path", "outputfoldername"] + candidates:
            v = getattr(general, attr, None)
            if isinstance(v, str) and v.strip():
                return os.path.abspath(v)

    return os.getcwd()


def _cds_tss_and_strand(cds) -> Tuple[int, int]:
    """Return TSS (on forward axis) and strand for a CDS (handles split genes)."""
    loc = getattr(cds, "location", None)
    if loc is None:
        return None, None
    strand = int(getattr(loc, "strand", 1) or 1)
    if isinstance(loc, CompoundLocation) and loc.parts:
        first, last = loc.parts[0], loc.parts[-1]
        tss = int(first.start) if strand == 1 else int(last.end) - 1
    else:
        tss = int(loc.start) if strand == 1 else int(loc.end) - 1
    return tss, strand


def _cds_label(cds, idx: int) -> str:
    q = getattr(cds, "qualifiers", {}) or {}
    for key in ("locus_tag", "gene", "protein_id"):
        vals = q.get(key)
        if vals:
            return str(vals[0])
    return f"CDS_{idx+1}"


def _signed_distance(hit_start: int, tss: int, gene_strand: int) -> int:
    """Negative = upstream relative to gene; Positive = downstream."""
    return (hit_start - tss) if gene_strand == 1 else (tss - hit_start)


def _load_hits_from_options(seq_record, options) -> List[TFBSHit]:
    """Pull hits from options.extrarecord (populated by run_tfbs_finder_for_record)."""
    try:
        ns = options.extrarecord.get(seq_record.id)
        if not ns or "TFBSFinderResults" not in getattr(ns, "extradata", {}):
            return []
        results: TFBSFinderResults = ns.extradata["TFBSFinderResults"]
        hits_json = [h.to_json() for h in results.hits_by_record.get(seq_record.id, [])]
        return [TFBSHit.from_json(h) for h in hits_json]
    except Exception as e:
        logging.warning("TFBS panel: could not load hits from options: %s", e)
        return []


def _csv_download_link(headers: List[str], rows: List[List[str]], filename: str) -> str:
    """Create a small CSV in-memory and return a HTML anchor with data URI for download."""
    try:
        buf = io.StringIO()
        w = csv.writer(buf, delimiter=",", lineterminator="\n")
        w.writerow(headers)
        for r in rows:
            r = [(x.replace("−", "-").replace("＋", "+") if isinstance(x, str) else x) for x in r]
            w.writerow(r)
        csv_text = buf.getvalue()
        data_uri = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_text)
        return f' <a download="{_safe_filename(filename)}" href="{data_uri}">Download CSV</a>'
    except Exception:
        return ""


def _append_global_csv(outdir: str, headers: List[str], rows: List[List[str]]) -> None:
    """Append rows to the single run-level CSV at outdir/GLOBAL_CSV_NAME; write header once."""
    try:
        os.makedirs(outdir, exist_ok=True)
        fpath = os.path.join(outdir, GLOBAL_CSV_NAME)
        write_header = not os.path.exists(fpath)
        with open(fpath, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter=",", lineterminator="\n")
            if write_header:
                w.writerow(headers)
            for r in rows:
                r = [(x.replace("−", "-").replace("＋", "+") if isinstance(x, str) else x) for x in r]
                w.writerow(r)
        if write_header:
            logging.info("TFBS: created %s", fpath)
    except Exception as e:
        logging.warning("TFBS: failed to write global CSV: %s", e)


def _ensure_tfbs_dir(options) -> Tuple[str, str]:
    """Returns (abs_dir, rel_href_prefix) for a 'tfbs' subfolder under the results directory."""
    outdir = _get_output_dir(options)
    tfbs_dir = os.path.join(outdir, "tfbs")
    os.makedirs(tfbs_dir, exist_ok=True)
    return tfbs_dir, "tfbs"  # POSIX for href


def _human_size(n: int) -> str:
    units = ["B", "KB", "MB", "GB"]
    size = float(n)
    for u in units:
        if size < 1024 or u == units[-1]:
            return f"{size:.1f} {u}"
        size /= 1024.0
    return f"{size:.1f} GB"


def _count_csv_rows(csv_path: str) -> int:
    try:
        with open(csv_path, "r", encoding="utf-8") as fh:
            return max(0, sum(1 for _ in fh) - 1)  # minus header
    except Exception:
        return 0


def _write_tsv_gz_from_rows(
    headers: List[str],
    rows: List[List[str]],
    abs_dir: str,
    rel_prefix: str,
    filename: str,
) -> Tuple[str, int, int]:
    """
    Write rows as a gzipped TSV at abs_dir/filename and return:
      (rel_href, n_rows, file_size_bytes)
    """
    abs_path = os.path.join(abs_dir, filename)
    with gzip.open(abs_path, "wt", newline="", encoding="utf-8") as gz:
        w = csv.writer(gz, delimiter="\t", lineterminator="\n")
        w.writerow(headers)
        for r in rows:
            w.writerow(r)
    try:
        size = os.path.getsize(abs_path)
    except Exception:
        size = 0
    rel_href = posixpath.join(rel_prefix, filename)
    return rel_href, len(rows), size


def _build_tables_for_cluster(seq_record, cluster, options) -> Tuple[str, str]:
    """
    Returns (per_cds_html, per_hit_download_html).
    - per_cds_html: aggregated motif list per CDS (inline; with small CSV download)
    - per_hit_download_html: no giant table; only a Download link to TSV.gz on disk
    Also appends all detailed hits to a single run-level CSV in the results directory.
    """
    hits = _load_hits_from_options(seq_record, options)
    if not hits:
        return "", ""

    # Restrict to CDS overlapping this cluster
    cluster_feat = utils.get_cluster_by_nr(seq_record, cluster["idx"])
    cstart, cend = int(cluster_feat.location.start), int(cluster_feat.location.end)
    cds_feats = [
        f for f in utils.get_cds_features(seq_record)
        if int(f.location.end) > cstart and int(f.location.start) < cend
    ]

    half = int(getattr(options, "tfbs_range", 1000))
    seqlen = len(seq_record.seq)

    # Build CDS windows (and clip display to cluster span)
    cds_windows: List[Dict[str, Any]] = []
    for i, cds in enumerate(cds_feats):
        tss, gstrand = _cds_tss_and_strand(cds)
        if tss is None:
            continue
        label = _cds_label(cds, i)
        a = max(0, tss - half)
        b = min(seqlen, tss + half + 1)  # half-open
        disp_a = max(a, cstart)
        disp_b = min(b, cend)
        cds_windows.append(
            {
                "label": label,
                "strand": gstrand,
                "tss": tss,
                "win_a": a,
                "win_b": b,
                "disp_a": disp_a,
                "disp_b": disp_b,
            }
        )
    if not cds_windows:
        return "", ""

    # Mappings
    tair_map = load_tair_symbols()
    validated_map = load_validated_map()

    # Accumulators
    per_hit_rows_csv: List[List[str]] = []   # detailed rows for per-BGC TSV.gz + global CSV
    per_cds_rows_csv: List[List[str]] = []   # aggregated rows for small in-memory CSV
    per_cds: Dict[str, Dict[str, Any]] = {}
    validated_hits_total = 0

    outdir = _get_output_dir(options)
    global_headers = [
        "Record", "Cluster", "ClusterStart", "ClusterEnd", "Product",
        "CDS", "CDS_strand", "TSS", "Window",
        "TF_locus", "TF_symbol", "Hit_start", "Hit_strand",
        "Distance_to_TSS", "Score", "Confidence", "Validated", "Evidence", "TF_link",
    ]
    global_rows: List[List[str]] = []
    product = cluster.get("type", cluster.get("product", "-"))

    # Build rows
    for cw in cds_windows:
        label, gstrand, tss, a, b, disp_a, disp_b = (
            cw["label"], cw["strand"], cw["tss"], cw["win_a"], cw["win_b"], cw["disp_a"], cw["disp_b"]
        )
        target_norm = norm_agi(label)

        for h in hits:
            mlen = len(h.consensus) if h.consensus else 1
            if a <= h.start <= max(a, b - mlen):
                dist = _signed_distance(h.start, tss, gstrand)

                tf_id_norm = norm_agi(h.name or "")
                tf_sym = tair_map.get(tf_id_norm, "")

                meta_id = validated_map.get((tf_id_norm, target_norm))
                is_validated = meta_id is not None
                ev_txt = evidence_summary(meta_id) if is_validated else ""
                if is_validated:
                    validated_hits_total += 1

                # ASCII for CSV
                strand_gene_csv = "+" if gstrand == 1 else "-"
                strand_hit_csv = "+" if h.strand == 1 else "-"

                window_str = f"{disp_a}-{disp_b}"

                # Per-BGC TSV.gz row
                per_hit_rows_csv.append(
                    [
                        label,
                        strand_gene_csv,
                        str(tss),
                        window_str,
                        tf_id_norm,
                        tf_sym,
                        str(h.start),
                        strand_hit_csv,
                        str(dist),
                        f"{h.score:.2f}/{h.max_score:.2f}",
                        str(h.confidence).capitalize(),
                        ev_txt,
                    ]
                )

                # Global CSV row (run-level)
                global_rows.append(
                    [
                        seq_record.id,
                        str(cluster["idx"]),
                        str(cstart),
                        str(cend),
                        str(product),
                        label,
                        strand_gene_csv,
                        str(tss),
                        window_str,
                        tf_id_norm,
                        tf_sym,
                        str(h.start),
                        strand_hit_csv,
                        str(dist),
                        f"{h.score:.2f}/{h.max_score:.2f}",
                        str(h.confidence).capitalize(),
                        "yes" if is_validated else "no",
                        ev_txt,
                        getattr(h, "link", "") or "",
                    ]
                )

                # Aggregate per CDS
                agg = per_cds.setdefault(
                    label,
                    {
                        "strand": gstrand,
                        "tss": tss,
                        "disp_a": disp_a,
                        "disp_b": disp_b,
                        "motifs": {},        # TF id -> link
                        "validated": set(),  # TF ids validated for this CDS
                        "count": 0,
                    },
                )
                if tf_id_norm not in agg["motifs"] or not agg["motifs"][tf_id_norm]:
                    agg["motifs"][tf_id_norm] = getattr(h, "link", "") or ""
                if is_validated:
                    agg["validated"].add(tf_id_norm)
                agg["count"] += 1

    # Totals for this cluster
    total_hits = len(per_hit_rows_csv)

    # Append to global CSV
    if global_rows:
        _append_global_csv(outdir, global_headers, global_rows)

    # Aggregated per-CDS table (inline)
    per_cds_html = ""
    if per_cds:
        rows_html: List[str] = []
        per_cds_rows_csv = []
        for label, agg in sorted(per_cds.items()):
            # HTML motifs with TAIR symbol and '*' if validated
            motifs_html = ", ".join(
                (
                    (lambda nm, lnk: (
                        f'<a href="{lnk}" target="_blank" rel="noopener">'
                        f'{nm}{" (" + tair_map.get(nm, "") + ")" if tair_map.get(nm, "") else ""}'
                        f'{" *" if nm in agg["validated"] else ""}'
                        f"</a>"
                    ) if lnk else
                    f'{nm}{" (" + tair_map.get(nm, "") + ")" if tair_map.get(nm, "") else ""}'
                    f'{" *" if nm in agg["validated"] else ""}')
                )(nm, lnk)
                for nm, lnk in sorted(agg["motifs"].items())
            )
            window_str_cds = f"{agg['disp_a']}-{agg['disp_b']}"
            strand_html = '+' if agg['strand'] == 1 else '−'

            rows_html.append(
                "<tr>"
                f"<td>{label}</td>"
                f"<td>{strand_html}</td>"
                f"<td>{agg['tss']}</td>"
                f"<td>{window_str_cds}</td>"
                f"<td>{motifs_html}</td>"
                f"<td>{agg['count']}</td>"
                "</tr>"
            )

            # CSV (plain text motifs; symbol + '*' if validated)
            motifs_txt = ", ".join(
                f"{nm}"
                f"{' (' + tair_map.get(nm, '') + ')' if tair_map.get(nm, '') else ''}"
                f"{' *' if nm in agg['validated'] else ''}"
                for nm, _ in sorted(agg["motifs"].items())
            )
            per_cds_rows_csv.append(
                [
                    label,
                    ("+" if agg["strand"] == 1 else "-"),
                    str(agg["tss"]),
                    window_str_cds,
                    motifs_txt,
                    str(agg["count"]),
                ]
            )

        dl_cds = _csv_download_link(
            headers=["CDS", "CDS strand", "TSS", "Window", "Motifs", "#Hits"],
            rows=per_cds_rows_csv,
            filename=f"tfbs_motifs_per_cds_cluster{cluster['idx']}_{_safe_filename(seq_record.id)}.csv",
        )

        per_cds_html = (
            "<details class='tfbs-block'>"
            "<summary><strong>Motifs per CDS (aggregated)</strong> "
            f"<span style='opacity:.7'>(validated: {validated_hits_total:,} / total: {total_hits:,})</span> "
            "<span style='opacity:.7'>* validated via connectTF.org</span>"
            f"{dl_cds}</summary>"
            "<div class='mt-2'>"
            "<table class='table table-sm'>"
            "<thead><tr>"
            "<th>CDS</th><th>CDS strand</th><th>TSS</th><th>Window</th>"
            "<th>Motifs</th><th>#Hits</th>"
            "</tr></thead><tbody>"
            + "".join(rows_html)
            + "</tbody></table>"
            "</div>"
            "</details>"
        )

    # Per-BGC detailed hits: write TSV.gz and show a download button
    per_hit_html = ""
    if per_hit_rows_csv:
        tfbs_dir, rel_prefix = _ensure_tfbs_dir(options)
        file_base = f"tfbs_hits_cluster{cluster['idx']}_{_safe_filename(seq_record.id)}.tsv.gz"
        dl_headers = [
            "CDS", "CDS_strand", "TSS", "Window", "TF_locus", "TF_symbol",
            "Hit_start", "Hit_strand", "Distance_to_TSS", "Score", "Confidence", "Evidence"
        ]
        rel_href, nrows, fsize = _write_tsv_gz_from_rows(
            dl_headers, per_hit_rows_csv, tfbs_dir, rel_prefix, file_base
        )
        per_hit_html = (
            "<p class='tfbs-download'>"
            f"<a class='as-button' href='{rel_href}' download>"
            f"Download full TFBS hits for this BGC (TSV.gz, {nrows:,} rows, {_human_size(fsize)})"
            "</a> "
            f"<span style='opacity:.7'>(validated: {validated_hits_total:,} / total: {total_hits:,})</span>"
            "</p>"
        )

    return per_cds_html, per_hit_html


def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    """Insert the TFBS panel into the cluster page and add a run-level CSV download if present."""
    try:
        per_cds_html, per_hit_html = _build_tables_for_cluster(seq_record, cluster, options)
    except Exception as e:
        logging.exception("TFBS panel: failed to render: %s", e)
        return details

    if not per_cds_html and not per_hit_html:
        return details  # keep the page clean if nothing to show

    if details is None:
        details = pq("<div>")
        details.addClass("details")

    header = pq("<h3>")
    header.text(f"Transcription factor binding sites (±{getattr(options, 'tfbs_range', 1000)} bp)")
    details.append(header)

    container = pq("<div>")
    container.addClass("tfbs-container")
    if per_cds_html:
        container.append(pq(per_cds_html))
    if per_hit_html:
        container.append(pq(per_hit_html))

    # Run-level CSV link (if available)
    outdir = _get_output_dir(options)
    global_csv = os.path.join(outdir, GLOBAL_CSV_NAME)
    if os.path.exists(global_csv):
        nrows = _count_csv_rows(global_csv)
        size = _human_size(os.path.getsize(global_csv))
        link_html = (
            "<p class='tfbs-download-all'>"
            f"<a class='as-button' href='{GLOBAL_CSV_NAME}' download>"
            f"Download TFBS hits for all BGCs (CSV, {nrows:,} rows, {size})"
            "</a></p>"
        )
        container.append(pq(link_html))

    details.append(container)
    return details