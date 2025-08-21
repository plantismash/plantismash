# License: GNU AGPL v3 or later

"""
TFBS detection (per-BGC) using PWMs with MOODS.

This implementation:
- iterates over each BGC (cluster) on the record,
- builds ±range promoter windows for CDS that overlap the BGC,
- clips each window to the BGC span,
- merges overlapping windows within that BGC,
- scans each merged interval exactly once,
- aggregates hits across all BGCs (the HTML/output module maps hits to clusters).

Public entry point: run_tfbs_finder(record, pvalue, start_overlap, matrix_path=PWM_PATH)
"""

from __future__ import annotations
from dataclasses import dataclass
from enum import IntEnum, auto
from typing import Any, Dict, List, Optional, Tuple
import os
import json
import logging
import tempfile

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation

from MOODS import tools as moods_tools
from MOODS import scan, parsers

from antismash import utils

# --------------------------------------------------------------------------
# Constants / cache
# --------------------------------------------------------------------------

PWM_PATH = utils.get_full_path(__file__, os.path.join("data", "Athaliana_motifs.filtered.json"))
_MATRIX_CACHE: Dict[str, List["Matrix"]] = {}   # cache parsed matrices per file path


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------

class Confidence(IntEnum):
    WEAK = auto()
    MEDIUM = auto()
    STRONG = auto()

    def __str__(self) -> str:
        return self.name.lower()


@dataclass
class Matrix:
    name: str
    pwm: List[List[float]]           # 4×N, PFM or log-odds
    max_score: float
    min_score: float
    description: str
    species: str
    link: str
    consensus: str
    _threshold: float = -1.0
    is_log_odds: bool = False        # if True, pwm already is log-odds

    @property
    def score_threshold(self) -> float:
        if self._threshold < 0:
            self._threshold = (self.min_score + self.max_score) / 2
        return self._threshold

    def get_score_confidence(self, score: float) -> Confidence:
        if score <= self.min_score:
            return Confidence.WEAK
        if score >= self.score_threshold:
            return Confidence.STRONG
        return Confidence.MEDIUM

    def to_json(self) -> Dict[str, Any]:
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    @staticmethod
    def from_json(name: str, data: Dict[str, Any]) -> "Matrix":
        return Matrix(
            name=name,
            pwm=data["pwm"],
            max_score=data.get("max_score", 0.0),
            min_score=data.get("min_score", 0.0),
            description=data.get("description", ""),
            species=data.get("species", ""),
            link=data.get("link", ""),
            consensus=data.get("consensus", ""),
            is_log_odds=data.get("is_log_odds", False),
        )


@dataclass
class TFBSHit:
    name: str
    start: int                  # absolute genomic coord on record (0-based)
    species: str
    link: str
    description: str
    consensus: str
    confidence: Confidence
    strand: int                 # +1 / -1
    score: float
    max_score: float

    def to_json(self) -> Dict[str, Any]:
        data = dict(vars(self))
        data["confidence"] = str(data["confidence"])
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "TFBSHit":
        d = dict(data)
        d["confidence"] = Confidence[d["confidence"].upper()]
        return TFBSHit(**d)


class TFBSFinderResults:
    schema_version = 1

    def __init__(self, record_id: str, pvalue: float, start_overlap: int,
                 hits_by_record: Dict[str, List[TFBSHit]]) -> None:
        self.record_id = record_id
        self.pvalue = pvalue
        self.start_overlap = start_overlap
        self.hits_by_record = hits_by_record

    def to_json(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "pvalue": self.pvalue,
            "start_overlap": self.start_overlap,
            "hits_by_record": {
                k: [hit.to_json() for hit in v] for k, v in self.hits_by_record.items()
            }
        }

    def get_hits_for_record(self, record_id: str,
                            confidence: Optional[Confidence] = None,
                            allow_better: bool = False) -> List[TFBSHit]:
        hits = self.hits_by_record.get(record_id, [])
        if confidence is None:
            return hits
        if allow_better:
            return [h for h in hits if h.confidence >= confidence]
        return [h for h in hits if h.confidence == confidence]

    def format_html(self) -> str:
        out = [f"<h3>TFBS Finder Results</h3>"]
        if not self.hits_by_record:
            out.append("<p>No transcription factor binding sites detected.</p>")
            return "\n".join(out)
        out += [
            "<table class='table table-sm'>",
            "<thead><tr><th>Motif</th><th>Start</th><th>Strand</th><th>Score</th>"
            "<th>Confidence</th><th>Species</th></tr></thead>",
            "<tbody>",
        ]
        for _, hits in self.hits_by_record.items():
            for h in hits:
                strand = "+" if h.strand == 1 else "−"
                out.append(
                    f"<tr><td>{h.name}</td><td>{h.start}</td><td>{strand}</td>"
                    f"<td>{h.score:.2f}/{h.max_score:.2f}</td>"
                    f"<td>{str(h.confidence).capitalize()}</td><td>{h.species}</td></tr>"
                )
        out += ["</tbody></table>"]
        return "\n".join(out)

    @staticmethod
    def from_json(previous: Dict[str, Any], record: SeqRecord) -> Optional["TFBSFinderResults"]:
        try:
            if previous.get("schema_version") != TFBSFinderResults.schema_version:
                return None
            if previous.get("record_id") != record.id:
                return None
            pvalue = float(previous["pvalue"])
            start_overlap = int(previous["start_overlap"])
            hits_by_record: Dict[str, List[TFBSHit]] = {}
            for k, hits in previous.get("hits_by_record", {}).items():
                hits_by_record[str(k)] = [TFBSHit.from_json(h) for h in hits]
            return TFBSFinderResults(
                record_id=previous["record_id"],
                pvalue=pvalue,
                start_overlap=start_overlap,
                hits_by_record=hits_by_record,
            )
        except Exception:
            return None


# --------------------------------------------------------------------------
# Helpers: windows, matrices, MOODS
# --------------------------------------------------------------------------

def _cds_tss_and_strand(cds) -> Tuple[Optional[int], Optional[int]]:
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


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping [a,b] inclusive intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [[intervals[0][0], intervals[0][1]]]
    for s, e in intervals[1:]:
        if s > merged[-1][1] + 1:
            merged.append([s, e])
        else:
            if e > merged[-1][1]:
                merged[-1][1] = e
    return [(int(s), int(e)) for s, e in merged]


def _collect_windows_for_cluster(record: SeqRecord,
                                 cluster_feature,
                                 half_window: int) -> Tuple[List[Tuple[int, int]], int]:
    """
    Build CDS promoter windows for CDS overlapping this cluster and clip to cluster span.
    Returns (raw_intervals_inclusive, cds_count_included).
    """
    seqlen = len(record.seq)
    cstart = int(cluster_feature.location.start)
    cend   = int(cluster_feature.location.end) - 1  # inclusive
    cds_count = 0
    raw: List[Tuple[int, int]] = []

    for cds in utils.get_cds_features(record):
        cds_start = int(cds.location.start)
        cds_end   = int(cds.location.end) - 1  # inclusive
        if cds_end < cstart or cds_start > cend:
            continue  # CDS does not overlap cluster

        tss, _ = _cds_tss_and_strand(cds)
        if tss is None:
            continue

        # Build ±window around TSS, then clip to cluster span and contig
        a = max(0, tss - half_window)
        b = min(seqlen - 1, tss + half_window)
        a = max(a, cstart)
        b = min(b, cend)
        if b >= a:
            raw.append((a, b))
            cds_count += 1

    return raw, cds_count


def _safe_bg_from_seq(seq: Seq) -> List[float]:
    s = str(seq).upper()
    nA = s.count("A")
    nC = s.count("C")
    nG = s.count("G")
    nT = s.count("T")
    total = nA + nC + nG + nT
    if total == 0:
        return [0.25, 0.25, 0.25, 0.25]
    eps = 1e-9
    arr = np.array([nA, nC, nG, nT], dtype=float) + eps
    arr /= arr.sum()
    return arr.tolist()  # A,C,G,T


def _matrix_to_log_odds(matrix: Matrix, background: List[float]) -> List[List[float]]:
    """
    Produce a 4×N log-odds matrix for a motif.
    If matrix.is_log_odds == True, assume matrix.pwm already is log-odds.
    Otherwise, treat matrix.pwm as PFM and convert using MOODS with the *given* background.
    """
    pwm = matrix.pwm if not isinstance(matrix.pwm, np.ndarray) else matrix.pwm.tolist()
    if not pwm or len(pwm) != 4 or any(len(r) != len(pwm[0]) for r in pwm):
        raise ValueError(f"{matrix.name}: PWM must be 4×N")

    # Heuristic: negatives → already log-odds
    if matrix.is_log_odds or any(val < 0 for row in pwm for val in row):
        lod = pwm
    else:
        # Convert PFM -> log-odds via MOODS; needs a temp file
        pfm_str = "\n".join(" ".join(f"{v:.6f}" for v in row) for row in pwm) + "\n"
        tmp = None
        try:
            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pfm", delete=False)
            tmp.write(pfm_str)
            tmp.flush()
            tmp.close()
            lod = parsers.pfm_to_log_odds(tmp.name, background, 1e-3)
            lod = lod.tolist() if isinstance(lod, np.ndarray) else lod
        finally:
            if tmp is not None:
                try:
                    os.unlink(tmp.name)
                except Exception:
                    pass

    if not lod or len(lod) != 4 or any(len(r) != len(lod[0]) for r in lod):
        raise ValueError(f"{matrix.name}: log-odds must be 4×N")
    return lod


def _load_matrices_cached(json_file: str) -> List[Matrix]:
    mats = _MATRIX_CACHE.get(json_file)
    if mats is not None:
        return mats
    with open(json_file, encoding="utf-8") as fh:
        data = json.load(fh)
    mats = []
    for name, values in data.items():
        try:
            m = Matrix.from_json(name, values)
            if len(m.pwm) != 4 or any(len(row) != len(m.pwm[0]) for row in m.pwm):
                raise ValueError("PWM must be 4×N")
            mats.append(m)
        except Exception as e:
            logging.error("Skipping motif %r due to parse/shape error: %s", name, e)
    _MATRIX_CACHE[json_file] = mats
    logging.debug("Loaded %d matrices from %s", len(mats), json_file)
    return mats


def _scan_segment_with_pwms(record: SeqRecord,
                            seg_a: int,
                            seg_b: int,
                            matrices: List[Matrix],
                            pvalue: float) -> List[Tuple[int, int, int, float]]:
    """
    Scan record.seq[seg_a:seg_b+1] once across all PWMs and return raw hits:
      List[(matrix_idx, absolute_start, strand(+1/-1), score)]
    Uses per-segment background and per-motif MOODS thresholds from p-value.
    """
    seq = record.seq[seg_a:seg_b+1]
    background = _safe_bg_from_seq(seq)  # A,C,G,T
    hits: List[Tuple[int, int, int, float]] = []
    rc_seq = str(Seq(str(seq)).reverse_complement())
    fwd_seq = str(seq)

    # Debug throttle (optional): limit number of motifs via env
    max_pwms = int(os.environ.get("TFBS_MAX_MOTIFS", "0") or "0")
    mats = matrices[:max_pwms] if max_pwms > 0 else matrices

    for idx, m in enumerate(mats):
        try:
            lod = _matrix_to_log_odds(m, background)
            motif_len = len(lod[0])

            thr = moods_tools.threshold_from_p(lod, background, pvalue)
            thresholds = [thr]

            fwd = scan.scan_dna(fwd_seq, [lod], background, thresholds, 7)[0]
            for mh in fwd:
                hits.append((idx, seg_a + mh.pos, 1, mh.score))

            rev = scan.scan_dna(rc_seq, [lod], background, thresholds, 7)[0]
            for mh in rev:
                abs_pos = seg_a + (len(fwd_seq) - mh.pos - motif_len)
                hits.append((idx, abs_pos, -1, mh.score))

        except Exception as e:
            logging.error("MOODS failed for %s: %s", m.name, e)

    return hits


def _filter_hits_to_objects(matrices: List[Matrix],
                            raw_hits: List[Tuple[int, int, int, float]]) -> List[TFBSHit]:
    out: List[TFBSHit] = []
    for mat_idx, start, strand, score in raw_hits:
        m = matrices[mat_idx]
        conf = m.get_score_confidence(score)
        out.append(
            TFBSHit(
                name=m.name,
                start=int(start),
                species=m.species,
                link=m.link,
                description=m.description,
                consensus=m.consensus,
                confidence=conf,
                strand=int(strand),
                score=float(score),
                max_score=float(m.max_score),
            )
        )
    return out


# --------------------------------------------------------------------------
# Public entry point (per-BGC scanning only)
# --------------------------------------------------------------------------

def run_tfbs_finder(record: SeqRecord,
                    pvalue: float,
                    start_overlap: int,
                    matrix_path: str = PWM_PATH) -> TFBSFinderResults:
    """
    Run TFBS scan **per BGC** on this record:
      - for each cluster on the record, build ±start_overlap windows around CDS TSS,
        clipped to the cluster span;
      - merge windows inside that cluster and scan each merged interval once;
      - aggregate/deduplicate hits across clusters.

    Returns TFBSFinderResults with hits_by_record = { record.id: [TFBSHit, ...] }.
    """
    logging.info("TFBS: %s starting (per-BGC mode)", record.id)

    matrices = _load_matrices_cached(matrix_path)
    if not matrices:
        logging.warning("TFBS: no matrices loaded (%s)", matrix_path)
        return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: []})

    clusters = utils.get_sorted_cluster_features(record)
    if not clusters:
        logging.info("TFBS: %s has no clusters; nothing to scan", record.id)
        return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: []})

    all_raw_hits: List[Tuple[int, int, int, float]] = []
    total_bp = 0
    total_int = 0

    for c in clusters:
        cidx = utils.get_cluster_number(c)
        cstart = int(c.location.start)
        cend   = int(c.location.end) - 1  # inclusive

        raw_windows, cds_count = _collect_windows_for_cluster(record, c, start_overlap)
        if not raw_windows:
            logging.info("TFBS: cluster #%d %d-%d: no CDS windows; skip",
                         cidx, cstart, cend)
            continue

        merged = _merge_intervals(raw_windows)
        bp = sum(b - a + 1 for a, b in merged)
        total_bp += bp
        total_int += len(merged)
        logging.info("TFBS: cluster #%d %d-%d: CDS windows=%d; merged intervals=%d; bp=%d",
                     cidx, cstart, cend, cds_count, len(merged), bp)

        # Optional debug throttle: limit intervals per cluster
        max_intervals = int(os.environ.get("TFBS_MAX_INTERVALS", "0") or "0")
        intervals = merged[:max_intervals] if max_intervals > 0 else merged

        scanned_bp = 0
        for j, (a, b) in enumerate(intervals, 1):
            seg_hits = _scan_segment_with_pwms(record, a, b, matrices, pvalue)
            all_raw_hits.extend(seg_hits)
            scanned_bp += (b - a + 1)
            if j % 20 == 0 or j == len(intervals):
                logging.info("TFBS: cluster #%d scanned %d/%d intervals (%.1f%%), bp=%d, hits=%d",
                             cidx, j, len(intervals), 100.0 * j / len(intervals),
                             scanned_bp, len(all_raw_hits))

    # Deduplicate across clusters (same motif, start, strand, score)
    if all_raw_hits:
        as_set = {}
        for mi, st, sd, sc in all_raw_hits:
            as_set[(mi, st, sd, round(sc, 4))] = (mi, st, sd, sc)
        all_raw_hits = list(as_set.values())

    logging.info("TFBS: scanned total ~%d bp across %d merged intervals × %d motifs; raw hits=%d",
                 total_bp, total_int, len(matrices), len(all_raw_hits))

    logging.warning("🧪 Filtering hits")
    hits = _filter_hits_to_objects(matrices, all_raw_hits)

    # Attach to record.annotations for downstream HTML (best effort)
    try:
        record.annotations.setdefault("tfbs_finder", {})[record.id] = [h.to_json() for h in hits]
    except Exception as e:
        logging.debug("TFBS: could not attach hits to record annotations: %s", e)

    logging.info("TFBS: %s finished with %d hit record(s)", record.id, len(hits))
    return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: hits})