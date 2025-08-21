# License: GNU AGPL v3 or later

"""
Genome-wide module to detect transcription factor binding sites (TFBSs)
using position weight matrices (PWMs).
"""

import logging
from dataclasses import dataclass
from enum import IntEnum, auto
from typing import Any, Dict, List, Optional, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation

import tempfile 
from MOODS import tools as moods_tools
from MOODS import scan, parsers 
import numpy as np 

from collections import Counter
from antismash import utils 
import json
import os 

PWM_PATH = utils.get_full_path(__file__, os.path.join("data", "Athaliana_motifs.filtered.json"))

class Confidence(IntEnum):
    WEAK = auto()
    MEDIUM = auto()
    STRONG = auto()

    def __str__(self) -> str:
        return self.name.lower()

@dataclass
class Matrix:
    name: str
    pwm: List[List[float]]
    max_score: float
    min_score: float
    description: str
    species: str
    link: str
    consensus: str
    _threshold: float = -1.
    is_log_odds: bool = False

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
    start: int
    species: str
    link: str
    description: str
    consensus: str
    confidence: Confidence
    strand: int
    score: float
    max_score: float

    def to_json(self) -> Dict[str, Any]:
        data = dict(vars(self))
        data["confidence"] = str(data["confidence"])
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "TFBSHit":
        data = dict(data)
        data["confidence"] = Confidence[data["confidence"].upper()]
        return TFBSHit(**data)

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
        """Simple HTML output for plantiSMASH final page"""
        output = [f"<h3>TFBS Finder Results</h3>"]
        if not self.hits_by_record:
            output.append("<p>No transcription factor binding sites detected.</p>")
            return "\n".join(output)

        output.append("<table class='table table-sm'>")
        output.append("<thead><tr><th>Motif</th><th>Start</th><th>Strand</th><th>Score</th><th>Confidence</th><th>Species</th></tr></thead>")
        output.append("<tbody>")

        for record_id, hits in self.hits_by_record.items():
            for hit in hits:
                strand = "+" if hit.strand == 1 else "−"
                output.append(
                    f"<tr><td>{hit.name}</td><td>{hit.start}</td><td>{strand}</td><td>{hit.score:.1f}/{hit.max_score:.1f}</td>"
                    f"<td>{str(hit.confidence).capitalize()}</td><td>{hit.species}</td></tr>"
                )

        output.append("</tbody></table>")
        return "\n".join(output)

    @staticmethod
    def from_json(previous: Dict[str, Any], record: SeqRecord) -> Optional["TFBSFinderResults"]:
        """Rebuild results from JSON; return None if schema/options/record mismatch."""
        try:
            if previous.get("schema_version") != TFBSFinderResults.schema_version:
                return None
            if previous.get("record_id") != record.id:
                return None

            pvalue = float(previous["pvalue"])
            start_overlap = int(previous["start_overlap"])

            hits_by_record: Dict[str, List[TFBSHit]] = {}
            for key, hits in previous.get("hits_by_record", {}).items():
                hits_by_record[str(key)] = [TFBSHit.from_json(h) for h in hits]

            return TFBSFinderResults(
                record_id=previous["record_id"],
                pvalue=pvalue,
                start_overlap=start_overlap,
                hits_by_record=hits_by_record,
            )
        except Exception:
            return None

### Utility functions ###

def load_matrices(json_file: str) -> List[Matrix]:
    with open(json_file, encoding="utf-8") as handle:
        data = json.load(handle)

    mats: List[Matrix] = []
    for name, values in data.items():
        try:
            m = Matrix.from_json(name, values)
            # quick shape sanity check: 4 x N
            if len(m.pwm) != 4 or any(len(row) != len(m.pwm[0]) for row in m.pwm):
                raise ValueError("PWM must be 4×N")
            mats.append(m)
        except Exception as e:
            logging.error("Skipping motif %r due to parse/shape error: %s", name, e)
            continue

    logging.debug("Loaded %d matrices from %s", len(mats), json_file)
    return mats


def get_bg_distribution(seq: Seq) -> Tuple[float, float, float, float]:
    s = str(seq).upper()
    gc = (s.count('G') + s.count('C')) / (2.0 * max(1, len(s)))
    at = 0.5 - gc
    return (at, gc, gc, at)  # A, C, G, T

def matrix_to_pfm_string(pwm_4xN: List[List[float]]) -> str:
    """Return a JASPAR-style PFM string (4 lines: A,C,G,T; N space-separated values)."""
    # use spaces (not tabs) and ensure trailing newline; MOODS is picky on parsing
    lines = [" ".join(f"{v:.6f}" for v in row) for row in pwm_4xN]
    return "\n".join(lines) + "\n"


def pwm_to_string(pwm: List[List[float]]) -> str:
    """Convert PWM (list of lists) to string format for MOODS."""
    return "\n".join("\t".join(f"{val:.6f}" for val in row) for row in pwm)

def _prepare_log_odds(pwm_like: List[List[float]], background: List[float], name: str) -> List[List[float]]:
    """Return a 4×N log-odds matrix for MOODS.
       - Accepts list-of-lists or numpy arrays.
       - If values look like PFMs (no negatives), converts via MOODS parsers.
       - Otherwise assumes they are already log-odds.
    """
    # Normalize type
    pwm = pwm_like
    if isinstance(pwm, np.ndarray):
        pwm = pwm.tolist()

    # Basic 4×N shape check on input
    if not pwm or len(pwm) != 4 or len(pwm[0]) == 0 or any(len(r) != len(pwm[0]) for r in pwm):
        raise ValueError(f"{name}: PWM must be 4×N")

    # Heuristic: log-odds usually contain negatives; PFMs typically do not
    looks_like_log_odds = any(val < 0 for row in pwm for val in row)

    if looks_like_log_odds:
        lod = pwm
    else:
        # Convert PFM → log-odds using MOODS; needs a temp file
        pfm_str = matrix_to_pfm_string(pwm)
        tmp = None
        try:
            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pfm", delete=False)
            tmp.write(pfm_str)
            tmp.flush()
            tmp.close()
            lod = parsers.pfm_to_log_odds(tmp.name, background, 1e-3)
            if isinstance(lod, np.ndarray):
                lod = lod.tolist()
        finally:
            if tmp is not None:
                try:
                    os.unlink(tmp.name)
                except Exception:
                    pass

    # Validate 4×N on output
    if not lod or len(lod) != 4 or len(lod[0]) == 0 or any(len(r) != len(lod[0]) for r in lod):
        raise ValueError(f"{name}: log-odds must be 4×N")

    return lod


def run_moods(seq: Seq,
              background: List[float],
              matrices: List[Matrix],
              pvalue: float,
              offset: int) -> List[Tuple[int, int, int, float]]:

    hits: List[Tuple[int, int, int, float]] = []
    seq = str(seq).upper()
    rc_seq = str(Seq(seq).reverse_complement())

    for idx, matrix in enumerate(matrices):
        logging.debug("🔍 Scanning for matrix: %s", matrix.name)

        # Prepare a 4×N log-odds matrix (convert if needed)
        try:
            lod = _prepare_log_odds(matrix.pwm, background, matrix.name)
        except Exception as e:
            logging.error("❌ preparing log-odds for %s failed: %s", matrix.name, e)
            continue

        motif_len = len(lod[0])
        logging.debug("✅ Log-odds for %s: 4×%d", matrix.name, motif_len)

        # Compute per-motif threshold from p-value (on the *same* log-odds)
        try:
            thr = moods_tools.threshold_from_p(lod, background, pvalue)
        except Exception as e:
            logging.error("❌ threshold_from_p failed for %s: %s", matrix.name, e)
            continue

        thresholds = [thr]

        # Forward strand
        try:
            fwd = scan.scan_dna(seq, [lod], background, thresholds, 7)[0]
        except Exception as e:
            logging.error("❌ scan_dna (forward) failed for %s: %s", matrix.name, e)
            fwd = []
        for mh in fwd:
            hits.append((idx, offset + mh.pos, 1, mh.score))

        # Reverse strand (report forward coordinates for starts)
        try:
            rev = scan.scan_dna(rc_seq, [lod], background, thresholds, 7)[0]
        except Exception as e:
            logging.error("❌ scan_dna (reverse) failed for %s: %s", matrix.name, e)
            rev = []
        for mh in rev:
            real_pos = offset + (len(seq) - mh.pos - motif_len)
            hits.append((idx, real_pos, -1, mh.score))

    return hits



def filter_hits(matrices: List[Matrix], seqlen: int, hits: List[Tuple[int, int, int, float]]) -> List[TFBSHit]:
    tfbs_hits = []
    for matrix_idx, start, strand, score in hits:
        matrix = matrices[matrix_idx]
        conf = matrix.get_score_confidence(score)
        tfbs_hits.append(TFBSHit(
            name=matrix.name,
            start=start,
            species=matrix.species,
            link=matrix.link,
            description=matrix.description,
            consensus=matrix.consensus,
            confidence=conf,
            strand=strand,
            score=score,
            max_score=matrix.max_score
        ))
    return tfbs_hits


### Core analysis ###

def _promoter_windows_around_gene_starts(record: SeqRecord, half_window: int) -> List[Tuple[int, int]]:
    """Return merged [start,end) windows of size ±half_window around CDS start sites.
       For + strand, start site is CDS.start; for – strand, it's CDS.end-1.
    """
    seqlen = len(record.seq)
    starts: List[int] = []

    for feat in getattr(record, "features", []):
        if getattr(feat, "type", None) != "CDS":
            continue
        loc = getattr(feat, "location", None)
        if loc is None:
            continue
        strand = int(getattr(loc, "strand", 1) or 1)

        # Handle split genes too
        if isinstance(loc, CompoundLocation) and loc.parts:
            first, last = loc.parts[0], loc.parts[-1]
            s = int(first.start) if strand == 1 else int(last.end) - 1
        else:
            s = int(loc.start) if strand == 1 else int(loc.end) - 1

        starts.append(s)

    # Build and merge intervals
    intervals = []
    for s in starts:
        a = max(0, s - half_window)
        b = min(seqlen, s + half_window + 1)  # half-open
        intervals.append((a, b))
    intervals.sort()

    merged: List[Tuple[int, int]] = []
    for a, b in intervals:
        if not merged or a > merged[-1][1]:
            merged.append([a, b])
        else:
            merged[-1][1] = max(merged[-1][1], b)

    return [(int(a), int(b)) for a, b in merged]

def run_tfbs_finder(record: SeqRecord, pvalue: float, start_overlap: int,
                    matrix_path: str = PWM_PATH,
                    region: Optional[Tuple[int, int]] = None) -> TFBSFinderResults:
    logging.warning("📥 Starting run_tfbs_finder for record: %s", record.id)

    logging.debug("🔢 Loading matrices from %s", matrix_path)
    matrices = load_matrices(matrix_path)
    logging.debug("✅ Loaded %d matrices", len(matrices))

    # Decide scanning areas
    if region:
        areas = [region]
        logging.debug("📏 Using specified region: %d-%d", region[0], region[1])
    else:
        if start_overlap and start_overlap > 0:
            areas = _promoter_windows_around_gene_starts(record, start_overlap)
            if not areas:
                areas = [(0, len(record.seq))]
                logging.debug("🔎 No CDS features found; scanning full sequence.")
        else:
            areas = [(0, len(record.seq))]

    total_bp = sum(b - a for a, b in areas)
    logging.debug("🗺️ Scanning %d window(s) totalling %d bp (of %d bp)",
                  len(areas), total_bp, len(record.seq))
    logging.debug("📊 p-value for thresholds: %.3g", pvalue)

    # Run MOODS per window (background per window)
    all_hits: List[Tuple[int, int, int, float]] = []
    for a, b in areas:
        sub_seq = record.seq[a:b]
        background = get_bg_distribution(sub_seq)  # (A,C,G,T) for this window
        logging.debug("  • Window %d-%d: bg=%s", a, b, background)
        window_hits = run_moods(sub_seq, background, matrices, pvalue, offset=a)
        all_hits.extend(window_hits)

    logging.warning("🧪 Filtering hits")
    hits_by_record = {
        record.id: filter_hits(matrices, len(record.seq), all_hits)
    }

    logging.warning("🏁 run_tfbs_finder finished for record: %s", record.id)

    # Attach JSON-serialisable TFBS hits to the record annotations for HTML panel
    try:
        record.annotations.setdefault("tfbs_finder", {})[record.id] = [
            hit.to_json() for hit in hits_by_record[record.id]
        ]
        logging.debug("🔗 Attached %d TFBS hits to record.annotations['tfbs_finder']",
                      len(hits_by_record[record.id]))
    except Exception as e:
        logging.warning("Could not attach TFBS hits to record: %s", e)


    return TFBSFinderResults(record.id, pvalue, start_overlap, hits_by_record)




