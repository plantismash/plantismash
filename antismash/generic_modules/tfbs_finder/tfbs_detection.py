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
from MOODS import tools, scan

from collections import Counter
from antismash import utils 
import json
import os 
from argparse import Namespace

PWM_PATH = utils.get_full_path(__file__, os.path.join("data", "Athaliana_subset.json"))

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
        return {key: val for key, val in vars(self).items() if not key.startswith("_")}

    @staticmethod
    def from_json(name: str, data: Dict[str, Any]) -> "Matrix":
        return Matrix(name, **data)

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
    def from_json(data: Dict[str, Any], record: SeqRecord) -> Optional["TFBSFinderResults"]:
        if data.get("record_id") != record.id:
            logging.warning("Record ID mismatch for TFBSFinderResults")
            return None
        hits = {k: [TFBSHit.from_json(h) for h in v] for k, v in data["hits_by_record"].items()}
        return TFBSFinderResults(data["record_id"], data["pvalue"], data["start_overlap"], hits)

### Utility functions ###

def load_matrices(json_file: str) -> List[Matrix]:
    with open(json_file) as handle:
        data = json.load(handle)
    return [Matrix.from_json(name, values) for name, values in data.items()]

def get_bg_distribution(seq: Seq) -> List[float]:
    counts = Counter(seq.upper())
    total = sum(counts[base] for base in "ACGT")
    if total == 0:
        return [0.25, 0.25, 0.25, 0.25]
    return [counts[base] / total for base in "ACGT"]

def run_moods(seq: Seq, background: List[float], matrices: List[Matrix], pvalue: float, offset: int) -> List[Tuple[int, int, int, float]]:
    seq = seq.upper()
    mood_pwms = [tools.parse_matrix(m.pwm) for m in matrices]
    tools.set_background(background)
    thresholds = [tools.threshold_from_p(m, background, pvalue) for m in mood_pwms]
    results = scan.scan_dna(seq, mood_pwms, thresholds)
    hits = []
    for i, matrix_hits in enumerate(results):
        for hit_start, score in matrix_hits:
            hits.append((i, offset + hit_start, 1, score))
    rc_seq = str(Seq(seq).reverse_complement())
    rc_results = scan.scan_dna(rc_seq, mood_pwms, thresholds)
    for i, matrix_hits in enumerate(rc_results):
        for hit_start, score in matrix_hits:
            hits.append((i, offset + len(seq) - hit_start - len(matrices[i].pwm), -1, score))
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

def run_tfbs_finder(record: SeqRecord, pvalue: float, start_overlap: int,
                    matrix_path: str = PWM_PATH) -> TFBSFinderResults:
    matrices = load_matrices(matrix_path)
    seq = record.seq
    hits_by_record = {}
    background = get_bg_distribution(seq)
    moods_results = run_moods(seq, background, matrices, pvalue, 0)
    hits_by_record[record.id] = filter_hits(matrices, len(seq), moods_results)
    return TFBSFinderResults(record.id, pvalue, start_overlap, hits_by_record)

### Module entry point for plantiSMASH ###

def run_analyses(seq_records: List[SeqRecord], options) -> None:
    for record in seq_records:
        if not record.seq or len(record.seq) == 0:
            continue
        results = run_tfbs_finder(record, options.tfbs_pvalue, options.tfbs_range)
        if record.id not in options.extrarecord:
            options.extrarecord[record.id] = Namespace()
        if not hasattr(options.extrarecord[record.id], "extradata"):
            options.extrarecord[record.id].extradata = {}
        options.extrarecord[record.id].extradata["TFBSFinderResults"] = results
