# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Module to detect transcription factor binding sites (TFBSs)
using position weight matrices (PWMs) in BGCs.
"""

from dataclasses import dataclass
from enum import IntEnum, auto
import itertools
import logging

from typing import Any, Dict, Iterator, List, Optional, Tuple

from collections import Counter
from Bio.Seq import Seq
from MOODS import tools, scan

from antismash.config import get_config
from antismash.common import json as jsonlib, path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.features import Feature, Region
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
)

PWM_PATH = path.get_full_path(__file__, 'data', 'PWMs.json')

class Confidence(IntEnum):
    """ Defined and sortable values for confidence terms """
    WEAK = auto()
    MEDIUM = auto()
    STRONG = auto()

    def __str__(self) -> str:
        return self.name.lower()

@dataclass
class Matrix:
    """ Class to store position weight matrices (Matrix) """
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
        """ The minimum threshold required for a hit's score to be considered good """
        if self._threshold < 0:
            self._threshold = (self.min_score + self.max_score) / 2
        return self._threshold

    def get_score_confidence(self, score: float) -> Confidence:
        """ Returns the confidence level for the given score """
        if score <= self.min_score:
            return Confidence.WEAK
        if score >= self.score_threshold:
            return Confidence.STRONG
        return Confidence.MEDIUM

    def to_json(self) -> Dict[str, Any]:
        """ Returns a JSON-ready representation of the instance """
        return {key: val for key, val in vars(self).items() if not key.startswith("_")}

    @staticmethod
    def from_json(name: str, data: Dict[str, Any]) -> "Matrix":
        """ Reconstructs an instance from a JSON representation """
        return Matrix(name, **data)

@dataclass
class TFBSHit:
    """ Class to store the transcription factor binding site hits """
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
        """ Returns a JSON-ready representation of the instance """
        data = dict(vars(self))
        data["confidence"] = str(data["confidence"])
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "TFBSHit":
        """ Reconstructs an instance from a JSON representation """
        data = dict(data)
        data["confidence"] = Confidence[data["confidence"].upper()]
        return TFBSHit(**data)

class TFBSFinderResults(ModuleResults):
    """ Results class for the TFBS finder analysis """
    schema_version = 1

    def __init__(self, record_id: str, pvalue: float, start_overlap: int,
                 hits_by_region: Dict[int, List[TFBSHit]]) -> None:
        super().__init__(record_id)
        self.pvalue = pvalue
        self.start_overlap = start_overlap
        self.hits_by_region = hits_by_region

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """
        hits_by_region = {region: [hit.to_json() for hit in hits] for region, hits in self.hits_by_region.items()}
        return {
            "schema_version": self.schema_version,
            "pvalue": self.pvalue,
            "start_overlap": self.start_overlap,
            "record_id": self.record_id,
            "hits_by_region": hits_by_region,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["TFBSFinderResults"]:
        """ Constructs a new results instance from a JSON format and the original record analysed. """
        if json["schema_version"] != TFBSFinderResults.schema_version:
            return None
        if record.id != json.get("record_id"):
            logging.warning("TFBS results are for different record, discarding previous results")
            return None
        hits_by_region = {int(name): [TFBSHit.from_json(hit) for hit in hits] for name, hits in json["hits_by_region"].items()}
        return TFBSFinderResults(json["record_id"], json["pvalue"], json["start_overlap"], hits_by_region)

def run_tfbs_finder(record: Record, pvalue: float, start_overlap: int, matrix_path: str = PWM_PATH) -> TFBSFinderResults:
    """Run TFBS finder on a given record"""
    hits_by_region = {}
    matrices = load_matrices(matrix_path)
    for region in record.get_regions():
        sequence = region.extract(record.seq)
        background = get_bg_distribution(sequence)
        moods_results = run_moods(sequence, background, matrices, pvalue, region.start)
        areas = get_valid_areas(region.start, region.end, iter(region.cds_children), start_overlap)
        hits_by_region[region.get_region_number()] = filter_hits(matrices, areas, moods_results)
        
    #returns the visualisation 
    return TFBSFinderResults(record.id, pvalue, start_overlap, hits_by_region)
