from antismash import utils
from . import repeatfinder

def specific_analysis(seq_record, options):
    clusters = utils.get_cluster_features(seq_record)

    # Step 1: run repeatfinder on full record BEFORE filtering
    find_repeats(seq_record)

    valid_clusters = []

    # Step 2: after repeats are annotated, do filtering
    for cluster in clusters:
        if 'product' not in cluster.qualifiers or 'cyclopeptide' not in cluster.qualifiers['product'][0]:
            valid_clusters.append(cluster)
            continue

        if cluster_has_burp_and_known_repeat(seq_record, cluster, options.require_internal_cyclopeptide_repeats):
            print(f"Keeping cyclopeptide cluster: {cluster.qualifiers['product'][0]}")
            valid_clusters.append(cluster)
        else:
            print(f"Skipping cyclopeptide cluster: {cluster.qualifiers['product'][0]}")

    # Step 3: remove only the *clusters* not matching
    seq_record.features = [f for f in seq_record.features if f not in clusters or f in valid_clusters]


def find_repeats(seq_record):
    repeatfinder.run_fbk(seq_record)
    fbk_output_to_result(seq_record)


def fbk_output_to_result(seq_record):
    for feat in seq_record.features:
        result = Result()
        has_ripp = False

        if "seq_met" in feat.qualifiers:
            result.feature_type = feat.qualifiers["seq_met"][0]
        else:
            result.feature_type = feat.type

        if "pattern" in feat.qualifiers and "table" in feat.qualifiers:
            if len(feat.qualifiers["table"][0]) > 5 and len(feat.qualifiers.get("top_kmer_hits", [])) >= 3:
                has_ripp = True
                result.pattern = feat.qualifiers["pattern"]
                result.instances = feat.qualifiers["table"]

        if "ripp_evidence" in feat.qualifiers and len(feat.qualifiers["ripp_evidence"]) > 0:
            has_ripp = True
            result.evidence = feat.qualifiers["ripp_evidence"]

        if has_ripp:
            result.sequence = feat.qualifiers.get("translation", [""])[0]
            result.position = (feat.location.start, feat.location.end)

            if 'gene' in feat.qualifiers:
                result.cds_id = feat.qualifiers['gene'][0]
            elif 'locus_tag' in feat.qualifiers:
                result.cds_id = feat.qualifiers['locus_tag'][0]
            elif 'db_xref' in feat.qualifiers:
                result.cds_id = feat.qualifiers['db_xref'][0]
            else:
                result.cds_id = f"{result.position[0]}-{result.position[1]}"

            # attach
            feat.qualifiers['cyclopeptide_analysis'] = [result.encode()]


def cluster_has_burp_and_known_repeat(seq_record, cluster, require_internal):
    has_burp = False
    has_internal_repeat = False
    has_external_repeat = False

    cluster_start = cluster.location.start
    cluster_end = cluster.location.end

    for feature in seq_record.features:
        if not (cluster_start <= feature.location.start <= cluster_end):
            continue

        domains = feature.qualifiers.get("domain", [])
        domain_record = feature.qualifiers.get("domain_record", [])
        if isinstance(domain_record, str):
            domain_record = [domain_record]

        is_burp = any("BURP" in d for d in domains) or any("BURP" in d for d in domain_record)

        has_repeat = (
            feature.qualifiers.get('has_repeat') == True or
            "cyclopeptide_analysis" in feature.qualifiers
        )

        if is_burp:
            has_burp = True
            if has_repeat:
                has_internal_repeat = True
        else:
            if has_repeat:
                has_external_repeat = True

    if not has_burp:
        return False

    if require_internal:
        return has_internal_repeat
    else:
        return has_internal_repeat or has_external_repeat


class Result:
    pattern = ""
    instances = []
    sequence = ""
    feature_type = ""
    position = None
    cds_id = None
    evidence = {}

    def __init__(self, qualifier=None):
        if qualifier is not None:
            self.decode(qualifier)

    def encode(self):
        qualifier = ""
        qualifier += self.pattern
        qualifier += "//" + "|".join(str(s) for s in self.instances)
        qualifier += "//" + self.sequence
        qualifier += "//" + self.feature_type
        qualifier += "//" + str(self.position[0]) + "|" + str(self.position[1])
        qualifier += "//" + self.cds_id
        qualifier += "//" + "|".join(self.evidence.keys())
        qualifier += "//" + "|".join([str(x) for x in self.evidence.values()])
        return qualifier

    def decode(self, qualifier):
        self.__init__()
        q = qualifier.split("//")
        self.pattern = q[0]
        self.instances = q[1].split("|")
        self.sequence = q[2]
        self.feature_type = q[3]
        self.position = tuple(q[4].split("|"))
        self.cds_id = q[5]
        keys = q[6].split("|")
        values = [list(v[2:-2].split(",")) for v in q[7].split("|")]
        self.evidence = dict(zip(keys, values))
