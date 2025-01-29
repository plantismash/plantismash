# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""General HMM detection module

"""
import logging
from os import path
from os import listdir
from antismash.generic_modules import Signature
from antismash import config
from antismash import utils
from Bio.SeqFeature import SeqFeature, FeatureLocation

name = "hmmdetection"
short_description = name.capitalize()
priority = 5000



class HmmSignature(Signature):
    """HMM signature"""
    def __init__(self, name, description, cutoff, hmm_file):
        self.hmm_file = utils.get_full_path(__file__, hmm_file)
        self.name = name
        super(HmmSignature, self).__init__(name, 'model',
              description, cutoff, utils.get_full_path(__file__, hmm_file))

short_description = name.capitalize()

# The tuple is the name of the binary and whether it is an optional requirement
_required_binaries = [
        ('hmmsearch', False),
    ]

#Define all profiles using details from hmmdetails.txt
#hmmdetails = [line.split("\t") for line in open(utils.get_full_path(__file__, "hmmdetails.txt"),"r").read().split("\n") if line.count("\t") == 3]
#_signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]

def check_prereqs():
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    cfg = config.get_config()
    cfg.cdh_use_binary = True

    #Check if cd-hit binary exist
    if cfg.enable_cdhit:
        if utils.locate_executable("cd-hit") is None:
            cfg.cdh_use_binary = False
            logging.info("Cannot locate cd-hit binary, using ps-cd-hit algorithm instead")

    for hmm_model in cfg.enabled_detection_models:
        #Check if hmmdetails.txt is readable and well-formatted
        dir_path = path.dirname(path.abspath(__file__))
        model_name = ""
        if hmm_model != "default":
            dir_path = path.join(dir_path, hmm_model)
            model_name = "(" + hmm_model + ")"

        #Check if hmmdetails.txt is readable and well-formatted
        lineno = 1
        for line in open(path.join(dir_path, "hmmdetails.txt"),"r"):
            if line.count("\t") != 3:
                failure_messages.append("Failed to use HMM profile %s from line %s due to misformatting:\n %r" %
                                        (model_name, lineno, line))
            lineno += 1

        #Check if cluster_rules.txt is readable and well-formatted
        lineno = 1
        for line in open(path.join(dir_path, "cluster_rules.txt"),"r"):
            if line.count("\t") != 3:
                failure_messages.append("Failed to use cluster rules %s from the line %s due to misformatting:\n %r" %
                                        (model_name, lineno, line))
            lineno += 1

    for sig in get_sig_profiles():
        if not path.exists(sig.path):
            failure_messages.append("Failed to find HMM profile %r" %
                                    sig.path)
    return failure_messages


def get_supported_detection_models():
    # TODO: This function might be removed, as plants is the only supported model
    "Get a list of all supported detection types"
    detection_types = ["default"]
    for fname in listdir(path.dirname(path.abspath(__file__))):
        if path.isdir(path.join(path.dirname(path.abspath(__file__)), fname)):
            detection_types.append(fname)
    return detection_types


def get_supported_cluster_types():
    "Get a list of all supported cluster types"
    clustertypes = [line.split("\t")[0] for line in open(utils.get_full_path(__file__, 'cluster_rules.txt'), "r")][1:]
    # TODO: Iterating across all directories is not needed, as plants is the only supported type 
    for fname in listdir(path.dirname(path.abspath(__file__))):
        # a way to avoid the pycache directory to be checked 
        if fname == "__pycache__":
            continue
        dir_path = path.join(path.dirname(path.abspath(__file__)), fname)
        if path.isdir(dir_path):
            clustertypes.extend([(fname + "/" + line.split("\t")[0]) for line in open(path.join(dir_path, "cluster_rules.txt"), "r")][1:])
    return clustertypes


def get_sig_profiles():
    cfg = config.get_config()
    _signature_profiles = []
    for hmm_model in cfg.enabled_detection_models:
        dir_path = path.dirname(path.abspath(__file__))
        prefix = ""
        if hmm_model != "default":
            dir_path = path.join(dir_path, hmm_model)
            prefix = hmm_model + "/"
        hmmdetails = [line.split("\t") for line in open(path.join(dir_path, "hmmdetails.txt"),"r").read().split("\n") if line.count("\t") == 3]
        _signature_profiles.extend([HmmSignature(prefix + details[0], details[1], int(details[2]), path.join(dir_path, details[3])) for details in hmmdetails])
    return _signature_profiles


def is_sec_met_feature_of_type(type_, feature):
    "Check if feature is related to a secondary metabolite gene cluster"
    return 'sec_met' in feature.qualifiers and feature.qualifiers['sec_met'][0].startswith(type_)

def find_clusters(seq_record, rulesdict, overlaps):
    #Functions that detects the gene clusters based on the identified core genes
    features = utils.get_cds_features(seq_record)
    clustertype = ""
    clusters = []
    cfg = config.get_config()
    clusternr = cfg.next_clusternr
    last_cutoff = 0
    cluster_cds = []

    for feature in features:
        within_cutoff = False
        if ('sec_met' in feature.qualifiers) and (len([feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat]) > 0):
            feature_start = min(feature.location.start, feature.location.end)
            feature_end = max(feature.location.start, feature.location.end)
            feature_type = [feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat][0].partition("Type: ")[2]
            feature_cutoff = max([rulesdict[value][1] for value in feature_type.split("-")])
            feature_extension = max([rulesdict[value][2] for value in feature_type.split("-")])
            if (cfg.enable_dynamic_cutoff):
                multiply_cutoff = get_dynamic_cutoff_multiplier(utils.get_gene_id(feature), overlaps)
                feature_cutoff = int(feature_cutoff * multiply_cutoff)
                feature_extension = int(feature_extension * multiply_cutoff)
            cluster = None

            if len(clusters) > 0:
                cluster = clusters[-1]
                cluster_start = cluster.location.start
                cluster_end = cluster.location.end
                # Check cutoff
                cutoff = max(last_cutoff, feature_cutoff)
                within_cutoff = feature_start <= cluster_end + cutoff
                within_gene_num_cutoff = (min([abs(overlaps[1][utils.get_gene_id(feature)] - overlaps[1][ncds]) for ncds in cluster_cds]) - 1 <= cfg.gene_num_cutoff)
                if (cfg.gene_num_cutoff_only):
                    within_cutoff = within_gene_num_cutoff
                else:
                    within_cutoff = within_cutoff or within_gene_num_cutoff

            if not within_cutoff:
                if len(clusters) > 0:
                    # Finalize the last extended cluster
                    cluster = clusters[-1]
                    cluster.location = FeatureLocation(max(0, cluster.location.start - cluster.qualifiers['extension'][0]), min(len(seq_record), cluster.location.end + cluster.qualifiers['extension'][0]))
                # Create new cluster
                new_cluster = SeqFeature(FeatureLocation(feature_start, feature_end), type="cluster")
                new_cluster.qualifiers['note'] = ["Cluster number: " + str(clusternr)]
                new_cluster.qualifiers['cutoff'] = [feature_cutoff]
                new_cluster.qualifiers['extension'] = [feature_extension]
                new_cluster.qualifiers['product'] = [feature_type]
                clusters.append(new_cluster)
                cluster = clusters[-1]
                cluster_cds = [utils.get_gene_id(feature)]
                clusternr += 1

            # Update cluster
            last_cutoff = feature_cutoff
            cluster.location = FeatureLocation(min(cluster.location.start, feature_start), max(cluster.location.end, feature_end))
            cluster.qualifiers['cutoff'] = [max(cluster.qualifiers['cutoff'][0], feature_cutoff)]
            cluster.qualifiers['extension']  = [max(cluster.qualifiers['extension'][0], feature_extension)]
            cluster.qualifiers['product'] =  ["-".join(list(set(cluster.qualifiers['product'][0].split('-')) | set(feature_type.split('-'))))]
            if "-" in cluster.qualifiers['product'][0]:
                cluster.qualifiers['product'] = ["-".join([ct for ct in cluster.qualifiers['product'][0].split('-') if ct != "other"])]
            if (utils.get_gene_id(feature) not in cluster_cds):
                cluster_cds.append(utils.get_gene_id(feature))

    if len(clusters) > 0:
        # Finalize the last extended cluster
        cluster = clusters[-1]
        cluster.location = FeatureLocation(max(0, cluster.location.start - cluster.qualifiers['extension'][0]), min(len(seq_record), cluster.location.end + cluster.qualifiers['extension'][0]))

    seq_record.features.extend(clusters)
    cfg.next_clusternr = clusternr


def filter_results(results, results_by_id, overlaps, feature_by_id):
    #Filter results by comparing scores of different models according to filterhmmdetails.txt
    cfg = config.get_config()
    filterhmm_list = []
    for hmm_model in cfg.enabled_detection_models:
        dir_path = path.dirname(path.abspath(__file__))
        prefix = ""
        if hmm_model != "default":
            dir_path = path.join(dir_path, hmm_model)
            prefix = hmm_model + "/"
        for line in open(path.join(dir_path, "filterhmmdetails.txt"),"r").read().split("\n"):
            filterhmms = [(prefix + hmm) for hmm in line.split(",")]
            if filterhmms not in filterhmm_list:
                filterhmm_list.append(filterhmms)
            for cds in list(results_by_id.keys()):
                cdsresults = results_by_id[cds]
                hmmhits = [hit.query_id for hit in cdsresults]
                #Check if multiple competing HMM hits are present
                competing_hits = set(hmmhits) & set(filterhmms)
                if len(competing_hits) > 1:
                    #Identify overlapping hits
                    overlapping_groups = []
                    for hit in cdsresults:
                        for otherhit in [cdsresult for cdsresult in cdsresults if hit != cdsresult]:
                            overlap = len(set(range(hit.hit_start, hit.hit_end)) & set(range(otherhit.hit_start, otherhit.hit_end)))
                            if overlap > 20:
                                added = "n"
                                for group in overlapping_groups:
                                    if hit in group and otherhit in group:
                                        added = "y"
                                        break
                                    elif hit in group and otherhit not in group:
                                        group.append(otherhit)
                                        added = "y"
                                        break
                                    elif hit not in group and otherhit in group:
                                        group.append(hit)
                                        added = "y"
                                        break
                                if added == "n":
                                    overlapping_groups.append([hit, otherhit])
                    #Remove worst-scoring of overlapping hits
                    for group in overlapping_groups:
                        highestscore = max([hit.bitscore for hit in group])
                        hit_with_highestscore = group[[hit.bitscore for hit in group].index(highestscore)]
                        to_delete = [hit for hit in group if hit != hit_with_highestscore]
                        for res in [res for res in results]:
                            if res in to_delete:
                                del results[results.index(res)]
                                del results_by_id[cds][results_by_id[cds].index(res)]
                                if len(results_by_id[cds]) < 1:
                                    del results_by_id[cds]

    #Filter multiple results of the same model within a gene
    for cds in list(results_by_id.keys()):
        best_hit_scores = {}
        to_delete = []
        for hit in results_by_id[cds]:
            if (hit.query_id not in best_hit_scores) or (best_hit_scores[hit.query_id] < hit.bitscore):
                best_hit_scores[hit.query_id] = hit.bitscore
        for hit in results_by_id[cds]:
            if (hit.bitscore < best_hit_scores[hit.query_id]):
                to_delete.append(hit)
        for hit in to_delete:
            del results[results.index(hit)]
            del results_by_id[cds][results_by_id[cds].index(hit)]
            if len(results_by_id[cds]) < 1:
                del results_by_id[cds]

    #Filter results of overlapping genes
    overlap_id_with_result = {}
    for cds in list(results_by_id.keys()):
        if overlaps[1][cds] not in list(overlap_id_with_result.keys()):
            overlap_id_with_result[overlaps[1][cds]] = [cds]
        elif cds not in overlap_id_with_result[overlaps[1][cds]]:
            overlap_id_with_result[overlaps[1][cds]].append(cds)
    for overlap_id in list(overlap_id_with_result.keys()):
        best_hit_scores = {}
        for cds in overlap_id_with_result[overlap_id]:
            for hit in results_by_id[cds]:
                feature = feature_by_id[hit.hit_id]
                if (hit.query_id not in best_hit_scores) or (best_hit_scores[hit.query_id] < abs(feature.location.end - feature.location.start)):
                    best_hit_scores[hit.query_id] = abs(feature.location.end - feature.location.start)
        for cds in overlap_id_with_result[overlap_id]:
            to_delete = []
            for hit in results_by_id[cds]:
                feature = feature_by_id[hit.hit_id]
                if (abs(feature.location.end - feature.location.start) < best_hit_scores[hit.query_id]):
                    to_delete.append(hit)
                else: # filter for filterhmmdetails.txt
                    for filterhmms in filterhmm_list:
                        if hit.query_id not in filterhmms:
                            continue
                        for similar_hit in filterhmms:
                            if similar_hit not in list(best_hit_scores.keys()):
                                continue
                            if (abs(feature.location.end - feature.location.start) < best_hit_scores[similar_hit]):
                                to_delete.append(hit)
                                break
            for hit in to_delete:
                del results[results.index(hit)]
                del results_by_id[cds][results_by_id[cds].index(hit)]
                if len(results_by_id[cds]) < 1:
                    del results_by_id[cds]

    return results, results_by_id

def create_rules_dict(enabled_clustertypes):
    "Create a cluster rules dictionary from the cluster rules file"
    rulesdict = {}
    first = True
    cfg = config.get_config()
    
    for hmm_model in cfg.enabled_detection_models:
        dir_path = path.dirname(path.abspath(__file__))
        prefix = ""
        if hmm_model != "default":
            dir_path = path.join(dir_path, hmm_model)
            prefix = hmm_model + "/"
        #TODO: We should move all user-customizable files into config subdirectory; the rulefiles are redundant also in hmm_detection_dblookup
        for line in open(path.join(dir_path, "cluster_rules.txt"),"r"):
            # skip the first line with the legend
            if first:
                first = False
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            key = prefix + parts.pop(0)
            if key not in enabled_clustertypes:
                continue
            rules = parts.pop(0)
            cutoff = int(float(parts.pop(0)) * 1000.00 * cfg.cutoff_multiplier)
            extension = int(float(parts.pop(0)) * 1000.00 * cfg.cutoff_multiplier)
            rulesdict[key] = (rules, cutoff, extension)
    return rulesdict

def apply_cluster_rules(results_by_id, feature_by_id, enabled_clustertypes, rulesdict, overlaps, options):
    "Apply cluster rules to determine if HMMs lead to secondary metabolite core gene detection"
    typedict = {}
    cfg = config.get_config()
    cds_with_hits = sorted(list(results_by_id.keys()), key = lambda gene_id: feature_by_id[gene_id].location.start)
    for cds in cds_with_hits:
        _type = "none"
        #if typedict[cds] exist (the case of in-advance assignment from neighboring genes), use that instead of "none"
        if cds in list(typedict.keys()):
            _type = typedict[cds]
        cdsresults = [res.query_id for res in results_by_id[cds]]
        for clustertype in [ct for ct in enabled_clustertypes if ct not in _type.split("-")]:
            prefix = ""
            if len(clustertype.split("/")) > 1:
                prefix = clustertype.split("/")[0] + "/"
            single_rules = [rule for rule in rulesdict[clustertype][0].split(" or ") if " & " not in rule and "cluster(" not in rule]
            combined_rules = [rule[1:-1] for rule in rulesdict[clustertype][0].split(" or ") if " & " in rule]
            cluster_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "cluster(" in rule]
            minimum_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "minimum(" in rule]
            if "-" in clustertype:
                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
            else:
                cutoff = rulesdict[clustertype][1]
            base_cutoff = cutoff
            if cfg.enable_dynamic_cutoff:
                cutoff = int(base_cutoff * get_dynamic_cutoff_multiplier(cds, overlaps))
            #Assign cluster type if a single argument rule matches
            #Example rule format: "Domain1"
            if len(set([(prefix + rule) for rule in single_rules]) & set(cdsresults)) >= 1:
                if not (_type != "none" and clustertype == "other"):
                    if _type == "none" or _type == "other" or _type == clustertype:
                        _type = clustertype
                    elif clustertype not in _type.split("-"):
                        _type = clustertype + "-" + _type
                if _type != "other":
                    continue
            #Assign cluster type if a combinatorial argument rule matches
            #Example rule format: "(Domain1 & Domain2)"
            for rule in combined_rules:
                required_matches = [(prefix + rl) for rl in rule.split(" & ")]
                if len(set(required_matches) & set(cdsresults)) == len(required_matches):
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type.split("-"):
                            _type = clustertype + "-" + _type
            if _type == clustertype and _type != "other":
                continue
            #Assign cluster type if distance-based combinatorial parameter matches
            #Example rule format: "cluster(Domain1,Domain2)"
            for rule in cluster_rules:
                #Find which needed domains are already found in present CDS
                required_matches = [(prefix + rl) for rl in rule.split(",")]
                cluster_results = list(set(required_matches) & set(cdsresults))
                missing_results = list(set(required_matches) - set(cdsresults))
                #If more than one, search nearby CDSs for the other needed domains
                if len(cluster_results) > 0:
                    locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
                    for othercds in cds_with_hits:
                        needed_hits_found = list(set([res.query_id for res in results_by_id[othercds]]) & set(missing_results))
                        if len(needed_hits_found) > 0:
                            feature = feature_by_id[othercds]
                            flocations = [feature.location.start, feature.location.end]
                            #If hit found in nearby CDS, add it to the set of found domains relevant to the present rule
                            if min([abs(max(locations) - min(flocations)), abs(max(locations) - min(flocations))]) < cutoff:
                                cluster_results.extend(needed_hits_found)
                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
                if len(list(set(required_matches) - set(cluster_results))) == 0:
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type.split("-"):
                            _type = clustertype + "-" + _type
                        break
            #Assign cluster type if distance-based combinatorial parameter matches with minimum number of domain hits
            #Example rule format: "minimum(4,[Domain1,Domain2], [Domain1])"
            #This denotes that at least 4 hits of either Domain1 or Domain2 should be found in the same region with at least 1 hit from Domain 1
            for rule in minimum_rules:
                #Find which needed domains are already found in present CDS
                min_number = int(rule.partition(",")[0])
                required_matches = [(prefix + rl) for rl in rule.partition("[")[2].partition("]")[0].split(",")]
                essential_matches = [rl for rl in rule.rpartition("[")[2].partition("]")[0].split(",")]
                if essential_matches == ['']:
                    essential_matches = []
                cluster_results = list(set(required_matches) & set(cdsresults))
                missing_essential = list(set(essential_matches) - set(calc_essential_found(essential_matches, cdsresults, prefix)))
                neighborcds = [] #neighboring cds that together met the requirements
                nrcds = 0
                if len(cluster_results) > 0:
                    nrcds = 1
                #If essentials found but havent fulfilled, or required found but havent fulfilled, search nearby CDSs for the other needed domains
                locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
                nearest_cds = cds
                if ((len(list(set(essential_matches) & set(missing_essential))) > 0) or (nrcds > 0 and nrcds < min_number)):
                    #scan neighboring genes forward
                    for othercds in cds_with_hits:
                        if (overlaps[1][othercds] <= overlaps[1][cds]) or (overlaps[1][othercds] in [overlaps[1][ncds] for ncds in neighborcds]):
                            continue
                        feature = feature_by_id[othercds]
                        flocations = [feature.location.start, feature.location.end]
                        cluster_start = min(locations)
                        cluster_end = max(locations)
                        cluster_cds = list(set([cds]) | set(neighborcds))
                        if cfg.enable_dynamic_cutoff:
                            #in min_rules, use the nearest neighbor to calculate dynamic cutoff
                            if len(neighborcds) > 0:
                                nearest_cds = neighborcds[-1]
                            cutoff = max(int(base_cutoff * get_dynamic_cutoff_multiplier(othercds, overlaps)), int(base_cutoff * get_dynamic_cutoff_multiplier(nearest_cds, overlaps)))
                        within_cutoff = ((min(flocations) >= cluster_start and max(flocations) <= cluster_end)
                                        or (cluster_start - cutoff <= max(flocations) <= cluster_start)
                                        or (cluster_end <= min(flocations) <= cluster_end + cutoff))
                        within_gene_num_cutoff = (min([abs(overlaps[1][othercds] - overlaps[1][ncds]) for ncds in cluster_cds]) - 1 <= cfg.gene_num_cutoff)
                        if (cfg.gene_num_cutoff_only):
                            within_cutoff = within_gene_num_cutoff
                        else:
                            within_cutoff = within_cutoff or within_gene_num_cutoff
                        if (within_cutoff):
                            othercds_results = [res.query_id for res in results_by_id[othercds]]
                            #check essential hits
                            essential_found = calc_essential_found(missing_essential, othercds_results, prefix)
                            if len(essential_found) > 0:
                                missing_essential = list(set(missing_essential) - set(essential_found))
                                #append to neighborcds for in-advance tagging
                                if not (othercds in neighborcds):
                                    neighborcds.append(othercds)
                                    locations.extend(flocations)
                            #check required hits
                            needed_hits_found = list(set(othercds_results) & set(required_matches))
                            if len(needed_hits_found) > 0:
                                nrcds += 1
                                #append to neighborcds for in-advance tagging
                                if not (othercds in neighborcds):
                                    neighborcds.append(othercds)
                                    locations.extend(flocations)
                        else:
                            break
                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
                #for the cds and all the corresponding neighbor cdss
                if nrcds >= min_number and len(missing_essential) == 0:
                    if cfg.min_domain_number > 1:
                        #check for min. unique domains
                        unique_domains = list(set(cdsresults))
                        for ncds in neighborcds:
                            unique_domains = list(set(unique_domains) | set([res.query_id for res in results_by_id[ncds]]))
                        if len(unique_domains) < cfg.min_domain_number:
                            logging.debug('Cluster dropped due to not passing min. unique domains filtering [%s:%s]' % (feature_by_id[cds].location.start, feature_by_id[cds].location.end))
                            continue
                    if cfg.enable_cdhit:
                        #do cdhit filtering
                        clcds = [feature_by_id[key] for key in neighborcds]
                        clcds.append(feature_by_id[cds])
                        cdhit_table, gene_to_cluster = utils.get_cdhit_table(clcds, options)
                        detected_cdh_cluster = [gene_to_cluster[cds]]
                        for othercds in neighborcds:
                            if not (gene_to_cluster[othercds] in detected_cdh_cluster):
                                detected_cdh_cluster.append(gene_to_cluster[othercds])
                            else:
                                nrcds -= 1
                        if not nrcds >= min_number:
                            logging.debug('Cluster dropped due to not passing cd-hit filtering [%s:%s]' % (feature_by_id[cds].location.start, feature_by_id[cds].location.end))
                            continue
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type.split("-"):
                            _type = clustertype + "-" + _type
                    for ncds in neighborcds:
                        _ntype = "none"
                        if ncds in list(typedict.keys()):
                            _ntype = typedict[ncds]
                        if not (_ntype != "none" and clustertype == "other"):
                            if _ntype == "none" or _ntype == "other" or _ntype == clustertype:
                                _ntype = clustertype
                            elif clustertype not in _ntype.split("-"):
                                _ntype = clustertype + "-" + _ntype
                        typedict[ncds] = _ntype
                    break
        #Save type to typedict
        typedict[cds] = _type
    return typedict

def detect_signature_genes(seq_record, enabled_clustertypes, options):
    "Function to be executed by module"
    logging.info('Detecting gene clusters using HMM library')
    feature_by_id = utils.get_feature_dict(seq_record)
    rulesdict = create_rules_dict(enabled_clustertypes)
    results = []
    sig_by_name = {}
    results_by_id = {}
    for sig in get_sig_profiles():
        sig_by_name[sig.name] = sig

    for feature in utils.get_cds_features(seq_record):
        prefix = "%s:" % seq_record.id.replace(":", "_")
        gene_id = utils.get_gene_id(feature)
        if (prefix + gene_id) in options.hmm_results:
            results_by_id[gene_id] = options.hmm_results[prefix + gene_id]
            for res in results_by_id[gene_id]:
                results.append(res)

    short_cds_buffer = []
    if options.ignore_short_aa:
        # Temporarily filter out cds with < prot_min_length AA length
        min_length_aa = 50
        if options.eukaryotic:
            min_length_aa = 100
        for f in seq_record.features:
            if f.type == "CDS" and len(f.qualifiers['translation'][0]) < min_length_aa and utils.get_gene_id(f) not in results_by_id:
                short_cds_buffer.append(f)
                seq_record.features.remove(f)

    #Get overlap tables (for overlap filtering etc)
    overlaps = utils.get_overlaps_table(seq_record)

    #Filter results by comparing scores of different models (for PKS systems)
    results_to_delete = [gene_id for gene_id in results_by_id]
    results, results_by_id = filter_results(results, results_by_id, overlaps, feature_by_id)

    #Update filtered results back to the options.hmm_results
    for gene_id in results_by_id:
        results_to_delete.remove(gene_id)
        prefix = "%s:" % seq_record.id.replace(":", "_")
        if (prefix + gene_id) in options.hmm_results:
            options.hmm_results[(prefix + gene_id)] = results_by_id[gene_id]
    for gene_id in results_to_delete:
        prefix = "%s:" % seq_record.id.replace(":", "_")
        if (prefix + gene_id) in options.hmm_results:
            del options.hmm_results[(prefix + gene_id)]

    #Use rules to determine gene clusters
    typedict = apply_cluster_rules(results_by_id, feature_by_id, enabled_clustertypes, rulesdict, overlaps, options)

    #Rearrange hybrid clusters name in typedict alphabetically
    fix_hybrid_clusters_typedict(typedict)

    #Find number of sequences on which each pHMM is based
    nseqdict = get_nseq()

    #Save final results to seq_record
    for cds in list(results_by_id.keys()):
        feature = feature_by_id[cds]

        # Add domain record to feature no matter what
        result = "; ".join(["%s (E-value: %s, bitscore: %s, seeds: %s)" % (
        res.query_id, res.evalue, res.bitscore, nseqdict.get(res.query_id, '?')) for res in results_by_id[cds]])
        feature.qualifiers['domain_record'] = result

        if typedict[cds] != "none":
            _update_sec_met_entry(feature, results_by_id[cds], typedict[cds], nseqdict)


    find_clusters(seq_record, rulesdict, overlaps)


    #Rearrange hybrid clusters name alphabetically
    fix_hybrid_clusters(seq_record)

    #Add details of gene cluster detection to cluster features
    store_detection_details(results_by_id, rulesdict, seq_record)

    # Re-add the short CDSs
    seq_record.features.extend(short_cds_buffer)
    utils.sort_features(seq_record)

    #If all-orfs option on, remove irrelevant short orfs
    if options.all_orfs:
        remove_irrelevant_allorfs(seq_record)

    #Display %identity
    if options.enable_cdhit:
        store_percentage_identities(seq_record, options)


def get_nseq():
    nseqdict = {}
    for hmm in get_sig_profiles():
        hmmfile = hmm.hmm_file
        for line in open(hmmfile, 'r'):
            if line.startswith('NSEQ '):
                nseqdict[hmm.name] = line[6:].strip()
                break
        if not hmm.name in nseqdict:
            nseqdict[hmm.name] = "?"

    return nseqdict

def overlaps(feature1, feature2):
    if (feature2.location.start <= feature1.location.start <= feature2.location.end) or (feature2.location.start <= feature1.location.end <= feature2.location.end):
        return True
    else:
        return False


def get_dynamic_cutoff_multiplier(gene_id, overlaps, radius = 10):
    """Find the dynamic multiplier for cutoff based on radius number of nearest
    non overlapping gene numbers density"""
    result = 1.00

    center_idx = overlaps[1][gene_id]
    left_idx = center_idx
    right_idx = center_idx

    while (right_idx - left_idx) <= radius:
        if left_idx <= 0 and right_idx >= (len(overlaps[2]) - 1):
            break
        if left_idx <= 0:
            right_idx += 1
            continue
        if right_idx >= (len(overlaps[2]) - 1):
            left_idx -= 1
            continue
        if (overlaps[2][center_idx][0] - overlaps[2][left_idx - 1][1]) > (overlaps[2][right_idx + 1][0] - overlaps[2][center_idx][1]):
            right_idx += 1
        else:
            left_idx -= 1

    result = float(overlaps[2][right_idx][1] - overlaps[2][left_idx][0]) / float(right_idx - left_idx + 1)
    result /= 1000.00

    return result


def calc_intergenic_distance(feature_1, feature_2):
    "calculate intergenic distance between 2 features, return 0 if overlapped"
    result = 0

    f1_loc = [feature_1.location.start, feature_1.location.end]
    f1_loc.sort()
    f2_loc = [feature_2.location.start, feature_2.location.end]
    f2_loc.sort()

    if (f1_loc[0] > f2_loc[1]):
        #f1 is downstream of f2
        result = f1_loc[0] - f2_loc[1]
    elif (f2_loc[0] > f1_loc[1]):
        #f2 is downstream of f1
        result = f2_loc[0] - f1_loc[1]

    return result


def calc_essential_found(essential_matches, cdsresults, prefix = None):
    "Given array of essential matches, returns array of essential matches forfilled by the cdsresults"
    essential_found = []
    for match in essential_matches:
        for domain in match.split("/"):
            if not prefix is None:
                domain = prefix + domain
            if domain in cdsresults:
                essential_found.append(match)
                break
    return essential_found


def remove_irrelevant_allorfs(seq_record):
    #Get features
    allfeatures = utils.get_cds_features(seq_record)
    #Remove auto-orf features without unique sec_met qualifiers; remove glimmer ORFs overlapping with sec_met auto-orfs not catched by Glimmer
    auto_orf_features = [feature for feature in allfeatures if 'note' in feature.qualifiers and "auto-all-orf" in feature.qualifiers['note']]
    other_features = [feature for feature in allfeatures if 'note' not in feature.qualifiers or "auto-all-orf" not in feature.qualifiers['note']]
    to_delete = []
    for autofeature in auto_orf_features:
        if "sec_met" not in autofeature.qualifiers:
            to_delete.append(autofeature)
        else:
            glimmer_has_sec_met = False
            for otherfeature in other_features:
                if overlaps(autofeature, otherfeature) and 'sec_met' in otherfeature.qualifiers:
                    to_delete.append(autofeature)
                    glimmer_has_sec_met = True
            if glimmer_has_sec_met == False:
                for otherfeature in other_features:
                    if overlaps(autofeature, otherfeature) and 'sec_met' not in otherfeature.qualifiers:
                        to_delete.append(otherfeature)
    featurenrs = []
    idx = 0
    for feature in seq_record.features:
        if feature in to_delete:
            featurenrs.append(idx)
        idx += 1
    featurenrs.reverse()
    for featurenr in featurenrs:
        del seq_record.features[featurenr]


def store_percentage_identities(seq_record, options):
    clusters = utils.get_cluster_features(seq_record)
    cfg = config.get_config()
    for cluster in clusters:
        features = [feature for feature in utils.get_cluster_cds_features(cluster, seq_record) if 'sec_met' in feature.qualifiers]
        cdhit_table, gene_to_cluster = utils.get_cdhit_table(features, options, float(cfg.cdh_display_cutoff))
        for cdhit_cluster in cdhit_table:
            if len(cdhit_cluster["genes"]) > 1:
                cl_features = [feature for feature in features if utils.get_gene_id(feature) in list(cdhit_cluster["genes"].keys())]
                pct_table = utils.get_pct_identity_table(cl_features)
                for cds in cl_features:
                    result = ",".join(["%s=%s" % (othercds, pct_table[utils.get_gene_id(cds)][othercds]) for othercds in list(pct_table[utils.get_gene_id(cds)].keys())])
                    for ann in cds.qualifiers['sec_met']:
                        if ann.startswith("Percentage identity"):
                            del ann
                    cds.qualifiers['sec_met'].append("Percentage identity: %s" % (result))


def fix_hybrid_clusters(seq_record):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        clustertypes = cluster.qualifiers['product'][0].split("-")
        clustertypes.sort()
        cluster.qualifiers['product'][0] = "-".join(clustertypes)


def fix_hybrid_clusters_typedict(typedict):
    for key in list(typedict.keys()):
        clustertypes = typedict[key].split("-")
        clustertypes.sort()
        typedict[key] = "-".join(clustertypes)


def store_detection_details(results_by_id, rulesdict, seq_record):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        type_combo = utils.get_cluster_type(cluster)
        if '-' in type_combo:
            clustertypes = type_combo.split('-')
        else:
            clustertypes = [type_combo]

        if not 'note' in cluster.qualifiers:
            cluster.qualifiers['note'] = []
        rule_string = "Detection rule(s) for this cluster type:"
        for clustertype in clustertypes:
            rule_string += " %s: (%s);" % (clustertype, rulesdict[clustertype][0])

        cluster.qualifiers['note'].append(rule_string)

def _update_sec_met_entry(feature, results, clustertype, nseqdict):
    result = "; ".join(["%s (E-value: %s, bitscore: %s, seeds: %s)" % (res.query_id, res.evalue, res.bitscore, nseqdict.get(res.query_id, '?'))  for res in results])

    if not 'sec_met' in feature.qualifiers:
        feature.qualifiers['sec_met'] = [
            "Type: %s" % clustertype,
            "Domains detected: %s" % (result),
            "Kind: biosynthetic"
        ]
    else:
        for ann in feature.qualifiers['sec_met']:
            if not ann.startswith("Domains detected"):
                continue
            ann += "Domains detected: %s" % (result)

__all__ = [ check_prereqs, detect_signature_genes]
