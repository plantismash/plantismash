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
import logging
import os
import sys
import shutil
from os import path
import subprocess
import string
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
from argparse import Namespace
import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO
from Bio.Alphabet import generic_protein
from antismash.config import get_config
from helperlibs.wrappers.io import TemporaryDirectory

from zipfile import ZipFile, ZIP_DEFLATED, LargeZipFile
from Bio.Seq import Seq
import re
try:
    from antismash.db import biosql
    USE_BIOSQL = True
except ImportError:
    USE_BIOSQL = False
import argparse

class Storage(dict):
    """Simple storage class"""
    def __init__(self, indict=None):
        if indict is None:
            indict = {}
        dict.__init__(self, indict)
        self.__initialized = True

    def __getattr__(self, attr):
        try:
            return self.__getitem__(attr)
        except KeyError:
            raise AttributeError(attr)

    def __setattr__(self, attr, value):
        if '_Storage__initialized' not in self.__dict__:
            return dict.__setattr__(self, attr, value)
        elif attr in self:
            dict.__setattr__(self, attr, value)
        else:
            self.__setitem__(attr, value)


def getArgParser():
    class AntiSmashParser(argparse.ArgumentParser):
        """Custom argument parser for antiSMASH
        """
        _showAll = False
        _displayGroup = {}

        def __init__(self, *args):
            """Initialisation method for the parser class"""
            kwargs = {}
            kwargs["add_help"] = False
            super(AntiSmashParser, self).__init__(*args, **kwargs)

        def add_argument_group(self, *args, **kwargs):
            if not args[0] in self._displayGroup:
                self._displayGroup[args[0]] = []
            if "basic" in kwargs:
                if kwargs["basic"]:
                    self._displayGroup[args[0]].extend(["basic"])
                del kwargs["basic"]
            if "param" in kwargs:
                self._displayGroup[args[0]].extend(kwargs["param"])
                del kwargs["param"]
            group = super(AntiSmashParser, self).add_argument_group(*args, **kwargs)
            return group

        def print_help(self, file=None, showAll=False):
            self._showAll = showAll
            super(AntiSmashParser, self).print_help(file)

        def format_help(self):
            """Custom help format"""
            help_text = """
########### antiSMASH ver. {version} #############

{usage}

{args}
--------
Options
--------
{opts}
""".format(version=get_version(), usage=self.format_usage(), args=self._get_args_text(), opts=self._get_opts_text())
            return help_text

        def format_usage(self):
            if self._showAll:
                formatter = self._get_formatter()
                formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
                return formatter.format_help()
            return "usage: {prog} [-h] [options ..] [sequence [sequence ..]]".format(prog=self.prog) + "\n"

        def _get_args_text(self):
            # fetch arg lists using formatter
            formatter = self._get_formatter()
            for action_group in self._action_groups:
                if action_group.title == "positional arguments":
                    formatter.start_section("arguments")
                    formatter.add_arguments(action_group._group_actions)
                    formatter.end_section()
                    break
            return formatter.format_help()

        def _get_opts_text(self):
            # fetch opt lists using formatter
            formatter = self._get_formatter()
            for action_group in self._action_groups:
                if action_group.title in ["optional arguments"]:
                    formatter.add_arguments(action_group._group_actions)
            for action_group in self._action_groups:
                if action_group.title not in ["optional arguments", "positional arguments"]:
                    show_opt = self._showAll
                    if (not show_opt):
                        if "basic" in self._displayGroup[action_group.title]:
                            show_opt = True
                        elif len(list(set(sys.argv) & set(self._displayGroup[action_group.title]))) > 0:
                            show_opt = True
                    if show_opt:
                        formatter.start_section(action_group.title)
                        if action_group.description is None:
                            action_group.description = ''
                        formatter.add_text(action_group.description)
                        formatter.add_arguments(action_group._group_actions)
                        formatter.end_section()
            return formatter.format_help()

    return AntiSmashParser()


def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features

def get_all_features_of_type_with_query(seq_record, feature_type, query_tag, query_value):
    """Return all features of type 'type' which contain a 'query_tag' with value 'query_value'
    Note: query has to be exact!"""

    features_with_type = get_all_features_of_type(seq_record, feature_type)
    features = []
    for feature_to_test in features_with_type:

        if (query_tag in feature_to_test.qualifiers) and (query_value in feature_to_test.qualifiers[query_tag]):
            features.append(feature_to_test)
    return features

def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")

def get_withincluster_cds_features(seq_record):
    features = get_cds_features(seq_record)
    clusters = get_cluster_features(seq_record)
    withinclusterfeatures = []
    for feature in features:
        for cluster in clusters:
            if not (cluster.location.start <= feature.location.start <= cluster.location.end or \
               cluster.location.start <= feature.location.end <= cluster.location.end):
                continue
            if feature not in withinclusterfeatures:
                withinclusterfeatures.append(feature)
    return withinclusterfeatures

def get_cluster_cds_features(cluster, seq_record):
    clustercdsfeatures = []
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            clustercdsfeatures.append(feature)
    return clustercdsfeatures

def get_cluster_aSDomain_features(cluster, seq_record):
    aSDomainfeatures = []
    for feature in seq_record.features:
        if feature.type != 'aSDomain':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            aSDomainfeatures.append(feature)
    return aSDomainfeatures


def get_cluster_domains(cluster, seq_record):
    "Get unique domains list of a gene cluster"
    result = []
    features = get_cluster_cds_features(cluster, seq_record)
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for note in feature.qualifiers.get('sec_met', []):
                if not note.startswith('Domains detected'):
                    continue
                note = note[18:]
                for domain in note.split(';'):
                    domain = domain.split(" (")[0]
                    if not domain in result:
                        result.append(domain)
    result.sort()
    return result

def sort_substrates(substrates, type):
    # Count occurrences
    count_dict = {}
    for substrate in substrates:
        if substrate in count_dict:
            count_dict[substrate] += 1
        else:
            count_dict[substrate] = 1

    # Sort by count (descending) then alphabetically
    sorted_substrates = sorted(count_dict.items(), key=lambda x: (-x[1], x[0]))

    # Prioritize elements based on conditions
    prioritized = []
    others = []
    for substrate, count in sorted_substrates:
        if ("terpenoid"  in substrate.lower() or "oleananes" in substrate.lower()) and "terpene" in type:
            substrate = "*" + substrate
            if count > 1:
                prioritized.append("{}-{}".format(substrate, count))
            else:
                prioritized.append(substrate)
        elif substrate in type:
            substrate = "*" + substrate
            if count > 1:
                prioritized.append("{}-{}".format(substrate, count))
            else:
                prioritized.append(substrate)
        else:
            if count > 1:
                others.append("{}-{}".format(substrate, count))
            else:
                others.append(substrate)

    # Combine prioritized and other elements
    sorted_substrates_times = prioritized + others

    return sorted_substrates_times if sorted_substrates_times else ["-"]


def get_cluster_substrates(cluster, seq_record):
    "Get unique predicted substrates type list of a gene cluster"
    result = []
    substrates = []
    type = get_cluster_type(cluster).split("-")
    features = get_cluster_cds_features(cluster, seq_record)
    for feature in features:
        if 'substrates' in feature.qualifiers:
            substrates += feature.qualifiers['substrates']

    return sort_substrates(substrates, type)


def features_overlap(a, b):
    "Check if two features have overlapping locations"
    astart = min([a.location.start, a.location.end])
    aend = max([a.location.start, a.location.end])
    bstart = min([b.location.start, b.location.end])
    bend = max([b.location.start, b.location.end])
    return (astart >= bstart and astart <= bend) or \
           (aend >= bstart and aend <= bend) or \
           (bstart >= astart and bstart <= aend) or \
           (bend >= astart and bend <= aend)

def get_secmet_cds_features(seq_record):
    features = get_cds_features(seq_record)
    secmet_features = []
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for annotation in feature.qualifiers['sec_met']:
                if annotation.startswith('Type: '):
                    secmet_features.append(feature)
                    break
    return secmet_features

def get_pksnrps_cds_features(seq_record):
    features = get_cds_features(seq_record)
    pksnrpscoregenes = []
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for annotation in feature.qualifiers['sec_met']:
                if annotation.startswith('NRPS/PKS Domain:'):
                    pksnrpscoregenes.append(feature)
                    break
    return pksnrpscoregenes

def get_nrpspks_domain_dict(seq_record):
    domaindict = {}
    features = get_cds_features(seq_record)
    for feature in features:
        domainlist = []
        if 'sec_met' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                hit_id =  domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[0]
                domstart = domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[2].partition("-")[0]
                domend = domain.partition("NRPS/PKS Domain: ")[2].partition("). ")[0].rpartition("-")[2]
                evalue = domain.partition("E-value: ")[2].partition(". Score:")[0]
                bitscore = domain.partition("Score: ")[2].partition(";")[0]
                domainlist.append([hit_id, int(domstart), int(domend), evalue, float(bitscore)])
            if len(domainlist) > 0:
                domaindict[get_gene_id(feature)] = domainlist
    return domaindict

def get_nrpspks_substr_spec_preds(seq_record):
    substr_spec_preds = Storage()
    substr_spec_preds.consensuspreds = {}
    substr_spec_preds.nrps_svm_preds = {}
    substr_spec_preds.nrps_code_preds = {}
    substr_spec_preds.minowa_nrps_preds = {}
    substr_spec_preds.pks_code_preds = {}
    substr_spec_preds.minowa_pks_preds = {}
    substr_spec_preds.minowa_cal_preds = {}
    substr_spec_preds.kr_activity_preds = {}
    substr_spec_preds.kr_stereo_preds = {}
    features = get_cds_features(seq_record)
    for feature in features:
        nrat, nra, nrcal, nrkr = 0, 0, 0, 0
        if 'sec_met' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                if "AMP-binding" in domain or "A-OX" in domain:
                    nra += 1
                    domainname = get_gene_id(feature) + "_A" + str(nra)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    nrps_svm_pred = predictionstext.partition(" (NRPSPredictor2 SVM)")[0]
                    nrps_code_pred = predictionstext.partition(" (NRPSPredictor2 SVM), ")[2].partition(" (Stachelhaus code)")[0]
                    minowa_nrps_pred = predictionstext.partition("(Stachelhaus code), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.nrps_svm_preds[domainname] = nrps_svm_pred
                    substr_spec_preds.nrps_code_preds[domainname] = nrps_code_pred
                    substr_spec_preds.minowa_nrps_preds[domainname] = minowa_nrps_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_AT" in domain:
                    nrat += 1
                    domainname = get_gene_id(feature) + "_AT" + str(nrat)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    pks_code_pred = predictionstext.partition(" (PKS signature)")[0]
                    minowa_pks_pred = predictionstext.partition("(PKS signature), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.pks_code_preds[domainname] = pks_code_pred
                    substr_spec_preds.minowa_pks_preds[domainname] = minowa_pks_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "CAL_domain" in domain:
                    nrcal += 1
                    domainname = get_gene_id(feature) + "_CAL" + str(nrcal)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    minowa_cal_pred = predictionstext.partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.minowa_cal_preds[domainname] = minowa_cal_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_KR" in domain:
                    nrkr += 1
                    domainname = get_gene_id(feature) + "_KR" + str(nrkr)
                    activityprediction = domain.partition("Predicted KR activity: ")[2].partition(";")[0]
                    stereoprediction = domain.partition("Predicted KR stereochemistry: ")[2].partition(";")[0]
                    substr_spec_preds.kr_activity_preds[domainname] = activityprediction
                    substr_spec_preds.kr_stereo_preds[domainname] = stereoprediction
    return substr_spec_preds

def get_smcog_annotations(seq_record):
    smcogdict = {}
    smcogdescriptions = {}
    features = get_cds_features(seq_record)
    for feature in features:
        if 'note' in feature.qualifiers:
            notes = feature.qualifiers['note']
            for note in notes:
                if "smCOG: " in note:
                    smcogid = note.partition("smCOG: ")[2].partition(":")[0]
                    smcog_descr = note.partition("smCOG: ")[2].partition(":")[2].partition("(Score:")[0]
                    smcogdict[get_gene_id(feature)] = smcogid
                    smcogdescriptions[smcogid] = smcog_descr
    return smcogdict, smcogdescriptions

def get_pfam_features(seq_record):
    "Return all CDS_motif features containing a 'PFAM-Id: ' note for a seq_record"
    pfam_features = []
    for feature in get_all_features_of_type(seq_record, "PFAM_domain"):
        if 'db_xref' not in feature.qualifiers:
            continue
        for xref in feature.qualifiers['db_xref']:
            if xref.startswith("PFAM: "):
                pfam_features.append(feature)
                break
    return pfam_features

def get_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    return get_all_features_of_type(seq_record, "cluster")

def get_sorted_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    clusters = get_all_features_of_type(seq_record, "cluster")
    if len(clusters) == 0:
        return []
    numberdict = {}
    for cluster in clusters:
        numberdict[get_cluster_number(cluster)] = cluster
    return [numberdict[clusternr] for clusternr in numberdict.keys()]

def get_structure_pred(cluster):
    "Return all structure prediction for a cluster feature"
    if 'note' in cluster.qualifiers:
        for note in cluster.qualifiers['note']:
            if "Monomers prediction: " in note:
                return note.partition("Monomers prediction: ")[2]
    if get_cluster_type(cluster) == 'ectoine':
        return 'ectoine'
    return "N/A"

def get_cluster_number(cluster):
    "Get the integer representation of the Cluster number qualifier"
    if 'note' in cluster.qualifiers and "Cluster number: " in cluster.qualifiers['note'][0]:
        clusternr = int(cluster.qualifiers['note'][0].partition("Cluster number: ")[2])
        return clusternr
    else:
        return 0

def get_cluster_type(cluster):
    "Get product type of a gene cluster"
    return cluster.qualifiers['product'][0]

def get_cluster_by_nr(seq_record, queryclusternr):
    "Return all cluster features for a seq_record of a certain type"
    clusters = get_all_features_of_type(seq_record, "cluster")
    for cluster in clusters:
        if "Cluster number: " in cluster.qualifiers['note'][0]:
            clusternr = int(cluster.qualifiers['note'][0].partition("Cluster number: ")[2])
            if clusternr == queryclusternr:
                return cluster

def get_cluster_features_of_type(seq_record, clustertype):
    "Return all cluster features for a seq_record of a certain type"
    clusters = get_all_features_of_type(seq_record, "cluster")
    return [cluster for cluster in clusters if clustertype in cluster.qualifiers['product'][0]]

def locate_executable(name):
    "Find an executable in the path and return the full path"
    # In windows, executables tend to end on .exe
    if sys.platform == 'win32':
        name += ".exe"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for p in os.environ["PATH"].split(os.pathsep):
        full_name = path.join(p, name)
        if path.isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None

def locate_file(name):
    "Find a file and return the full path"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.R_OK):
            logging.debug("Found file %r", name)
            return name
    return None

# Ignore the pylint warning about input being redifined, as we're just
# following the subprocess names here.
# pylint: disable=redefined-builtin
def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
        raise
# pylint: enable=redefined-builtin


def run_hmmsearch(query_hmmfile, target_sequence, cutoff = None):
    "Run hmmsearch"
    config = get_config()
    command = ["hmmsearch", "--cpu", str(config.cpus),
               query_hmmfile, '-']
    if cutoff != None:
        if cutoff < 0:
            tc_exist = False
            ga_exist = False
            nc_exist = False
            for line in open(query_hmmfile,"r"):
                firstWord = line.split()[0]
                if firstWord == "TC":
                    tc_exist = True
                if firstWord == "GA":
                    ga_exist = True
                if firstWord == "NC":
                    nc_exist = True
            if tc_exist:
                command.insert(1, "--cut_tc")
            elif ga_exist:
                command.insert(1, "--cut_ga")
            elif nc_exist:
                command.insert(1, "--cut_nc")
            else:
                cutoff = 20
        if cutoff > 0:
            command.insert(1, str(cutoff))
            command.insert(1, "-T")
            command.insert(1, str(cutoff))
            command.insert(1, "--domT")
    try:
        out, err, retcode = execute(command, input=target_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmsearch returned %d: %r while searching %r', retcode,
                        err, query_hmmfile)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results


def run_hmmscan(target_hmmfile, query_sequence, opts=None, query_sequence_path=None):
    "Run hmmscan"
    config = get_config()
    command = ["hmmscan", "--cpu", str(config.cpus), "--nobias"]
    if opts is not None:
        command.extend(opts)
    if query_sequence_path:
        command.extend([target_hmmfile, query_sequence])
        try:
            out, err, retcode = execute(command)
        except OSError:
            return []
    else:
        command.extend([target_hmmfile, '-'])
        try:
            out, err, retcode = execute(command, input=query_sequence)
        except OSError:
            return []
    if retcode != 0:
        logging.debug('hmmscan returned %d: %r while scanning %r' , retcode,
                        err, query_sequence)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results


def run_hmmpfam2(query_hmmfile, target_sequence):
    "Run hmmpfam2"
    config = get_config()
    command = ["hmmpfam2", "--cpu", str(config.cpus),
               query_hmmfile, '-']
    try:
        out, err, retcode = execute(command, input=target_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmpfam2 returned %d: %r while searching %r', retcode,
                        err, query_hmmfile)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer2-text'))
    return results


def run_hmmpress(hmmfile):
    "Run hmmpress"
    command = ['hmmpress', hmmfile]
    try:
        out, err, retcode = execute(command)
    except OSError as e:
        retcode = 1
        err = str(e)
    return out, err,  retcode


def hmmlengths(hmmfile):
    hmmlengthsdict = {}
    openedhmmfile = open(hmmfile,"r")
    filetext = openedhmmfile.read()
    filetext = filetext.replace("\r","\n")
    hmms = filetext.split("//")[:-1]
    for i in hmms:
        namepart = i.split("NAME  ")[1]
        name = namepart.split("\n")[0]
        lengthpart = i.split("LENG  ")[1]
        length = lengthpart.split("\n")[0]
        hmmlengthsdict[name] = int(length)
    return hmmlengthsdict

def cmp_feature_location(a, b):
    "Compare two features by their start/end locations"
    ret = cmp(a.location.start, b.location.start)
    if ret != 0:
        return ret
    return cmp(a.location.end, b.location.end)

def sort_features(seq_record):
    "Sort features in a seq_record by their position"
    #Check if all features have a proper location assigned
    for feature in seq_record.features:
        if feature.location is None:
            if feature.id != "<unknown id>":
                logging.error("Feature '%s' has no proper location assigned", feature.id)
            elif "locus_tag" in feature.qualifiers:
                logging.error("Feature '%s' has no proper location assigned", feature.qualifiers["locus_tag"][0])
            else:
                logging.error("File contains feature without proper location assignment")
            sys.exit(0) #FIXME: is sys.exit(0) really what we want to do here?
    #Sort features by location
    seq_record.features.sort(cmp=cmp_feature_location)

def fix_locus_tags(seq_record, config):
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    if 'next_locus_tag' not in config:
        config.next_locus_tag = 1

    cds_list = get_cds_features(seq_record)
    for feature in cds_list:
        if get_gene_id(feature) == "no_tag_found":
            feature.qualifiers['locus_tag'] = ['AUTOORF_%05d' % config.next_locus_tag]
            config.next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars  = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        if 'locus_tag' in feature.qualifiers:
            for char in feature.qualifiers['locus_tag'][0]:
                if char in illegal_chars:
                    feature.qualifiers['locus_tag'][0] = feature.qualifiers['locus_tag'][0].replace(char, "_")
        if 'gene' in feature.qualifiers:
            for char in feature.qualifiers['gene'][0]:
                if char in illegal_chars:
                    feature.qualifiers['gene'][0] = feature.qualifiers['gene'][0].replace(char, "_")
        if 'protein_id' in feature.qualifiers:
            for char in feature.qualifiers['protein_id'][0]:
                if char in illegal_chars:
                    feature.qualifiers['protein_id'][0] = feature.qualifiers['protein_id'][0].replace(char, "_")

def _shorten_ids(idstring, options):
    contigstrmatch = re.search(r"onti?g?(\d+)\b", idstring)
    if not contigstrmatch:
        # if there is a substring "[Ss]caf(fold)XXX" use this number
        contigstrmatch = re.search(r"caff?o?l?d?(\d+)\b", idstring)
    if not contigstrmatch:
        # if there is a substring " cXXX" use this number
        contigstrmatch = re.search(r"\bc(\d+)\b", idstring)
    if contigstrmatch:
        contig_no = int(contigstrmatch.group(1))
    else:
        # if the contig number cannot be parsed out, just count the contigs from 1 to n
        contig_no = options.orig_record_idx

    return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])


def fix_record_name_id(seq_record, options):
    "Fix a seq record's name and id to be <= 16 characters, the GenBank limit; if record name is too long, add c000X prefix"

    if seq_record.id == "unknown.1":
        seq_record.id = "unk_seq_{ctg:05d}".format(ctg=options.orig_record_idx)
        logging.warn('Invalid sequence id "unknown.1", replaced by %s', seq_record.id)

    if seq_record.name == "unknown":
        seq_record.name = "unk_seq_{ctg:05d}".format(ctg=options.orig_record_idx)
        logging.warn('Invalid sequence name "unknown", replaced by %s', seq_record.name)

    if len(seq_record.id.partition(".")[0]) > 16:
        oldid = seq_record.id

        #Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that is just too long because of the version number behind the dot
        if (seq_record.id[-2] == "." and
                seq_record.id.count(".") == 1 and
                len(seq_record.id.partition(".")[0]) <= 16 and
                seq_record.id.partition(".")[0] not in options.all_record_ids):
            seq_record.id = seq_record.id.partition(".")[0]
            options.all_record_ids.append(seq_record.id)
        else: #Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(oldid, options) not in options.all_record_ids:
                seq_record.id = _shorten_ids(oldid, options)
                options.all_record_ids.append(seq_record.id)
            else:
                x = 0
                while "%s_%i" % (seq_record.id[:16][:-4], x) in options.all_record_ids:
                    x += 1
                seq_record.id = "%s_%i" % (seq_record.id[:16][:-4], x)
                options.all_record_ids.append(seq_record.id)

        logging.warn('Fasta header too long: renamed "%s" to "%s"', oldid, seq_record.id)
        if seq_record.id not in options.extrarecord:
            options.extrarecord[seq_record.id] = Namespace()
        if "extradata" not in options.extrarecord[seq_record.id]:
            options.extrarecord[seq_record.id].extradata = {}
        if "orig_id" not in options.extrarecord[seq_record.id].extradata:
            options.extrarecord[seq_record.id].extradata["orig_id"] = oldid
        #seq_record.id = "%s...%s" % (seq_record.id[:12], seq_record.id[-1])


    if len(seq_record.name) > 16:

        seq_record.name = _shorten_ids(seq_record.name, options)

        #seq_record.name = "%s...%s" % (seq_record.name[:12], seq_record.name[-1])

    if 'accession' in seq_record.annotations and \
       len(seq_record.annotations['accession']) > 16:
        acc = seq_record.annotations['accession']

        seq_record.annotations['accession'] = _shorten_ids(acc, options)
        # seq_record.annotations['accession'] = "%s...%s" % (acc[:12], acc[-1])

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|}/ '''
    for char in seq_record.id:
        if char in illegal_chars:
            seq_record.id = seq_record.id.replace(char,"")
    for char in seq_record.name:
        if char in illegal_chars:
            seq_record.id = seq_record.name.replace(char,"")

def ascii_string(inputstring):
    return "".join([char for char in inputstring if char in (string.ascii_letters + string.digits + string.punctuation + string.whitespace)])

def get_gene_acc(feature):
    "Get the gene accesion number; if not defined, use content of locus_tag or gene qualifier instead"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]

    # I needed to include this for clusterblast to work as well with the data from db (which all has protein_id) as with user data (which hasn't)
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_accession"

def get_gene_accs(feature):
    "Get array of the protein_id tags"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id']
    elif 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag']
    elif 'gene' in feature.qualifiers:
        return feature.qualifiers['gene']
    return ["no_accession"]

def get_ncbi_gi(feature):
    """Get the NCBI gi from feature
    returns gi if found, None if not present"""

    gi = None
    if 'db_xref' in feature.qualifiers:
        for db_xref in feature.qualifiers['db_xref']:
            if 'GI:' in db_xref:
                gi = db_xref[3:]
    return gi

def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"

def get_gene_annotation(feature):
    "Get the gene annotation from the product qualifier"
    if 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    return "unannotated orf"

def get_gene_accession(feature):
    "Get the gene accession number from protein_id, locus_tag or gene name, in that order"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_tag_found"

def get_full_path(current_file, file_to_add):
    "Get the full path of file_to_add in the same directory as current_file"
    return path.join(path.dirname(path.abspath(current_file)), file_to_add)

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id

def get_feature_dict_protein_id(seq_record):
    """Get a dictionary mapping features to their protein_ids"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_ids = get_gene_accs(feature)

        # As the HMMer result file does only contain result lines for 1 protein_id;
        # if an entry has 2 (or more) protein_id entries, lookups failed in the original version
        # Therefore I changed this method to use the protein_id list instead of just
        # extracting the [0] element - althouth we loose some memory here due to redundancy.
        for gene_id in gene_ids:
            feature_by_id[gene_id] = feature
    return feature_by_id

def get_multifasta(seq_record):
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = get_cds_features(seq_record)
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers['translation'][0]
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping", gene_id)
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta


def get_overlaps_table(seq_record):
    """Given seq_record, returns array of overlapping genes
    and the corresponding gene_id-to-indexes"""
    overlaps = []
    overlap_by_id = {}
    overlaps_loc = []
    features = get_cds_features(seq_record)
    if len(features) < 1:
        return overlaps, overlap_by_id, overlaps_loc
    features.sort(key=lambda feature: feature.location.start)
    i = 0
    j = i + 1
    cds_queue = []
    while i < len(features):
        if j >= len(features):
            break
        cds = features[i]
        ncds = features[j]
        if (cds.location.end <= ncds.location.start + 1):
            overlaps.append([])
            overlaps_loc.append((cds.location.start, cds.location.end))
            cds_queue.append(cds)
            for cds in cds_queue:
                overlap_by_id[get_gene_id(cds)] = len(overlaps) - 1
                overlaps[-1].append(cds)
                overlaps_loc[-1] = (min(overlaps_loc[-1][0], cds.location.start), max(overlaps_loc[-1][1], cds.location.end))
            cds_queue = []
            i = j
        else:
            if (cds.location.end < ncds.location.end):
                cds_queue.append(cds)
                i = j
            else:
                cds_queue.append(ncds)
        j += 1
    overlaps.append([])
    overlaps_loc.append((features[i].location.start, features[i].location.end))
    cds_queue.append(features[i])
    for cds in cds_queue:
        overlap_by_id[get_gene_id(cds)] = len(overlaps) - 1
        overlaps[-1].append(cds)
        overlaps_loc[-1] = (min(overlaps_loc[-1][0], cds.location.start), max(overlaps_loc[-1][1], cds.location.end))
    return overlaps, overlap_by_id, overlaps_loc


def get_specific_multifasta(features):
    """Extract multi-protein FASTA from provided features"""
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers['translation'][0]

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping", gene_id)
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def get_aa_translation(seq_record, feature):
    """Obtain content for translation qualifier for specific CDS feature in sequence record"""
    fasta_seq = feature.extract(seq_record.seq).ungap('-').translate(to_stop=True)
    if len(fasta_seq) == 0:
        logging.debug("Retranslating %s with stop codons", feature.id)
        fasta_seq = feature.extract(seq_record.seq).ungap('-').translate()
    if "*" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("*","X"), generic_protein)
    if "-" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)

    return fasta_seq

def get_aa_sequence(feature, to_stop=False):
    """Extract sequence from specific CDS feature in sequence record"""
    fasta_seq = feature.qualifiers['translation'][0]
    if "*" in fasta_seq:
        if to_stop:
            fasta_seq = fasta_seq.split('*')[0]
        else:
            fasta_seq = fasta_seq.replace("*","X")
    if "-" in fasta_seq:
        fasta_seq = fasta_seq.replace("-","")
    return fasta_seq

def writefasta(names, seqs, filename):
    "Write sequence to a file"
    try:
        with open(filename, "w") as out_file:
            logging.info("Writing to FASTA file {} with {} sequences.".format(filename, len(names)))
            for i in range(len(names)):
                out_file.write(">{}\n{}\n".format(names[i], seqs[i]))
        logging.info("FASTA file {} successfully written and closed.".format(filename))
    except Exception as e:
        logging.error("Failed to write to FASTA file {}. Error: {}".format(filename, e))


def sortdictkeysbyvaluesrev(indict):
    items = [(value, key) for key, value in indict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]

def sortdictkeysbyvaluesrevv(indict):
    values = indict.values()
    values.sort()
    values.reverse()
    return values


def get_all_geo_rec_ids(seq_records):
    "List all GEO rec ids present in the seq_records"
    result = []
    for seq_record in seq_records:
        for feature in get_cds_features(seq_record):
            if "geo" in feature.qualifiers:
                for line in feature.qualifiers["geo"]:
                    reg = re.compile("(Record: )(?P<rec_id>.+)(;)(?P<ref_id>.+)(\()(?P<evalue>[0-9eE\-\.]+)(\);)(?P<data>.*)(;)(?P<dist>.*)")
                    match = reg.search(line)
                    if match:
                        rec_id = match.group("rec_id")
                        if rec_id not in result:
                            result.append(rec_id)
    return result


def parse_geo_feature(feature):
    "Parse 'geo' annotation from feature qualifiers"
    result = []

    if "geo" in feature.qualifiers:
        for line in feature.qualifiers["geo"]:
            reg = re.compile("(Record: )(?P<rec_id>.+)(;)(?P<ref_id>.+)(\()(?P<evalue>[0-9eE\-\.]+)(\);)(?P<data>.*)(;)(?P<dist>.*)")
            match = reg.search(line)
            if match:
                rec_id = match.group("rec_id")
                ref_id = match.group("ref_id")
                evalue = float(match.group("evalue"))
                data = {}
                reg = re.compile("(?P<sample_id>.+)(=)(?P<perc_val>[0-9]+)(\()(?P<real_val>[0-9eE\-\.]+)(\))")
                for d in match.group("data").split(","):
                    dm = reg.search(d)
                    if dm:
                        data[dm.group("sample_id")] = (float(dm.group("real_val")), int(dm.group("perc_val")))
                dist = {}
                reg = re.compile("(?P<gene_2>.+)(=)(?P<dist_val>[0-9eE\-\.]+)")
                for d in match.group("dist").split(","):
                    dm = reg.search(d)
                    if dm:
                        dist[dm.group("gene_2")] = float(dm.group("dist_val"))
                result.append({ "rec_id":rec_id, "ref_id":ref_id, "evalue": evalue, "data": data, "dist": dist })

    return result


def get_geotable_json(features):
    "Get JSON for displaying GeoBlast table"
    result = {}
    data_buffer = []
    options = get_config()

    cluster_genes = []
    for feature in features:
        cluster_genes.append(get_gene_id(feature))

    #build dendogram json
    dist_matrix = {}
    for feature in features:
        if "geo" in feature.qualifiers:
            gene_id = get_gene_id(feature)
            geo_array = parse_geo_feature(feature)
            gene_type = "other"
            if "sec_met" in feature.qualifiers:
                gene_type = "biosynthetic"
            for go in geo_array:
                if go["rec_id"] not in result:
                    result[go["rec_id"]] = {
                        "data" : {
                            "feature_names": [],
                            "nodes" : {}
                        },
                        "metadata" : {
                            "feature_names": ["Type"],
                            "nodes" : {}
                        }
                    }
                if gene_id not in result[go["rec_id"]]["data"]["nodes"]:
                    result[go["rec_id"]]["data"]["nodes"][gene_id] = {
                        "count" : 1,
                        "distance" : 0,
                        "features" : [],
                        "objects": [gene_id]
                    }
                if gene_id not in result[go["rec_id"]]["metadata"]["nodes"]:
                    result[go["rec_id"]]["metadata"]["nodes"][gene_id] = [
                        gene_type
                    ]
                samples_order = []
                for geo in options.geo_dataset:
                    if geo["info"]["id"] == go["rec_id"]:
                        samples_order = geo["info"]["samples"]
                def cmp_samples(a, b):
                    if (a in samples_order):
                        if (b in samples_order):
                            return samples_order.index(a) - samples_order.index(b)
                        return 1
                    if (b in samples_order):
                        if (a in samples_order):
                            return samples_order.index(a) - samples_order.index(b)
                        return -1
                    return 0
                for sample_id in sorted([s_id for s_id in go["data"]], cmp = cmp_samples):
                    if sample_id not in result[go["rec_id"]]["data"]["feature_names"]:
                        result[go["rec_id"]]["data"]["feature_names"].append(sample_id)
                    data_buffer.append((go["rec_id"], sample_id, gene_id, go["data"][sample_id][0], go["data"][sample_id][1], go["evalue"]))
                #build dist matrix
                if go["rec_id"] not in dist_matrix:
                    dist_matrix[go["rec_id"]] = {}
                if gene_id not in dist_matrix[go["rec_id"]]:
                    dist_matrix[go["rec_id"]][gene_id] = {}
                for gene_2 in go["dist"]:
                    if gene_2 in cluster_genes:
                        dist_matrix[go["rec_id"]][gene_id][gene_2] = go["dist"][gene_2]

    for data in data_buffer:
        if len(result[data[0]]["data"]["nodes"][data[2]]["features"]) == 0:
            for i in xrange(0, len(result[data[0]]["data"]["feature_names"])):
                result[data[0]]["data"]["nodes"][data[2]]["features"].append(None)
        result[data[0]]["data"]["nodes"][data[2]]["features"][result[data[0]]["data"]["feature_names"].index(data[1])] = data[3]

    # build hierarchical cluster
    for rec_id in result:
        if "data" in result[rec_id] and rec_id in dist_matrix:
            data = result[rec_id]["data"]
            dm = dist_matrix[rec_id]
            result[rec_id]["data"] = build_inchlib_dendrogram(data, dm)

    return result


def build_inchlib_dendrogram(data, dist_matrix, linkage = "complete"):
    """Given 'data' section of inchlib input populated only with leaves, construct
    the whole dendrogram tree based on a supplied distance matrix.
    Available linkage methods are :
        'complete' = complete-linkage clustering
        'single' = single-linkage clustering
        'mean' = average-linkage clustering
    """

    def get_sp(dm):
        sd = 100
        sp = None
        for node_1 in dm:
            for node_2 in dm[node_1]:
                if node_1 != node_2 and (dm[node_1][node_2] < sd or sp == None):
                    sp = (node_1, node_2)
                    sd = dm[node_1][node_2]
        return sp

    def merge(data, dm, sp, nc, linkage):
        nid = "node_%s" % (nc)
        # add the new node to dendogram
        node_1 = data["nodes"][sp[0]]
        node_2 = data["nodes"][sp[1]]
        new_node = {
            "count" : node_1["count"] + node_2["count"],
            "distance" : dm[sp[0]][sp[1]],
            "left_child" : sp[0],
            "right_child" : sp[1]
        }
        node_1["parent"] = nid
        node_2["parent"] = nid
        data["nodes"][nid] = new_node
        # update distance
        new_node = { nid : 0 }
        for node in dm:
            if node not in [sp[0], sp[1]]:
                if linkage == "single":
                    new_node[node] = min(dm[sp[0]][node], dm[sp[1]][node])
                elif linkage == "mean":
                    new_node[node] = int((dm[sp[0]][node] + dm[sp[1]][node]) / 2)
                else: #complete
                    new_node[node] = max(dm[sp[0]][node], dm[sp[1]][node])
        dm[nid] = new_node
        for node in dm:
            if node not in [sp[0], sp[1], nid]:
                dm[node][nid] = dm[nid][node]
        for node in dm:
            dm[node].pop(sp[0], None)
            dm[node].pop(sp[1], None)
        dm.pop(sp[0], None)
        dm.pop(sp[1], None)
        return data, dm

    #Construct dendrogram
    nc = 1
    while len(dist_matrix) > 1:
        shortest_pair = get_sp(dist_matrix)
        data, dm = merge(data, dist_matrix, shortest_pair, nc, linkage)
        nc += 1

    return data


def get_genefinding_basedir(options):
    "Identify the basedir for glimmer files"
    basedir = ''
    if 'glimmer' in options:
        if 'basedir' in options.glimmer:
            basedir = options.glimmer.basedir
    return basedir

def zip_dir(dir_path, archive, prefix_to_remove=""):
    """Recursively add a directory's contents to a zip archive"""
    entries = os.listdir(dir_path)
    for entry in entries:
        entry = path.join(dir_path, entry)
        if path.isdir(entry):
            zip_dir(entry, archive, prefix_to_remove)
        else:
            arcname = entry.replace(prefix_to_remove + os.sep, "", 1)
            archive.write(entry, arcname)

def zip_path(dir_path, name):
    """Create a zip file for a given path"""
    with TemporaryDirectory(change=True):
        try:
            archive = ZipFile(name, 'w', ZIP_DEFLATED)
            if path.isdir(dir_path):
                zip_dir(dir_path, archive, path.dirname(dir_path))
            else:
                archive.write(dir_path)
            archive.close()
        except LargeZipFile:
            archive = ZipFile(name, 'w', ZIP_DEFLATED, True)
            if path.isdir(dir_path):
                zip_dir(dir_path, archive, path.dirname(dir_path))
            else:
                archive.write(dir_path)
            archive.close()

        # small hack to make development testing easier
        if path.exists(path.join(dir_path, name)):
            os.remove(path.join(dir_path, name))
        shutil.move(name, dir_path)

def log_status(status, level='running'):
    """Write status to the status file"""
    options = get_config()
    if 'statusfile' not in options:
        return
    with open(options.statusfile, 'w') as statusfile:
        statusfile.write("%s: %s\n" % (level, status))


def get_git_version():
    """Get the sha1 of the current git version"""
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        out, _, retcode = execute(args)
        return out.strip()
    except OSError:
        pass

    return ""


def get_version():
    """Get the current version string"""
    import antismash
    version = antismash.__version__
    git_version = get_git_version()
    if git_version != '':
        version += "-%s" % git_version

    return version


def get_cluster_cdhit_table(cluster, seq_record, options):
    "Get cd-hit clusters of a Cluster record"
    withinclusterfeatures = []
    for feature in get_cds_features(seq_record):
        if not (cluster.location.start <= feature.location.start <= cluster.location.end or cluster.location.start <= feature.location.end <= cluster.location.end):
            continue
        if not 'sec_met' in feature.qualifiers:
            continue
        if feature not in withinclusterfeatures:
            withinclusterfeatures.append(feature)
    cdhit_table, gene_to_cluster = get_cdhit_table(withinclusterfeatures, options)
    return cdhit_table


def parse_cdhit_file(file_path):
    "Parsing CD-HIT result file into array"
    clusters = []
    gene_to_cluster = {}

    cdhit_file = open(file_path, "r")
    re_content = re.compile(".+\>(?P<gene_id>.+)(\.\.\.) (\*|at (?P<identity_val>[.*0-9]+)).*")
    parsing_content = False
    for line in cdhit_file.readlines():
        if line[0] == '>':
            new_cluster = True
        else:
            match = re_content.match(line)
            if match and (new_cluster or len(clusters) > 0):
                gene_id = match.group("gene_id")
                identity_val = -1.00
                if new_cluster:
                    clusters.append({"rep_gene" : "", "genes" : {}})
                    new_cluster = False
                if match.group("identity_val") == None:
                    clusters[len(clusters) - 1]["rep_gene"] = gene_id
                    identity_val = 100.00
                else:
                    identity_val = float(match.group("identity_val"))
                clusters[len(clusters) - 1]["genes"][gene_id] = identity_val
                gene_to_cluster[gene_id] = len(clusters) - 1
    cdhit_file.close()

    return clusters, gene_to_cluster


def get_cdhit_table(cds_array, options, cutoff = None):
    "Do cdhit analysis on a cds arrays"
    config = get_config()
    if cutoff is None:
        cutoff = float(config.cdh_cutoff)

    #for 0.2 <= cutoff < 0.4, or cd-hit binary not exist, use psi-cd-hit algorithm
    if 0.2 <= cutoff < 0.4 or (0.2 <= cutoff <= 1.0 and not config.cdh_use_binary):
        cds_sorted = sorted(cds_array, key = lambda feature: len(get_aa_sequence(feature)))
        clusters = []
        gene_to_cluster = {}
        while len(cds_sorted) > 0:
            rep_cds = cds_sorted.pop()
            added_cds = []
            clusters.append({"rep_gene" : get_gene_id(rep_cds), "genes" : {}})
            clusters[len(clusters) - 1]["genes"][get_gene_id(rep_cds)] = 100.0
            gene_to_cluster[get_gene_id(rep_cds)] = len(clusters) - 1
            for check_cds in cds_sorted:
                identity_val = get_percentage_identity(rep_cds, check_cds)
                if identity_val >= (cutoff * 100):
                    clusters[len(clusters) - 1]["genes"][get_gene_id(check_cds)] = identity_val
                    gene_to_cluster[get_gene_id(check_cds)] = len(clusters) - 1
                    added_cds.append(check_cds)
            for del_cds in added_cds:
                del cds_sorted[cds_sorted.index(del_cds)]
        return clusters, gene_to_cluster

    #for 0.4 <= cutoff <= 1.0, use cd-hit binary
    elif 0.4 <= cutoff <= 1.0:
        with TemporaryDirectory() as temp_dir:
            file_prefix = path.join(temp_dir, "temp")

            #temporarily assign numbers as gene_id to prevent error from long gene_id
            multifasta_text = ""
            temp_feature_id = []
            for feature in cds_array:
                i = str(len(temp_feature_id))
                temp_feature_id.append(get_gene_id(feature))
                multifasta_text += ">%s\n%s\n" % (i, get_aa_sequence(feature))

            #get temp fasta file to be fed into the cdhit shell command
            fasta_file = open(file_prefix + ".fasta", "w")
            fasta_file.write(multifasta_text)
            fasta_file.close()

            #prep some parameters
            word_length = 1
            if 0.4 <= cutoff <= 0.5:
                word_length = 2
            elif 0.5 <= cutoff <= 0.6:
                word_length = 3
            elif 0.6 <= cutoff <= 0.7:
                word_length = 4
            elif 0.7 <= cutoff <= 1.0:
                word_length = 5

            #execute cd-hit and load result file to an array
            command = ["cd-hit", "-i", fasta_file.name, "-o", file_prefix, "-M", options.cdh_memory,
                        "-c", str(cutoff), "-n", str(word_length), "-T", str(config.cpus), "-B", str(1)]
            error = False
            retcode = 0

            try:
                out, err, retcode = execute(command)
            except OSError:
                error = True

            if retcode != 0:
                logging.debug('cd-hit returned %d: %r while scanning' , retcode, err)
                error = True

            if not error:
                cdhit_table, gene_to_cluster = parse_cdhit_file(file_prefix + ".clstr")
                # restore gene_ids
                for i in [i for i in gene_to_cluster]:
                    gene_to_cluster[temp_feature_id[int(i)]] = gene_to_cluster[i]
                    del gene_to_cluster[i]
                for cluster in cdhit_table:
                    cluster["rep_gene"] = temp_feature_id[int(cluster["rep_gene"])]
                    for i in [i for i in cluster["genes"]]:
                        cluster["genes"][temp_feature_id[int(i)]] = cluster["genes"][i]
                        del cluster["genes"][i]
                return cdhit_table, gene_to_cluster
            else:
                sys.exit(1)
    else:
        logging.debug('cd-hit cutoff is on wrong value (use 0.2 - 1.0)')
        sys.exit(1)


def get_percentage_identity(feature_1, feature_2):
    """Get percentage identity value from 2 features using blastp alignment tools"""
    result = 0.00
    config = get_config()

    with TemporaryDirectory() as temp_dir:
        fasta_file1 = open(path.join(temp_dir, "_tempfasta_1"), "w")
        fasta_file1.write(">%s\n%s" % (get_gene_id(feature_1), feature_1.qualifiers['translation'][0]))
        fasta_file1.close()
        fasta_file2 = open(path.join(temp_dir, "_tempfasta_2"), "w")
        fasta_file2.write(">%s\n%s" % (get_gene_id(feature_2), feature_2.qualifiers['translation'][0]))
        fasta_file2.close()
        command = ["blastp", "-query", fasta_file1.name, "-subject", fasta_file2.name, "-outfmt", "6 pident", "-evalue", "0.0001", "-num_threads",  str(config.cpus)]
        try:
            out, err, retcode = execute(command)
            if (len(out.split("\n")) > 1):
                result = round(float(out.split("\n")[0]), 2)
        except OSError:
            result = 0.00

    return result


def get_pct_identity_table(features):
    """Generate percentage identity tables from array of features"""
    pct_values = {}
    for cds in features:
        cur_gene = get_gene_id(cds)
        if not cur_gene in pct_values.keys():
            pct_values[cur_gene] = {}
        for othercds in [feature for feature in features if feature != cds]:
            other_gene = get_gene_id(othercds)
            if not other_gene in pct_values[cur_gene].keys():
                if other_gene in pct_values.keys():
                    if cur_gene in pct_values[other_gene].keys():
                        pct_values[cur_gene][other_gene] = pct_values[other_gene][cur_gene]
                        continue
                pct_values[cur_gene][other_gene] = get_percentage_identity(cds, othercds)
    return pct_values


if USE_BIOSQL:
    def check_if_dbrecord_exists(name, options):
        """Open database connection; check if record acc exists; close db connection"""

        if "BioSQLconfig" not in options:
            logging.warning("Parameters for database access not defined in default.cfg. Skipping database operations")

        # Set up database object
        myDB = biosql.aSDB(options)

        # connect to namespace for full genomes (dbgenomesnamespace)
        try:
            myDB.connect(namespace=options.BioSQLconfig.dbgenomenamespace)
        except Exception, e:
            logging.exception("Could not connect to database %s, namespace %s : %s",
                              options.BioSQLconfig.dbdb, options.BioSQLconfig.dbgenomenamespace, e)
            return False

        entryid = myDB.fetch_entryid_by_name(name=name)
        myDB.close()

        if entryid:
            logging.debug('check_if_dbrecord_exists: found record id %s for query %s', entryid, name)
            return True
        else:
            return False
