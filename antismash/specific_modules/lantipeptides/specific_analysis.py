# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2012 Daniyal Kazempour
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# Copyright (C) 2012,2013 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
'''
More detailed lantipeptide analysis using HMMer-based leader peptide
cleavage site prediction as well as prediction of number of lanthionine
bridges and molcular mass.
'''

import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from antismash import utils

known_precursor_domains = (
    'Antimicr18',
    'Gallidermin',
    'L_biotic_A',
    'lacticin_l',
    'leader_d',
    'leader_abc',
    'leader_eh',
    'mature_ha',
    'lacticin_mat',
    'mature_b',
    'mature_d',
    'mature_ab',
    'mature_a',
    'TIGR03731',
    'LD_lanti_pre',
    'strep_PEQAXS',
)

class Lantipeptide(object):
    '''
    Class to calculate and store lantipeptide information
    '''
    def __init__(self, start, end, score, lantype):
        self.start = start
        self.end = end
        self.score = score
        self.lantype = lantype
        self._leader = ''
        self._core = ''
        self._weight = -1
        self._monoisotopic_weight = -1
        self._alt_weights = []
        self._lan_bridges = -1
        self._aminovinyl = False
        self._chlorinated = False
        self._oxygenated = False
        self._lac = False

    @property
    def core(self):
        return self._core

    @core.setter
    def core(self, seq):
        seq = seq.replace('X', '')
        self.core_analysis_monoisotopic = ProteinAnalysis(seq, monoisotopic=True)
        self.core_analysis = ProteinAnalysis(seq, monoisotopic=False)
        self._core = seq

    @property
    def leader(self):
        return self._leader

    @leader.setter
    def leader(self, seq):
        self._leader = seq

    def __repr__(self):
        return "Lantipeptide(%s..%s, %s, %r, %r, %s, %s(%s))" % (self.start, self.end, self.score, self.lantype, self._core, self._lan_bridges, self._monoisotopic_weight, self._weight)

    @property
    def number_of_lan_bridges(self):
        '''
        function determines the number of lanthionine bridges in the core peptide
        '''
        if self._lan_bridges > -1:
            return self._lan_bridges

        aas = self.core_analysis.count_amino_acids()
        no_cys = aas['C']
        no_thr_ser = aas['T'] + aas['S']
        self._lan_bridges = min(no_cys, no_thr_ser)
        if self._aminovinyl:
            self._lan_bridges -= 1
        return self._lan_bridges

    def _calculate_mw(self):
        '''
        (re)calculate the monoisotopic mass and molecular weight
        '''
        if not self._core:
            raise ValueError()

        aas = self.core_analysis.count_amino_acids()
        no_thr_ser = aas['T'] + aas['S']

        mol_mass = self.core_analysis.molecular_weight()
        mods = 18.02 * no_thr_ser
        if self._aminovinyl:
            mods += 46
        if self._chlorinated:
            mods -= 34
        if self._oxygenated:
            mods -= 16
        if self._lac:
            mods -= 2
        self._weight = mol_mass - mods

        # every unbridged Ser or Thr might not be dehydrated
        self._alt_weights = []
        for i in range(1, no_thr_ser - aas['C'] + 1):
            self._alt_weights.append(self._weight + 18.02 * i)

        monoisotopic_mass = self.core_analysis_monoisotopic.molecular_weight()
        mods = 18 * no_thr_ser
        if self._aminovinyl:
            mods += 46
        if self._chlorinated:
            mods -= 34
        if self._oxygenated:
            mods -= 16
        if self._lac:
            mods -= 2
        self._monoisotopic_weight = monoisotopic_mass - mods

    @property
    def monoisotopic_mass(self):
        '''
        function determines the weight of the core peptide and substracts
        the weight which is reduced, due to dehydratation
        '''
        if self._monoisotopic_weight > -1:
            return self._monoisotopic_weight

        self._calculate_mw()
        return self._monoisotopic_weight

    @property
    def molecular_weight(self):
        '''
        function determines the weight of the core peptide and substracts
        the weight which is reduced, due to dehydratation
        '''
        if self._weight > -1:
            return self._weight

        self._calculate_mw()
        return self._weight

    @property
    def alternative_weights(self):
        '''
        function determines the possible alternative weights assuming one or
        more of the Ser/Thr residues aren't dehydrated
        '''
        if self._alt_weights != []:
            return self._alt_weights

        self._calculate_mw()
        return self._alt_weights

    @property
    def aminovinyl_group(self):
        '''
        Check if lantipeptide contains an aminovinyl group
        '''
        return self._aminovinyl

    @aminovinyl_group.setter
    def aminovinyl_group(self, value):
        '''
        Define if lantipeptide contains an aminovinyl group and trigger
        recalculation of the molecular weight if needed
        '''
        self._aminovinyl = value
        if self._core:
            self._calculate_mw()
            # recalculate the number of lan bridges
            self._lan_bridges = -1

    @property
    def chlorinated(self):
        '''
        Check if lantipeptide is chlorinated
        '''
        return self._chlorinated

    @chlorinated.setter
    def chlorinated(self, value):
        '''
        Define if lantipeptide is chlorinated and trigger
        recalculation of the molecular weight if needed
        '''
        self._chlorinated = value
        if self._core:
            self._calculate_mw()

    @property
    def oxygenated(self):
        '''
        Check if lantipeptide is oxygenated
        '''
        return self._oxygenated

    @oxygenated.setter
    def oxygenated(self, value):
        '''
        Define if lantipeptide is oxygenated and trigger
        recalculation of the molecular weight if needed
        '''
        self._oxygenated = value
        if self._core:
            self._calculate_mw()

    @property
    def lactonated(self):
        '''
        Check if lantipeptide starts with a lactone
        '''
        return self._lac

    @lactonated.setter
    def lactonated(self, value):
        self._lac = value
        if self._core:
            self._calculate_mw()

def predict_cleavage_site(query_hmmfile, target_sequence):
    '''
    Function extracts from HMMER the start position, end position, score and lant-type
    of the HMM alignment
    '''
    hmmer_res = utils.run_hmmpfam2(query_hmmfile, target_sequence)
    resvec = None
    for res in hmmer_res:
        for hits in res:
            lanti_type = hits.description
            for hsp in hits:
                resvec = Lantipeptide(hsp.query_start-1, hsp.query_end, hsp.bitscore, lanti_type)
                return resvec
    return resvec


def predict_class_from_gene_cluster(seq_record, cluster):
    '''
    Predict the lantipeptide class from the gene cluster
    '''
    found_domains = []
    for feature in utils.get_cds_features(seq_record):
        if feature.location.start < cluster.location.start or \
           feature.location.end > cluster.location.end:
            continue

        if not 'sec_met' in feature.qualifiers:
            continue

        for entry in feature.qualifiers['sec_met']:
            if entry.startswith('Domains detected:'):
                entry = entry[17:]
                domains = entry.split(';')
                for domain in domains:
                    found_domains.append(domain.split()[0])

    if 'Lant_dehyd_N' in found_domains or 'Lant_dehyd_C' in found_domains:
        return 'Class-I'
    if 'DUF4135' in found_domains:
        return 'Class-II'
    if 'Pkinase' in found_domains:
        # this could be class 3 or class 4, but as nobody has seen class 4
        # in vivo yet, we'll ignore that
        return 'Class-III'

    # Ok, no biosynthetic enzymes found, let's try the prepeptide
    if 'Gallidermin' in found_domains:
        return 'Class-I'

    return None


thresh_dict = {
        'Class-I' : -15,
        'Class-II' : -7,
        'Class-III' : 5,
    }


def run_lantipred(seq_record, query, lant_class):
    hmmer_profiles = {'Class-I': 'class1.hmm',
                      'Class-II':'class2.hmm',
                      'Class-III': 'class3.hmm', }

    query_sequence = utils.get_aa_sequence(query, to_stop=True)
    lan_a_fasta = ">%s\n%s" % (utils.get_gene_id(query), query_sequence)

    #run sequence against profiles and parse them in a vector containing START, END, SCORE and LANTYPE
    profile = utils.get_full_path(__file__, hmmer_profiles[lant_class])
    result = predict_cleavage_site(profile, lan_a_fasta)

    if result is None:
        logging.debug('%r: No cleavage site predicted' % utils.get_gene_id(query))
        return

    if thresh_dict[lant_class] > result.score:
        logging.debug('%r: Score %0.2f below threshold %0.2f for class %r' %
                      (utils.get_gene_id(query), result.score,
                       thresh_dict[lant_class], lant_class))
        return

    #extract now (that class is known and thus the END component) the core peptide
    result.leader = query_sequence[:result.end]
    result.core = query_sequence[result.end:]
    if result.core.find('C') < 0:
        logging.debug('%r: No Cysteine residues found in core, false positive' %
                      utils.get_gene_id(query))
        return
    if not 'sec_met' in query.qualifiers:
        query.qualifiers['sec_met'] = []


    if ";".join(query.qualifiers['sec_met']).find(';Kind: biosynthetic') < 0:
        query.qualifiers['sec_met'].append('Kind: biosynthetic')

    return result


def find_lan_a_features(seq_record, cluster):
    lan_a_features = []
    for feature in utils.get_cds_features(seq_record):
        if feature.location.start < cluster.location.start or \
           feature.location.end > cluster.location.end:
            continue

        aa_seq = utils.get_aa_sequence(feature)
        if len(aa_seq) < 80:
            lan_a_features.append(feature)
            continue

        if not 'sec_met' in feature.qualifiers:
            continue

        domain = None
        for entry in feature.qualifiers['sec_met']:
            if entry.startswith('Domains detected:'):
                domain = entry.split()[2]
                break

        if domain is None:
            continue

        if domain not in known_precursor_domains:
            continue

        lan_a_features.append(feature)

    return lan_a_features


def find_flavoprotein(seq_record, cluster):
    "Look for an epiD-like flavoprotein responsible for aminovinylcystein"
    for feature in utils.get_cds_features(seq_record):
        if feature.location.start < cluster.location.start or \
           feature.location.end > cluster.location.end:
            continue

        if not 'sec_met' in feature.qualifiers:
            continue

        domain = None
        for entry in feature.qualifiers['sec_met']:
            if entry.startswith('Domains detected:'):
                domain = entry.split()[2]
                break

        if domain is None:
            continue

        if domain in 'Flavoprotein':
            return True

    return False


def find_halogenase(seq_record, cluster):
    "Look for a halogenase"
    #return False
    for feature in utils.get_cds_features(seq_record):
        if feature.location.start < cluster.location.start or \
           feature.location.end > cluster.location.end:
            continue

        if not 'sec_met' in feature.qualifiers:
            continue

        domain = None
        for entry in feature.qualifiers['sec_met']:
            if entry.startswith('Domains detected:'):
                domain = entry.split()[2]
                break

        if domain is None:
            continue

        if domain in 'Trp_halogenase':
            return True

    return False


def find_p450_oxygenase(seq_record, cluster):
    "Look for a p450 oxygenase"
    #return False
    for feature in utils.get_cds_features(seq_record):
        if feature.location.start < cluster.location.start or \
           feature.location.end > cluster.location.end:
            continue

        if not 'sec_met' in feature.qualifiers:
            continue

        domain = None
        for entry in feature.qualifiers['sec_met']:
            if entry.startswith('Domains detected:'):
                domain = entry.split()[2]
                break

        if domain is None:
            continue

        if domain in 'p450':
            return True

    return False


def find_short_chain_dehydrogenase(seq_record, cluster):
    "Look for an eciO-like short-chain dehydrogenase responsible for N-terminal lactone"
    for feature in utils.get_cds_features(seq_record):
        if feature.location.start < cluster.location.start or \
           feature.location.end > cluster.location.end:
            continue

        if not 'sec_met' in feature.qualifiers:
            continue

        domain = None
        for entry in feature.qualifiers['sec_met']:
            if entry.startswith('Domains detected:'):
                domain = entry.split()[2]
                break

        if domain is None:
            continue

        if domain in ('adh_short', 'adh_short_C2'):
            return True

    return False


def result_vec_to_features(orig_feature, res_vec):
    start = orig_feature.location.start
    end = orig_feature.location.start + (res_vec.end * 3)
    strand = orig_feature.location.strand
    loc = FeatureLocation(start, end, strand=strand)
    leader_feature = SeqFeature(loc, type='CDS_motif')
    leader_feature.qualifiers['note'] = ['leader peptide']
    leader_feature.qualifiers['note'].append('predicted leader seq: %s' % res_vec.leader)
    leader_feature.qualifiers['locus_tag'] = [ utils.get_gene_id(orig_feature) ]

    start = end
    end = orig_feature.location.end
    loc = FeatureLocation(start, end, strand=strand)
    core_feature = SeqFeature(loc, type='CDS_motif')
    core_feature.qualifiers['note'] = ['core peptide']
    core_feature.qualifiers['note'].append('monoisotopic mass: %0.1f' % res_vec.monoisotopic_mass)
    core_feature.qualifiers['note'].append('molecular weight: %0.1f' % res_vec.molecular_weight)
    if res_vec.alternative_weights:
        weights = map(lambda x: "%0.1f" % x, res_vec.alternative_weights)
        core_feature.qualifiers['note'].append('alternative weights: %s' % "; ".join(weights))
    core_feature.qualifiers['note'].append('number of bridges: %s' % res_vec.number_of_lan_bridges)
    core_feature.qualifiers['note'].append('predicted core seq: %s' % res_vec.core)
    core_feature.qualifiers['note'].append('predicted class: %s' % res_vec.lantype)
    core_feature.qualifiers['note'].append('score: %0.2f' % res_vec.score)
    if res_vec.aminovinyl_group:
        core_feature.qualifiers['note'].append('predicted additional modification: AviCys')
    if res_vec.chlorinated:
        core_feature.qualifiers['note'].append('predicted additional modification: Cl')
    if res_vec.oxygenated:
        core_feature.qualifiers['note'].append('predicted additional modification: OH')
    if res_vec.lactonated:
        core_feature.qualifiers['note'].append('predicted additional modification: Lac')
    core_feature.qualifiers['locus_tag'] = [ utils.get_gene_id(orig_feature) ]

    return [leader_feature, core_feature]


def specific_analysis(seq_record, options):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        if 'product' not in cluster.qualifiers or \
           'lantipeptide' not in cluster.qualifiers['product'][0]:
            continue

        lan_as = find_lan_a_features(seq_record, cluster)
        flavoprotein_found = find_flavoprotein(seq_record, cluster)
        halogenase_found = find_halogenase(seq_record, cluster)
        oxygenase_found = find_p450_oxygenase(seq_record, cluster)
        dehydrogenase_found = find_short_chain_dehydrogenase(seq_record, cluster)
        lant_class = predict_class_from_gene_cluster(seq_record, cluster)
        if lant_class is None:
            return
        for lan_a in lan_as:
            result_vec = run_lantipred(seq_record, lan_a, lant_class)
            if result_vec is None:
                continue
            if flavoprotein_found:
                result_vec.aminovinyl_group = True
            if halogenase_found:
                result_vec.chlorinated = True
            if oxygenase_found:
                result_vec.oxygenated = True
            if dehydrogenase_found and result_vec.core.startswith('S'):
                result_vec.lactonated = True

            new_features = result_vec_to_features(lan_a, result_vec)
            seq_record.features.extend(new_features)

