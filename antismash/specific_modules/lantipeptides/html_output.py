# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011-2013 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from pyquery import PyQuery as pq
from antismash import utils

def will_handle(product):
    if product.find('lantipeptide') > -1:
        return True

    return False

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    """Generate details div"""

    cluster_rec = utils.get_cluster_by_nr(seq_record, cluster['idx'])
    if cluster_rec is None:
        return details

    leader_peptides = _find_leader_peptides(cluster_rec, seq_record)
    core_peptides = _find_core_peptides(cluster_rec, seq_record)

    if details is None:
        details = pq('<div>')
        details.addClass('details')

        header = pq('<h3>')
        header.text('Detailed annotation')
        details.append(header)

    if len(core_peptides) == 0:
        details_text = pq('<div>')
        details_text.addClass('details-text')
        details_text.text('No core peptides found.')
        details.append(details_text)
        return details

    details_text = pq('<dl>')
    details_text.addClass('details-text')

    i = 0
    for cp in core_peptides:
        leader = leader_peptides[i]
        leader_seq = _get_leader_peptide_sequence(leader)
        core_seq = _get_core_peptide_sequence(cp)
        dt = pq('<dt>')
        dt.text('%s leader / core peptide, putative %s' % (utils.get_gene_id(cp),
                _get_core_peptide_class(cp)))
        details_text.append(dt)

        dd = pq('<dd>')
        core_seq = core_seq.replace('S', '<span class="dha">Dha</span>')
        core_seq = core_seq.replace('T', '<span class="dhb">Dhb</span>')
        core_seq = core_seq.replace('C', '<span class="cys">C</span>')
        seq = "%s - %s" % (leader_seq, core_seq)
        dd.html(seq)
        details_text.append(dd)
        i += 1

    details.append(details_text)

    legend = pq('<div>')
    legend.addClass('legend')
    legend_header = pq('<h4>')
    legend_header.text('Legend:')
    legend.append(legend_header)

    legend_text = pq('<div>')
    legend_text.html('<span class="dha">Dha</span>: Didehydroalanine<br>'
                '<span class="dhb">Dhb</span>: Didehydrobutyrine')
    legend.append(legend_text)
    details.append(legend)

    return details

def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    """Generate sidepanel div"""
    cluster_rec = utils.get_cluster_by_nr(seq_record, cluster['idx'])
    if cluster_rec is None:
        return sidepanel

    if sidepanel is None:
        sidepanel = pq('<div>')
        sidepanel.addClass('sidepanel')

    core_peptides = _find_core_peptides(cluster_rec, seq_record)
    if len(core_peptides) == 0:
        return sidepanel

    details = pq('<div>')
    details.addClass('more-details')
    details_header = pq('<h3>')
    details_header.text('Prediction details')
    details.append(details_header)
    details_list = pq('<dl>')
    details_list.addClass('prediction-text')

    for cp in core_peptides:
        dt = pq('<dt>')
        dt.text(utils.get_gene_id(cp))
        details_list.append(dt)
        dd = pq('<dd>')
        mass = _get_monoisotopic_mass(cp)
        mol_weight = _get_molecular_weight(cp)
        bridges = _get_number_bridges(cp)
        pred_class = _get_core_peptide_class(cp)
        score = _get_core_peptide_score(cp)
        dd.html('Putative %s<br>Score: %0.2f<br>Monoisotopic mass: %s Da<br>'\
                'Molecular weight: %s Da<br>Number of bridges: %s' %\
                (pred_class, score, mass, mol_weight, bridges))
        for mod in _get_core_peptide_extra_modifications(cp):
            dd.html('%s<br>Additional modifications: %s' % (dd.html(), mod))
        _alt_weights = _get_alternative_weights(cp)
        if _alt_weights:
            inner_dl = pq('<dl>')
            inner_dt = pq('<dt>')
            inner_dt.text('Alternative weights')
            inner_dl.append(inner_dt)
            inner_dd = pq('<dd>')
            inner_dd.addClass('alt-weight-desc')
            inner_dd.text('(assuming N unmodified Ser/Thr residues)')
            inner_dl.append(inner_dd)
            i = 1
            for weight in _alt_weights:
                inner_dd = pq('<dd>')
                weight_span = pq('<span>')
                weight_span.text('%0.1f Da' % weight)
                weight_span.addClass('alt-weight')
                n_span = pq('<span>')
                n_span.text('N = %d' % i)
                n_span.addClass('alt-weight-n')
                inner_dd.append(weight_span)
                inner_dd.append(n_span)
                inner_dl.append(inner_dd)
                i += 1
            dd.append(inner_dl)
        details_list.append(dd)

    details.append(details_list)
    sidepanel.append(details)

    cross_refs = pq("<div>")
    refs_header = pq('<h3>')
    refs_header.text('Database cross-links')
    cross_refs.append(refs_header)
    links = pq("<div>")
    links.addClass('prediction-text')

    a = pq("<a>")
    a.attr('href', 'http://bioinfo.lifl.fr/norine/form2.jsp')
    a.attr('target', '_new')
    a.text("Look up in NORINE database")
    links.append(a)
    cross_refs.append(links)
    sidepanel.append(cross_refs)

    return sidepanel

def _find_core_peptides(cluster, seq_record):
    """Find CDS_motifs containing lantipeptide core peptide annotations"""
    motifs = []
    for motif in utils.get_all_features_of_type(seq_record, 'CDS_motif'):
        if motif.location.start < cluster.location.start or \
           motif.location.end > cluster.location.end:
            continue

        if not motif.qualifiers.has_key('note'):
            continue

        if not 'core peptide' in motif.qualifiers['note']:
            continue

        motifs.append(motif)

    return motifs

def _find_leader_peptides(cluster, seq_record):
    """Find CDS_motifs containing lantipeptide leader peptide annotations"""
    motifs = []
    for motif in utils.get_all_features_of_type(seq_record, 'CDS_motif'):
        if motif.location.start < cluster.location.start or \
           motif.location.end > cluster.location.end:
            continue

        if not motif.qualifiers.has_key('note'):
            continue

        if not 'leader peptide' in motif.qualifiers['note']:
            continue

        motifs.append(motif)

    return motifs

def _get_monoisotopic_mass(motif):
    """Get monoisotopic mass of a core peptide motif"""
    for note in motif.qualifiers['note']:
        if not note.startswith('monoisotopic mass:'):
            continue
        return float(note.split(':')[-1])

def _get_molecular_weight(motif):
    """Get molecular weight of a core peptide motif"""
    for note in motif.qualifiers['note']:
        if not note.startswith('molecular weight:'):
            continue
        return float(note.split(':')[-1])

def _get_alternative_weights(motif):
    """Get alternative weights assuming less Ser/Thr residues are dehydrated"""
    alternative_weights = []
    for note in motif.qualifiers['note']:
        if not note.startswith('alternative weights:'):
            continue
        weights = note.split(':')[-1].split(';')
        alternative_weights = map(lambda x: float(x), weights)
    return alternative_weights

def _get_number_bridges(motif):
    """Get number of bridges of a core peptide motif"""
    for note in motif.qualifiers['note']:
        if not note.startswith('number of bridges:'):
            continue
        return int(note.split(':')[-1])

def _get_leader_peptide_sequence(motif):
    """Get AA sequence of a core peptide motif"""
    for note in motif.qualifiers['note']:
        if not note.startswith('predicted leader seq:'):
            continue
        return note.split(':')[-1].strip()

def _get_core_peptide_sequence(motif):
    """Get AA sequence of a core peptide motif"""
    for note in motif.qualifiers['note']:
        if not note.startswith('predicted core seq:'):
            continue
        return note.split(':')[-1].strip()

def _get_core_peptide_class(motif):
    """Get the predicted class of the core peptide"""
    for note in motif.qualifiers['note']:
        if not note.startswith('predicted class:'):
            continue
        return note.split(':')[-1].strip().replace('-', ' ')

def _get_core_peptide_score(motif):
    """Get the score of the core peptide prediction"""
    for note in motif.qualifiers['note']:
        if not note.startswith('score:'):
            continue
        return float(note.split(':')[-1].strip())

def _get_core_peptide_extra_modifications(motif):
    """Get additional modifications on the core peptide"""
    modifications = []
    for note in motif.qualifiers['note']:
        if not note.startswith('predicted additional modification:'):
            continue
        modifications.append(note.split(':')[-1].strip())
    return modifications
