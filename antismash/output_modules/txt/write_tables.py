# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2014 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011-2014 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# Copyright (C) 2014 Tilmann Weber
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""TXT output format module

"""

from antismash import utils
from antismash.output_modules.html.generator import get_detection_rules
from antismash.specific_modules.nrpspks.html_output import _get_monomer_prediction
from antismash.specific_modules.lantipeptides.html_output import _find_core_peptides
import logging


def write_genome(txt, info, options):
    "Write genome table to TXT"
    #TXT columns: genome accession, BGC IDs
    geneclusters = ["%s_c%s" % (info.seq_record.id.partition(".")[0], nr) for nr in info.clusternrs]
    txt.write("\t".join(["genome accession", "geneclusters"]) + "\n")
    txt.write(info.seq_record.id.partition(".")[0] + "\t" + ";".join(geneclusters) + "\n")

def write_transltable(txt, info, options):
    "Write translation table of original sequence IDs (e.g., from FASTA file) to <16 char sequence IDs"
    txt.write("\t".join(["old genome accession", "new genome accession"]) + "\n")
    new_accession = info.seq_record.id.partition(".")[0]
    try:
        old_accession = options.extrarecord[info.seq_record.id].extradata["orig_id"]
    except KeyError:
        old_accession = new_accession
    txt.write("{old_acc}\t{new_acc}\n".format(old_acc=old_accession, new_acc=new_accession))

def write_BGC(txt, info, options):
    "Write BGC table to TXT"
    #TXT columns: BGC ID, BGC_type, detection_rules_used, BGC_range, genes, subclusters,
    # NRPSs_PKSs, signature_genes, RiPPs, pred_structure, monomers
    txt.write("\t".join(["BGC ID", "BGC type", "detection rules used", "BGC_range", "genes", "subclusters", "NRPSs/PKSs", "signature_genes", "RiPPs", "predicted structure", "monomers"]) + "\n")
    for BGCnr in info.clusternrs:
        #Retrieve all data that will be written out
        BGC_ID = "%s_c%s" % (info.seq_record.id.partition(".")[0], BGCnr)
        cluster_feature = utils.get_cluster_by_nr(info.seq_record, BGCnr)
        cluster_gene_features = utils.get_cluster_cds_features(cluster_feature, info.seq_record)
        BGC_type = info.clustertypes[BGCnr].replace("-",";")
        detection_rules_used = '"' + ";".join(get_detection_rules(cluster_feature)) + '"'
        BGC_range = ";".join([str(cluster_feature.location.start), str(cluster_feature.location.end)])
        genes = ";".join(info.accessions[BGCnr])
        if 'subclusterblast' in cluster_feature.qualifiers:
            subclusters = ";".join([qual.partition("\t")[2] for qual in cluster_feature.qualifiers['subclusterblast']])
        else:
            subclusters = ""
        #TODO The subclusterblast module should probably be changed for the precalcs to provide a list here of the 100% hits instead of all hits
        NRPSs_PKSs = ";".join([utils.get_gene_acc(cds).partition(".")[0] for cds in cluster_gene_features if 'sec_met' in cds.qualifiers
            and len([qual for qual in cds.qualifiers['sec_met'] if qual.startswith('NRPS/PKS Domain:')]) > 0])
        signature_genes = ";".join([utils.get_gene_acc(cds).partition(".")[0] for cds in cluster_gene_features if 'sec_met' in cds.qualifiers])
        if len(_find_core_peptides(cluster_feature, info.seq_record)) != 0:
            ripp_list = []
            for peptide in _find_core_peptides(cluster_feature, info.seq_record):
                for cds in cluster_gene_features:
                    if utils.features_overlap(cds, peptide):
                        ripp_list.append(utils.get_gene_acc(cds).partition(".")[0])
                        break
#            RiPPs = ";".join([[utils.get_gene_acc(cds).partition(".")[0] for cds in cluster_gene_features
#                if utils.features_overlap(cds, peptide)][0] for peptide in
 #               _find_core_peptides(cluster_feature, info.seq_record)])
            RiPPs = ";".join(ripp_list)
        else:
            RiPPs = "-"
        if 'structure' in cluster_feature.qualifiers:
            pred_structure = ";".join(cluster_feature.qualifiers['structure'])
        else:
            pred_structure = "N/A"
        monomers = utils.get_structure_pred(cluster_feature)
        #Write data to TXT
        txt.write("\t".join([BGC_ID, BGC_type, detection_rules_used, BGC_range, genes,
            subclusters, NRPSs_PKSs, signature_genes, RiPPs, pred_structure, monomers]) + "\n")

#def write_subclusters(txt, info, options):

def write_signature_gene_info(txt, info, options):
    "Write signature gene table to TXT"
    #TXT columns: signature_gene, pHMM_hit, e-value, bit score, nr of seeds
    txt.write("\t".join(["signature gene", "pHMM hits", "e-value", "bit score", "number of seeds"]) + "\n")
    for BGCnr in info.clusternrs:
        #Retrieve all data that will be written out
        cluster_feature = utils.get_cluster_by_nr(info.seq_record, BGCnr)
        cluster_gene_features = utils.get_cluster_cds_features(cluster_feature, info.seq_record)
        signature_genes = [cds for cds in cluster_gene_features if 'sec_met' in cds.qualifiers]
        for cds in signature_genes:
            if len([qual for qual in cds.qualifiers['sec_met'] if qual.startswith('Domains detected: ')]) == 0:
                continue
            gene_ID = utils.get_gene_acc(cds).partition(".")[0]
            domdetect_qual = [qual for qual in cds.qualifiers['sec_met'] if qual.startswith('Domains detected: ')][0]
            if ";" in domdetect_qual:
                domains = domdetect_qual.partition("Domains detected: ")[2].split(";")
            else:
                domains = [domdetect_qual.partition("Domains detected: ")[2]]
            for domain in domains:
                domain_name = domain.partition(" (")[0].replace(" ", "")
                evalue = domain.partition("E-value: ")[2].partition(",")[0]
                bitscore = domain.partition("bitscore: ")[2].partition(",")[0]
                nr_seeds = domain.partition("seeds: ")[2].partition(")")[0]
                txt.write("\t".join([gene_ID, domain_name, evalue, bitscore, nr_seeds]) + "\n")

def write_gene(txt, info, options):
    "Write gene table to TXT"
    #TXT columns: gene ID, gene start, gene end, gene strand, smCOG, locus_tag/geneID, annotation
    txt.write("\t".join(["gene ID", "gene start", "gene end", "gene strand", "smCOG", "locus_tag", "annotation"]) + "\n")
    for BGCnr in info.clusternrs:
        #Retrieve all data that will be written out
        cluster_feature = utils.get_cluster_by_nr(info.seq_record, BGCnr)
        cluster_gene_features = utils.get_cluster_cds_features(cluster_feature, info.seq_record)
        for cds in cluster_gene_features:
            gene_id = utils.get_gene_acc(cds).partition(".")[0]
            cds_start = str(cds.location.start)
            cds_end = str(cds.location.end)
            if cds.strand == 1:
                cds_strand = "+"
            else:
                cds_strand = "-"
            smCOG = "" ##Not used for now
            locus_tag = utils.get_gene_id(cds).partition(".")[0]
            annotation = utils.get_gene_annotation(cds)
            txt.write("\t".join([gene_id, cds_start, cds_end, cds_strand, smCOG, locus_tag, annotation]) + "\n")

def write_NRPS_PKS(txt, info, options):
    "Write NRPS/PKS table to TXT"
    #TXT columns: NRPS/PKS ID, annotation, aSDomain, score, evalue, domain type, subtype, range, activity, NRPSPredictor2, Stachelhaus, Minowa, pkssignature, consensus
    txt.write("\t".join(["Cluster_ID", "NRPSPKS_ID", "annotation", "aSDomain", "score", "evalue", "domain_type", "subtype", "domain_start", "domain_end", "KR activity", "KR stereochemistry", "NRPSPredictor2", "Stachelhaus", "Minowa", "pkssignature", "consensus"]) + "\n")
    for BGCnr in info.clusternrs:
        #Retrieve all data that will be written out
        cluster_feature = utils.get_cluster_by_nr(info.seq_record, BGCnr)
        cluster_gene_features = utils.get_cluster_cds_features(cluster_feature, info.seq_record)
        cluster_id = "{seq_id}_c{cluster_nr}".format(seq_id=info.seq_record.id, cluster_nr=BGCnr)
        NRPSs_PKSs = [cds for cds in cluster_gene_features if 'sec_met' in cds.qualifiers
            and len([qual for qual in cds.qualifiers['sec_met'] if qual.startswith('NRPS/PKS Domain:')]) > 0]
        for cds in NRPSs_PKSs:
            enzyme_ID = utils.get_gene_acc(cds).partition(".")[0]
            if len([qual for qual in cds.qualifiers['sec_met'] if "NRPS/PKS subtype: " in qual]) > 0:
                enzyme_annotation = [qual for qual in cds.qualifiers['sec_met'] if qual.startswith("NRPS/PKS subtype")][0].partition("NRPS/PKS subtype: ")[2]
            else:
                logging.warn("No enzyme annotation for %s" % enzyme_ID)
                enzyme_annotation = ""
            aSDomains = [dom for dom in utils.get_cluster_aSDomain_features(cluster_feature, info.seq_record) if utils.features_overlap(cds, dom) and utils.get_gene_id(cds) in [dom.qualifiers['locus_tag'], dom.qualifiers['locus_tag'][0]]]
            for aSDomain in aSDomains:
                domtype = aSDomain.qualifiers['domain'][0]
                if "domain_subtype" in aSDomain.qualifiers:
                    subtype = aSDomain.qualifiers['domain_subtype'][0]
                else:
                    subtype = ""
                aSDomain_ID = aSDomain.qualifiers['asDomain_id'][0]
                score = str(aSDomain.qualifiers['score'][0])
                evalue = str(aSDomain.qualifiers['evalue'][0])
                dom_start = str(aSDomain.location.start)
                dom_end = str(aSDomain.location.end)
                kr_activity = ""
                kr_stereochemistry = ""
                NRPSPredictor2 = ""
                Stachelhaus = ""
                Minowa = ""
                pkssignature = ""
                consensus = ""
                if aSDomain.qualifiers.has_key('specificity'):
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("KR activity")]) > 0:
                        kr_activity = [qual.partition("KR activity: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("KR activity")][0]
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("KR stereochemistry")]) > 0:
                        kr_stereochemistry = [qual.partition("KR stereochemistry: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("KR stereochemistry")][0]
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("NRPSpredictor2")]) > 0:
                        NRPSPredictor2 = [qual.partition("NRPSpredictor2 SVM: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("NRPSpredictor2")][0]
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("Stachelhaus")]) > 0:
                        Stachelhaus = [qual.partition("Stachelhaus code: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("Stachelhaus")][0]
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("Minowa")]) > 0:
                        Minowa = [qual.partition("Minowa: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("Minowa")][0]
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("PKS signature")]) > 0:
                        pkssignature = [qual.partition("PKS signature: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("PKS signature")][0]
                    if len([qual for qual in aSDomain.qualifiers['specificity'] if qual.startswith("consensus")]) > 0:
                        consensus = [qual.partition("consensus: ")[2] for qual in aSDomain.qualifiers['specificity'] if qual.startswith("consensus")][0]

                txt.write("\t".join([cluster_id, enzyme_ID, enzyme_annotation, aSDomain_ID, score, evalue, domtype, subtype, dom_start, dom_end, kr_activity, kr_stereochemistry, NRPSPredictor2, Stachelhaus, Minowa, pkssignature, consensus]) + "\n")

def write_smCOG(txt, info, options):
    pass
#    WILL BE CALCULATED SEPARATELY

def write_RiPP(txt, info, options):
    "Write RiPP table to TXT"
    #TXT columns: RiPP ID, annotation, core peptide, mol weight, monoisotopic_mass, alt mol weights, nr bridges
    txt.write("\t".join(["RiPP ID", "annotation", "core peptide", "molecular weight", "monoisotopic_mass", "alternative molecular weights", "number of bridges"]) + "\n")
    for BGCnr in info.clusternrs:
        #Retrieve all data that will be written out
        cluster_feature = utils.get_cluster_by_nr(info.seq_record, BGCnr)
        cluster_gene_features = utils.get_cluster_cds_features(cluster_feature, info.seq_record)
        RiPP_features = _find_core_peptides(cluster_feature, info.seq_record)
        RiPPs = []
        for peptide in RiPP_features:
            for cds in cluster_gene_features:
                if utils.features_overlap(cds, peptide):
                    RiPPs.append(utils.get_gene_acc(cds).partition(".")[0])
                    break
        idx = 0
        for RiPP in RiPP_features:
            RiPP_ID = RiPPs[idx]
            note_quals = RiPP.qualifiers['note']
            annotation = [qual.partition("predicted class: ")[2] for qual in note_quals if "predicted class:" in qual][0]
            core_peptide = [qual.partition("predicted core seq: ")[2] for qual in note_quals if "predicted core seq:" in qual][0]
            mol_weight = [qual.partition("molecular weight: ")[2] for qual in note_quals if "molecular weight: " in qual][0]
            monoiso_mass = [qual.partition("monoisotopic mass: ")[2] for qual in note_quals if "monoisotopic mass: " in qual][0]
            if "alternative weights" in note_quals:
                alt_mol_weights = [qual.partition("alternative weights: ")[2].replace(" ", "") for qual in note_quals if "alternative weights:" in qual][0]
            else:
                alt_mol_weights = ""
            nr_bridges = [qual.partition("number of bridges: ")[2] for qual in note_quals if "number of bridges: " in qual][0]
            txt.write("\t".join([RiPP_ID, annotation, core_peptide, mol_weight, monoiso_mass, alt_mol_weights, nr_bridges]) + "\n")
            idx += 1
