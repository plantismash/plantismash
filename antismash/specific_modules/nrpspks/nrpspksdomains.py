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
from collections import defaultdict
from antismash import utils
from antismash.lib.hmmscanparser import parse_hmmscan_results
from Bio.SeqFeature import SeqFeature, FeatureLocation


def annotate_pksnrps(pksnrpsvars, seq_record, options):
    withinclustergenes = utils.get_withincluster_cds_features(seq_record)
    if len(withinclustergenes) == 0:
        logging.debug('No genes within a sec_met cluster found for %r' % seq_record.id)
        return pksnrpsvars
    run_nrpspks_specific_hmmer(seq_record, withinclustergenes, pksnrpsvars)
    name_nrpspks(seq_record, pksnrpsvars, withinclustergenes, options)
    pksnrpsvars.pksnrpscoregenes = utils.get_pksnrps_cds_features(seq_record)
    return pksnrpsvars

def filter_nonterminal_docking_domains(seq_record, pksnrpsvars):
    dockingdomains = ['NRPS-COM_Nterm', 'NRPS-COM_Cterm', 'PKS_Docking_Cterm', 'PKS_Docking_Nterm']
    hitgenes = pksnrpsvars.domaindict.keys()
    feature_by_id = utils.get_feature_dict(seq_record)
    for hitgene in hitgenes:
        to_remove = []
        cdsfeature = feature_by_id[hitgene]
        cds_seq = utils.get_aa_sequence(cdsfeature)
        hitgenelength = len(cds_seq)
        x = 0
        for hit in pksnrpsvars.domaindict[hitgene]:
            if hit[0] in dockingdomains:
                if not (hitgenelength - max(hit[1], hit[2]) < 50 or min(hit[1], hit[2]) < 50):
                    to_remove.append(x)
            x += 1
        to_remove.reverse()
        for idx in to_remove:
            del pksnrpsvars.domaindict[hitgene][idx]
        if pksnrpsvars.domaindict[hitgene] == []:
            del pksnrpsvars.domaindict[hitgene]

def run_nrpspks_specific_hmmer(seq_record, withinclustergenes, pksnrpsvars):
    nrpspksfasta = utils.get_specific_multifasta(withinclustergenes)
    #Analyse for abMotifs
    abmotif_opts = ["-E", "0.25"]
    abmotif_results = utils.run_hmmscan(utils.get_full_path(__file__, "abmotifs.hmm"), nrpspksfasta, abmotif_opts)
    mhmmlengthsdict = utils.hmmlengths(utils.get_full_path(__file__, "abmotifs.hmm"))
    pksnrpsvars.motifdict = parse_hmmscan_results(abmotif_results, mhmmlengthsdict)
    #Analyse for C/A/PCP/E/KS/AT/ATd/DH/KR/ER/ACP/TE/TD/COM/Docking/MT/CAL domains
    nrpspksdomain_opts = ["--cut_tc"]
    nrpspksdomain_results = utils.run_hmmscan(utils.get_full_path(__file__, "nrpspksdomains.hmm"), nrpspksfasta, nrpspksdomain_opts)
    hmmlengthsdict = utils.hmmlengths(utils.get_full_path(__file__, "nrpspksdomains.hmm"))
    pksnrpsvars.domaindict = parse_hmmscan_results(nrpspksdomain_results, hmmlengthsdict)
    filter_nonterminal_docking_domains(seq_record, pksnrpsvars)
    #Analyse KS domains & PKS/NRPS protein domain composition to detect NRPS/PKS types
    kshmmlengthsdict = utils.hmmlengths(utils.get_full_path(__file__, "ksdomains.hmm"))
    ksdomain_results = utils.run_hmmscan(utils.get_full_path(__file__, "ksdomains.hmm"), nrpspksfasta, nrpspksdomain_opts)
    pksnrpsvars.ksdomaindict = parse_hmmscan_results(ksdomain_results, kshmmlengthsdict)


def name_nrpspks(seq_record, pksnrpsvars, withinclustergenes, options):
    pksnrpsvars.nrpspkstypedict = {}
    for feature in withinclustergenes:
        k = utils.get_gene_id(feature)
        if not pksnrpsvars.domaindict.has_key(k):
            continue
        if pksnrpsvars.domaindict[k] == []:
            continue
        #structure of domaindict: domaindict[genename] = [[name,start,end,evalue,score],[name,start,end,evalue,score], etc.]
        domainlist = []
        nrKSdomains = 0
        for i in pksnrpsvars.domaindict[k]:
            domainlist.append(i[0])
            if i[0] == "PKS_KS":
                nrKSdomains += 1
        modKSscore = 0
        traKSscore = 0
        eneKSscore = 0
        iterKSscore = 0
        if pksnrpsvars.ksdomaindict.has_key(k):
            for i in pksnrpsvars.ksdomaindict[k]:
                if i[0] == "Trans-AT-KS":
                    traKSscore += 1
                if i[0] == "Modular-KS":
                    modKSscore += 1
                if i[0] == "Enediyne-KS":
                    eneKSscore += 1
                if i[0] == "Iterative-KS":
                    iterKSscore += 1
        if pksnrpsvars.domaindict.has_key(k):
            for i in pksnrpsvars.domaindict[k]:
                if "Cglyc" in domainlist and "Epimerization" in domainlist and "AMP-binding" in domainlist and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
                    nrpspkstype = "Glycopeptide NRPS"
                elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist) and "AMP-binding" in domainlist and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
                    nrpspkstype = "NRPS"
                elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist) or "AMP-binding" in domainlist and ("PKS_KS" in domainlist or "PKS_AT" in domainlist):
                    nrpspkstype = "Hybrid PKS-NRPS"
                elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" not in domainlist and "Trans-AT_docking" in domainlist and traKSscore > modKSscore and traKSscore > iterKSscore and traKSscore > eneKSscore:
                    nrpspkstype = "Type I Trans-AT PKS"
                elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist and iterKSscore > modKSscore and iterKSscore > traKSscore and iterKSscore > eneKSscore and nrKSdomains < 3:
                    nrpspkstype = "Type I Iterative PKS"
                elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist and eneKSscore > modKSscore and eneKSscore > traKSscore and eneKSscore > iterKSscore and nrKSdomains < 3:
                    nrpspkstype = "Type I Enediyne PKS"
                elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist and ((modKSscore > eneKSscore and modKSscore > traKSscore and modKSscore > iterKSscore) or nrKSdomains > 3):
                    nrpspkstype = "Type I Modular PKS"
                elif ("Condensation_LCL" not in domainlist and "Condensation_DCL" not in domainlist and "Condensation_Starter" not in domainlist and "Cglyc" not in domainlist and "Condensation_Dual" not in domainlist and "AMP-binding" not in domainlist) and "PKS_KS" in domainlist and "PKS_AT" in domainlist:
                    nrpspkstype = "PKS-like protein"
                elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist or "AMP-binding" in domainlist) and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
                    nrpspkstype = "NRPS-like protein"
                else:
                    nrpspkstype = "PKS/NRPS-like protein"
            if feature.qualifiers.has_key("sec_met"):
                feature.qualifiers['sec_met'].append("NRPS/PKS subtype: " + nrpspkstype)
            else:
                feature.qualifiers['sec_met'] = ["NRPS/PKS subtype: " + nrpspkstype]
            
            
            #Write motifs to seq_record
            motifFeatures = []
            if pksnrpsvars.motifdict.has_key(k):
                motifs = pksnrpsvars.motifdict[k]
                counter = 1
                for motif in motifs:
                    if feature.location.strand == 1:
                        start = feature.location.start + ( 3 * motif[1] )
                        end = feature.location.start + ( 3* motif[2] )
                    else:
                        end = feature.location.end - ( 3 * motif[1] )
                        start = feature.location.end - ( 3 * motif[2] )
                    loc = FeatureLocation(start, end, strand=feature.strand)
                    motifFeature = SeqFeature(loc, type=options.FeatureTags.pksnrpsmotifs_tag)
                    quals = defaultdict(list)
                    
                    quals['label'].append(str(motif[0]))
                    if feature.qualifiers.has_key('locus_tag'):       
                        quals['locus_tag'] = feature.qualifiers['locus_tag']
                    else:
                        quals['locus_tag'] = [k]
                    quals['motif'] = [motif[0]]
                    quals['asDomain_id'] = ['nrpspksmotif_'+'_'.join(quals['locus_tag'])+'_'+'{:04d}'.format(counter)]
                    counter += 1

                    quals['evalue'] = [str("{:.2E}".format(float(motif[3])))]
                    quals['score'] = [str(motif[4])]
                    quals['aSTool'] = ["pksnrpsmotif"]
                    quals['detection'] = ["hmmscan"]
                    quals['database'] = ["abmotifs"]
                    if feature.qualifiers.has_key('transl_table'):
                        [transl_table] = feature.qualifiers['transl_table']
                    else:
                        transl_table = 1
                    quals['translation'] = [str(motifFeature.extract(seq_record).seq.translate(table=transl_table))]
        
                    quals['note'].append("NRPS/PKS Motif: " + motif[0] + " (e-value: " + str(motif[3]) + ", bit-score: " + str(motif[4]) + ")") 
        
                    motifFeature.qualifiers = quals
                    
                    motifFeatures.append(motifFeature)
            nrpspksdomains = pksnrpsvars.domaindict[k]

            for domain in nrpspksdomains:
                if feature.qualifiers.has_key("sec_met"):
                    feature.qualifiers['sec_met'].append("NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain[0], str(domain[1]), str(domain[2]), str(domain[3]), str(domain[4])))
                else:
                    feature.qualifiers['sec_met'] = ["NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain[0], str(domain[1]), str(domain[2]), str(domain[3]), str(domain[4]))]
                

                
        seq_record.features.extend(motifFeatures)
        pksnrpsvars.nrpspkstypedict[k] = nrpspkstype
