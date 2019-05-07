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


import os
from antismash import utils
import itertools
import logging
from helperlibs.wrappers.io import TemporaryDirectory


def analyse_biosynthetic_order(pksnrpsvars, seq_record, options):
    #Find NRPS/PKS gene clusters
    nrpspksclusters = list(set(utils.get_cluster_features_of_type(seq_record, "nrps") + utils.get_cluster_features_of_type(seq_record, "pks")))
    #Predict biosynthetic gene order in gene cluster using starter domains, thioesterase domains, gene order and docking domains
    if not 'docking' in options:
        options.docking = {}
    for genecluster in nrpspksclusters:
        clusterpksnrpsgenes = find_clusterpksnrpsgenes(genecluster, pksnrpsvars.pksnrpscoregenes)
        if len(clusterpksnrpsgenes) > 0:
            pksgenes, clusterpksgenes, nrpsgenes, clusternrpsgenes, hybridgenes, clusterhybridgenes = find_cluster_modular_enzymes(clusterpksnrpsgenes, pksnrpsvars)
            #If more than three PKS genes, use dock_dom_analysis if possible to identify order
            if pksgenes > 3 and pksgenes < 11 and nrpsgenes == 0 and hybridgenes == 0:
                geneorder = perform_docking_domain_analysis(options, clusterpksgenes, utils.get_cluster_number(genecluster), seq_record, pksnrpsvars)
                options.docking[utils.get_cluster_number(genecluster)] = True
            else:
                geneorder = find_colinear_order(clusterpksnrpsgenes, seq_record, pksnrpsvars.domainnamesdict)
                options.docking[utils.get_cluster_number(genecluster)] = False
            generate_substrates_order(utils.get_cluster_number(genecluster), geneorder, pksnrpsvars, seq_record)

def find_clusterpksnrpsgenes(genecluster, pksnrpscoregenes):
    clusterpksnrpsgenes = []
    for gene in pksnrpscoregenes:
        if gene.location.start in xrange(genecluster.location.start, genecluster.location.end) or gene.location.end in xrange(genecluster.location.start, genecluster.location.end):
            clusterpksnrpsgenes.append(gene)
    return clusterpksnrpsgenes

def find_cluster_modular_enzymes(clusterpksnrpsgenes, pksnrpsvars):
    clusterpksnrpsgenenames = [utils.get_gene_id(feature) for feature in clusterpksnrpsgenes]
    pksgenes = 0
    clusterpksgenes = []
    nrpsgenes = 0
    clusternrpsgenes = []
    hybridgenes = 0
    clusterhybridgenes = []
    for j in clusterpksnrpsgenenames:
        k = pksnrpsvars.nrpspkstypedict[j]
        if "PKS" in k and "NRPS" not in k:
            pksgenes += 1
            clusterpksgenes.append(j)
        elif "PKS" not in k and "NRPS" in k:
            nrpsgenes += 1
            clusternrpsgenes.append(j)
        elif "PKS/NRPS" in k:
            if ("PKS_KS" in pksnrpsvars.domainnamesdict[j] or "PKS_AT" in pksnrpsvars.domainnamesdict[j]) and ("AMP-binding" not in pksnrpsvars.domainnamesdict[j] and "A-OX" not in pksnrpsvars.domainnamesdict[j] and "Condensation" not in pksnrpsvars.domainnamesdict[j]):
                pksgenes += 1
                clusterpksgenes.append(j)
            elif ("PKS_KS" not in pksnrpsvars.domainnamesdict[j] and  "PKS_AT" not in pksnrpsvars.domainnamesdict[j]) and ("AMP-binding" in pksnrpsvars.domainnamesdict[j] or "A-OX" in pksnrpsvars.domainnamesdict[j] or "Condensation" in pksnrpsvars.domainnamesdict[j]):
                nrpsgenes += 1
                clusternrpsgenes.append(j)
        elif "PKS" in k and "NRPS" in k:
            hybridgenes += 1
            clusterhybridgenes.append(j)
    return pksgenes, clusterpksgenes, nrpsgenes, clusternrpsgenes, hybridgenes, clusterhybridgenes

def generate_substrates_order(genecluster, geneorder, pksnrpsvars, seq_record):
    #Generate substrates order from predicted gene order and consensus predictions
    prediction = ""

    for f in utils.get_cluster_features(seq_record):
	cluster_info = f.qualifiers

    for k in geneorder:
        if len(prediction) == 0 or prediction[-1] != "(":
            prediction += "("
        domains = pksnrpsvars.domainnamesdict[k]
        nra = 0
        nrat = 0
        nrcal = 0
	nrtransat = 0
        domainnr = 0
        consensuspred_list=[]

        for l in domains:
            if 'transatpks' not in cluster_info['product'][0]:
                if "PKS_AT" in l:
                    if domainnr > 0:
                        prediction += "-"
                    nrat += 1
                    prediction = prediction + pksnrpsvars.consensuspreds[k + "_AT" + str(nrat)]
                    consensuspred_list.append(pksnrpsvars.consensuspreds[k + "_AT" + str(nrat)])
                    domainnr += 1
            elif 'transatpks' in cluster_info['product'][0]:
	        if "PKS_KS" in l:
		    if domainnr > 0:
                        prediction += "-"
		    nrtransat += 1
		    prediction = prediction + pksnrpsvars.consensuspreds[k + "_KS" + str(nrtransat)]
		    consensuspred_list.append(pksnrpsvars.consensuspreds[k + "_KS" + str(nrtransat)])
		    domainnr += 1
            if "AMP-binding" in l or "A-OX" in l:
                if domainnr > 0:
                    prediction += "-"
                nra += 1
                prediction = prediction + pksnrpsvars.consensuspreds[k + "_A" + str(nra)]
                consensuspred_list.append(pksnrpsvars.consensuspreds[k + "_A" + str(nra)])
                domainnr += 1
            if "CAL_domain" in l:
                if domainnr > 0:
                    prediction += "-"
                nrcal += 1
                prediction = prediction + pksnrpsvars.consensuspreds[k + "_CAL" + str(nrcal)]
                consensuspred_list.append(pksnrpsvars.consensuspreds[k + "_CAL" + str(nrcal)])
                domainnr += 1
        if pksnrpsvars.consensuspred_gene_dict.has_key(k):
            logging.warn ("WARNING: Consensus specificity prediction already defined for %s; possibly duplicate genename? Overwriting entries for %s" % (k,k))
        pksnrpsvars.consensuspred_gene_dict[k]=consensuspred_list
        if prediction[-3:] == "+ (":
            prediction = prediction[:-1]
        elif prediction[-1] != "(":
            prediction += ") + "
    prediction = prediction[:-3]
    pksnrpsvars.compound_pred_dict[genecluster] = prediction

def find_first_and_last_genes(clusterpksgenes, domainnamesdict):
    #Find first and last genes based on starter module and TE / TD
    startergene = ""
    endinggene = ""
    for k in clusterpksgenes:
        if "Thioesterase" in domainnamesdict[k] or "TD" in domainnamesdict[k]:
            if endinggene == "":
                endinggene = k
            else:
                endinggene = ""
        if len(domainnamesdict[k]) >=2 and  "PKS_AT" == domainnamesdict[k][0] and "ACP" == domainnamesdict[k][1]:
            if startergene == "":
                startergene = k
            else:
                startergene = ""
    if startergene == "":
        for k in clusterpksgenes:
            if len(domainnamesdict[k]) >=3 and "PKS_KS" == domainnamesdict[k][0] and "PKS_AT" == domainnamesdict[k][1] and "ACP" == domainnamesdict[k][2]:
                if startergene == "":
                    startergene = k
                else:
                    startergene = ""
                    break
    return startergene, endinggene

def extractpositions(fasta_file, positions, refsequencename, querysequencename):
    fasta_handle = open(fasta_file, 'r')
    fasta_lines = fasta_handle.readlines()
    fasta_handle.close()
    names_seqs = parse_fasta_string_list(fasta_lines)
    muscle_names, muscle_seqs = zip(*names_seqs)
    return extract_positions_from_refseq(muscle_seqs, muscle_names, positions,
                                         refsequencename, querysequencename)

def parse_fasta_string_list(string_list):
    """Parse a string list containing fasta sequences"""
    seq_list = []
    header = None
    seq = []
    for line in string_list:
        line = line.strip()
        # skip empty lines
        if line == "":
            continue
        if line.startswith('>'):
            if header and seq:
                seq_list.append((header, "".join(seq)))
            header = line[1:68]
            seq = []
        else:
            seq.append(line)
    if header and seq:
        seq_list.append((header, "".join(seq)))
    if len(seq_list) == 0:
        raise ParserError
    return seq_list

def extract_positions_from_refseq(muscle_seqs, muscle_names, positions,
                                  refsequencename, querysequencename):
    "Count residues in ref sequence and put positions in list"
    residues = []
    refseqnr = muscle_names.index(refsequencename)
    #Extract activity signature
    refseq = muscle_seqs[refseqnr]
    poslist = []
    wanted_pos = 0
    ref_seq_pos = 0
    while refseq != "":
        i = refseq[0]
        if ref_seq_pos in positions and i != "-":
            poslist.append(wanted_pos)
        if i != "-":
            ref_seq_pos += 1
        wanted_pos += 1
        refseq = refseq[1:]
    #Extract positions from query sequence
    query_seqnr = muscle_names.index(querysequencename)
    query_seq = muscle_seqs[query_seqnr]
    for j in poslist:
        residues.append(query_seq[j])
    return residues


def extract_nterminus(da_dir, clusterpksgenes, seq_record, startergene, feature_by_id):
    #Extract N-terminal 50 residues of each non-starting protein, scan for docking domains using hmmsearch, parse output to locate interacting residues
    ntermintresdict = {}
    ntermnames = []
    ntermseqs = []
    nterm_file = os.path.join(da_dir, 'nterm.fasta')
    for k in clusterpksgenes:
        if k != startergene:
            ntermnames.append(k)
            seq = str(utils.get_aa_sequence(feature_by_id[k]))
            ntermseqs.append(seq[:50])
    ntermfasta = "input.fasta"
    z = 0
    for k in ntermnames:
        utils.writefasta([ntermnames[z]], [ntermseqs[z]],
                   ntermfasta)
        utils.execute(["muscle", "-profile", "-quiet", "-in1", nterm_file, "-in2", "input.fasta", "-out", "muscle.fasta"])
        intresidues = extractpositions("muscle.fasta", [2, 15],
                                       "EryAIII_5_6_ref",
                                       ntermnames[z])
        ntermintresdict[ntermnames[z]] = intresidues
        z += 1
    return ntermintresdict


def extract_cterminus(da_dir, clusterpksgenes, seq_record, endinggene, feature_by_id):
    #Extract C-terminal 100 residues of each non-ending protein, scan for docking domains using hmmsearch, parse output to locate interacting residues
    ctermintresdict = {}
    ctermnames = []
    ctermseqs = []
    cterm_file = os.path.join(da_dir, 'cterm.fasta')
    for k in clusterpksgenes:
        if k != endinggene:
            ctermnames.append(k)
            seq = str(utils.get_aa_sequence(feature_by_id[k]))
            ctermseqs.append(seq[-100:])
    ctermfasta = "input.fasta"
    z = 0
    for k in ctermnames:
        utils.writefasta([ctermnames[z]], [ctermseqs[z]],
                   ctermfasta)
        utils.execute(["muscle", "-profile", "-quiet", "-in1", cterm_file, "-in2", "input.fasta", "-out", "muscle.fasta"])
        intresidues = extractpositions("muscle.fasta", [55, 64],
                                       "EryAII_ref", ctermnames[z])
        ctermintresdict[ctermnames[z]] = intresidues
        z += 1
    return ctermintresdict

def find_possible_orders(clusterpksgenes, startergene, endinggene):
    genes_to_order = []
    z = 0
    for k in clusterpksgenes:
        if k == startergene or k == endinggene:
            pass
        else:
            genes_to_order.append(k)
        z += 1
    possible_orders = list(itertools.permutations(genes_to_order,len(genes_to_order)))
    return possible_orders

def rank_biosynthetic_orders(ntermintresdict, ctermintresdict, startergene, endinggene, possible_orders):
    #If docking domains found in all, check for optimal order using interacting residues
    hydrophobic = ["A","V","I","L","F","W","Y","M"]
    positivecharge = ["H","K","R"]
    negativecharge = ["D","E"]
    possible_orders_scoredict = {}
    for k in possible_orders:
        score = 0
        interactions = []
        z = 0
        for l in k[:-1]:
            interactions.append([l,k[z + 1]])
            z += 1
        for l in interactions:
            res1a = ctermintresdict[l[0]][0]
            res1b = ntermintresdict[l[1]][0]
            res2a = ctermintresdict[l[0]][1]
            res2b = ntermintresdict[l[1]][1]
            if (res1a in hydrophobic and res1b in hydrophobic) or (res1a in positivecharge and res1b in negativecharge) or (res1a in negativecharge and res1b in positivecharge):
                score += 1
            if (res1a in positivecharge and res1b in positivecharge) or (res1a in negativecharge and res1b in negativecharge):
                score = score - 1
            if (res2a in hydrophobic and res2b in hydrophobic) or (res2a in positivecharge and res2b in negativecharge) or (res2a in negativecharge and res2b in positivecharge):
                score += 1
            if (res2a in positivecharge and res2b in positivecharge) or (res2a in negativecharge and res2b in negativecharge):
                score = score - 1
        possible_orders_scoredict[k] = score
    ranked_orders = utils.sortdictkeysbyvaluesrev(possible_orders_scoredict)
    ranked_orders_part = []
    ranked_orders2 = []
    a = 0
    for i in ranked_orders:
        if a == 0:
            score = possible_orders_scoredict[i]
            ranked_orders_part.append(i)
        elif a == (len(ranked_orders) - 1):
            ranked_orders_part.append(i)
            ranked_orders2 = ranked_orders2 + ranked_orders_part
        else:
            if possible_orders_scoredict[i] == score:
                ranked_orders_part.append(i)
            else:
                ranked_orders_part.reverse()
                ranked_orders2 = ranked_orders2 + ranked_orders_part
                score = possible_orders_scoredict[i]
                ranked_orders_part = []
                ranked_orders_part.append(i)
        a += 1
    ranked_orders = ranked_orders2[:1000]
    geneorders = ranked_orders
    geneorders2 = []
    for l in geneorders:
        geneorder = []
        if startergene != "":
            geneorder.append(startergene)
        for m in l:
            geneorder.append(m)
        if endinggene != "":
            geneorder.append(endinggene)
        geneorders2.append(geneorder)
    geneorders = geneorders2
    return geneorders, possible_orders_scoredict

def write_gene_orders_to_html(options, geneorders, possible_orders_scoredict, genecluster, startergene, endinggene):
    dockingdomain_outputfolder = options.outputfoldername + os.sep + "docking_analysis"
    if not os.path.exists(dockingdomain_outputfolder):
        os.mkdir(dockingdomain_outputfolder)
    dockhtmlfile = open(dockingdomain_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_" + str(genecluster) + ".html","w")
    if len(geneorders) == 1000:
        dockhtmlfile.write('<html>\n<head>\n<LINK href="style.css" rel="stylesheet" type="text/css">\n</head>\n<body>\nDocking domain analysis.  Score for 1000 highest scoring gene orders:<br><br><table border=1>\n')
    else:
        dockhtmlfile.write('<html>\n<head>\n<LINK href="style.css" rel="stylesheet" type="text/css">\n</head>\n<body>\nDocking domain analysis. Scores for all possible gene orders:<br><br><table border=1>\n')
    dockhtmlfile.write('<tr><td><b>Gene order</b></td><td><b>Score</b></td></tr>\n')
    for l in geneorders:
        string = "<tr><td>"
        for m in l:
            string = string + m + ","
        if startergene != "" and endinggene != "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l[1:-1])])
        elif startergene == "" and endinggene != "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l[:-1])])
        elif startergene != "" and endinggene == "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l[1:])])
        elif startergene == "" and endinggene == "":
            string = string[:-1] + "</td><td>" + str(possible_orders_scoredict[tuple(l)])
        dockhtmlfile.write(string + "</td></tr>\n")
    dockhtmlfile.write('\n</table></body></html>')
    dockhtmlfile.close()


def perform_docking_domain_analysis(options, clusterpksgenes, genecluster, seq_record, pksnrpsvars):
    feature_by_id = utils.get_feature_dict(seq_record)
    #log("Predicting PKS gene order by docking domain sequence " \
    #    "analysis", stdout=True)
    startergene, endinggene = find_first_and_last_genes(clusterpksgenes, pksnrpsvars.domainnamesdict)
    with TemporaryDirectory(change=True):
        dockinganalysis_dir = utils.get_full_path(__file__, "docking_analysis")
        ntermintresdict = extract_nterminus(dockinganalysis_dir, clusterpksgenes,
                                            seq_record, startergene, feature_by_id)
        ctermintresdict = extract_cterminus(dockinganalysis_dir, clusterpksgenes,
                                            seq_record, endinggene, feature_by_id)
    possible_orders = find_possible_orders(clusterpksgenes, startergene, endinggene)
    geneorders, possible_orders_scoredict = rank_biosynthetic_orders(ntermintresdict, ctermintresdict, startergene, endinggene, possible_orders)
    write_gene_orders_to_html(options, geneorders, possible_orders_scoredict, genecluster, startergene, endinggene)
    #log("Predicting PKS gene order by docking domain sequence " \
    #    "analysis succeeded.", stdout=True)
    #Write html outfile with docking domain analysis output
    pksnrpsvars.dockingdomainanalysis.append(genecluster)
    return geneorders[0]

def find_colinear_order(clusterpksnrpsgenes, seq_record, domainnamesdict):
    feature_by_id = utils.get_feature_dict(seq_record)
    #If NRPS genes, mixed NRPS/PKS genes, PKS genes without detected docking domains, or clusters with a 1-3 PKS genes, assume colinearity
    direction = 0
    for feature in clusterpksnrpsgenes:
        k = utils.get_gene_id(feature)
        if feature_by_id[k].strand == 1:
            direction += 1
        elif feature_by_id[k].strand == -1:
            direction = direction - 1
    if direction < 0:
        clusterpksnrpsgenes.reverse()
    #Reverse if first gene encodes a multidomain protein with a TE/TD domain
    if "Thioesterase" in domainnamesdict[utils.get_gene_id(clusterpksnrpsgenes[0])] or "TD" in domainnamesdict[utils.get_gene_id(clusterpksnrpsgenes[0])]:
        if len(domainnamesdict[utils.get_gene_id(clusterpksnrpsgenes[0])]) > 1:
            clusterpksnrpsgenes.reverse()
    geneorder = [utils.get_gene_id(feature) for feature in clusterpksnrpsgenes]
    return geneorder
