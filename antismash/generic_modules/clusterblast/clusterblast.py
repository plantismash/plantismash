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
import sys
import os
import time
from os import path
from antismash import utils
from multiprocessing import Process
from antismash import config
from helperlibs.wrappers.io import TemporaryDirectory

def runblast(query, target):
    command = ["blastp", "-db", target, "-query", query, "-outfmt", "6", "-max_target_seqs", "10000", "-evalue", "1e-05", "-out", query.split(".")[0] + ".out"]
    utils.execute(command)


def run_diamond(query, target, tempdir, options):
    command = [
        "diamond", "blastp",
        "--db", target,
        "--threads", str(options.cpus),
        "--query", query,
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--daa", "matches.daa",
        "--tmpdir", tempdir
    ]
    return utils.execute(command)


def convert_to_tabular(tempdir):
    daa_path = path.join(tempdir, "matches.daa")
    output_path = path.join(tempdir, "input.out")

    if not path.exists(daa_path) or os.stat(daa_path).st_size == 0:
        logging.error("Error: DIAMOND did not generate a valid .daa file. Skipping conversion.")
        return None, None, 1  # Simulate a failed process return code

    command = [
        "diamond", "view",
        "-a", daa_path,
        "-o", output_path
    ]
    return utils.execute(command)



def make_blastdb(inputfile, dbname):
    command = ["makeblastdb", "-in", inputfile, "-out", dbname, "-dbtype", "prot"]
    utils.execute(command)

def load_geneclusters(searchtype):
    #Load gene cluster database into memory
    options = config.get_config()
    if options.clusterblastdir =="":
        options.clusterblastdir = path.dirname(utils.get_full_path(__file__, ''))
        options.subclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'knownclusterblast')
    else:
        options.subclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'knownclusterblast')

    if searchtype == "general" and options.taxon == "plants":
        logging.info("ClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.clusterblastdir, "plantgeneclusters.txt")
    elif searchtype == "general":
        logging.info("ClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.clusterblastdir, "geneclusters.txt")
    elif searchtype == "subclusters":
        logging.info("SubClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.subclusterblastdir, "subclusters.txt")
    elif searchtype == "knownclusters":
        logging.info("KnownClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.knownclusterblastdir, "knownclusters.txt")
    geneclustersfile = open(geneclustersfile,"r")
    filetext = geneclustersfile.read()
    lines = [line for line in filetext.split("\n") if "\t" in line]
    clusters = {}
    for i in lines:
        tabs = i.split("\t")
        accession = tabs[0]
        clusterdescription = tabs[1]
        clusternr = tabs[2]
        clustertype = tabs[3]
        clustername = accession + "_" + clusternr
        clustertags = tabs[4].split(";")
        clusterprots = tabs[5].split(";")
        clusters[clustername] = [clusterprots,clusterdescription,clustertype,clustertags]
    return clusters

def load_geneclusterproteins(accessiondict, searchtype):
    options = config.get_config()
    if options.clusterblastdir =="":
        options.clusterblastdir = path.dirname(utils.get_full_path(__file__, ''))
        options.subclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'knownclusterblast')
    else:
        options.subclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'knownclusterblast')
    #Load gene cluster database proteins info into memory
    if searchtype == "general"  and options.taxon == "plants":
        logging.info("ClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.clusterblastdir, "plantgeneclusterprots.fasta")
    elif searchtype == "general":
        logging.info("ClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.clusterblastdir, "geneclusterprots.fasta")
    elif searchtype == "subclusters":
        logging.info("SubClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.subclusterblastdir, "subclusterprots.fasta")
    elif searchtype == "knownclusters":
        logging.info("KnownClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.knownclusterblastdir, "knownclusterprots.fasta")
    gclusterprotsfile = open(gclusterprotsfile,"r")
    filetext = gclusterprotsfile.read()
    filetext = filetext.replace("\r","\n")
    lines = filetext.split("\n")
    proteinlocations = {}
    proteinstrands = {}
    proteinannotations = {}
    proteintags = {}
    for i in lines:
        if len(i) > 0 and i[0] == ">":
            tabs = i.split("|")
            protein = tabs[6]
            locustag = tabs[4]
            if locustag in accessiondict:
                locustag = "h_" + locustag
            proteintags[protein] = locustag
            location = tabs[2]
            proteinlocations[protein] = location
            strand = tabs[3]
            proteinstrands[protein] = strand
            annotation = tabs[5]
            proteinannotations[protein] = annotation
    return proteinlocations, proteinstrands, proteinannotations, proteintags

def load_clusterblast_database(seq_record, searchtype="general"):
    options = config.get_config()
    accessiondict = {}
    for cds in utils.get_cds_features(seq_record):
        accessiondict[utils.get_gene_acc(cds)] = utils.get_gene_accession(cds)
    clusters = load_geneclusters(searchtype)
    proteinlocations, proteinstrands, proteinannotations, proteintags = load_geneclusterproteins(accessiondict, searchtype)
    return clusters, proteinlocations, proteinstrands, proteinannotations, proteintags

def find_overlapping_groups(cdsfeatures):
    #Identify groups of genes with overlaps
    overlapping_groups = []
    for cdsfeature in cdsfeatures:
        overlaps = False
        for othercdsfeature in cdsfeatures:
            if utils.features_overlap(cdsfeature, othercdsfeature) and (cdsfeature.strand == othercdsfeature.strand):
                # consider only overlaps on the same strand
                overlaps = True
                added = False
                overlapping_groups2 = []
                for group in overlapping_groups:
                    if othercdsfeature in group:
                        group.append(cdsfeature)
                    overlapping_groups2.append(group)
                overlapping_groups = overlapping_groups2
                if not added:
                    overlapping_groups.append([cdsfeature, othercdsfeature])
                    added = True
                break
        if not overlaps:
            overlapping_groups.append([cdsfeature])
    return overlapping_groups

def filter_overlap(cdsfeatures):
    #For groups of overlapping CDSs (e.g., alternative transcripts?), only use the longest one
    uniquecdsfeatures = []
    overlapping_groups = find_overlapping_groups(cdsfeatures)
    for group in overlapping_groups:
        lengths = [len(utils.get_aa_sequence(feature)) for feature in group]
        longest_idx = lengths.index(max(lengths))
        uniquecdsfeatures.append(group[longest_idx])
    return uniquecdsfeatures

def create_blast_inputs(genecluster, seq_record):
    options = config.get_config()
    #Create input fasta files for BLAST search
    if options.taxon == "plants":
        queryclusterprots = filter_overlap(utils.get_cluster_cds_features(genecluster, seq_record))
    else:
        queryclusterprots = utils.get_cluster_cds_features(genecluster, seq_record)
    queryclusternames = []
    queryclusterseqs = []
    queryclusterprotsnames = []
    for cds in queryclusterprots:
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        fullname = "|".join(["input", "c" + str(utils.get_cluster_number(genecluster)), \
                             str(cds.location.start).replace(">","").replace("<","") + "-" + \
                             str(cds.location.end).replace(">","").replace("<",""), \
                             strand, utils.get_gene_acc(cds), utils.get_gene_annotation(cds)])
        if fullname not in queryclusternames:
            queryclusterseqs.append(str(utils.get_aa_sequence(cds)))
            queryclusternames.append(fullname)
            queryclusterprotsnames.append(utils.get_gene_acc(cds))

    return queryclusternames, queryclusterseqs, queryclusterprotsnames

def run_internal_blastsearch():
    #Run and parse BLAST search
    make_blastdb("internal_input.fasta", "internal_input.fasta")
    runblast("internal_input.fasta", "internal_input.fasta")
    blastoutput = open("internal_input.out","r").read()
    return blastoutput

def uniqueblasthitfilter(blastlines):
    #Filter for best blast hits (of one query on each subject)
    query_subject_combinations = []
    blastlines2 = []
    for i in blastlines:
        tabs = i.split("\t")
        query = tabs[0]
        subject = tabs[1]
        query_subject_combination = query + "_" + subject
        if query_subject_combination in query_subject_combinations:
            pass
        else:
            query_subject_combinations.append(query_subject_combination)
            blastlines2.append(i)
    return blastlines2

def tresholdblasthitfilter(blastlines, minseqcoverage, minpercidentity, seqlengths, seq_record):
    #Filters blastlines to get rid of hits that do not meet criteria
    blastlines2 = []
    for i in blastlines:
        tabs = i.split("\t")
        query = tabs[0]
        perc_ident = int(float(tabs[2]) + 0.5)
        alignmentlength = float(tabs[3])
        if query.split("|")[4] in seqlengths:
            perc_coverage = (float(tabs[3]) / seqlengths[query.split("|")[4]]) * 100
        else:
            feature_by_id = utils.get_feature_dict_protein_id(seq_record)
            seqlength = len(utils.get_aa_sequence(feature_by_id[query.split("|")[4]]))
            perc_coverage = (float(tabs[3]) / seqlength) * 100
        if perc_ident > minpercidentity and (perc_coverage > minseqcoverage):
            blastlines2.append(i)
    return blastlines2

def blastparse(blasttext, minseqcoverage, minpercidentity, seqlengths, seq_record):
    options = config.get_config()
    geneclustergenes = [utils.get_gene_acc(cds) for cds in utils.get_withincluster_cds_features(seq_record)]
    blastdict = {}
    querylist = []
    hitclusters = []
    blastlines = blasttext.split("\n")[:-1]
    blastlines = uniqueblasthitfilter(blastlines)
    blastlines = tresholdblasthitfilter(blastlines, minseqcoverage, minpercidentity, seqlengths, seq_record)
    #Goes through the blastlines. For each query, creates a querydict and hitlist, and adds these to the blastdict when finding the next query
    firstquery = "y"
    percid_per_cluster = {}
    for i in blastlines:
        tabs = i.split("\t")
        query = tabs[0]
        subject = tabs[1].split("|")[4]
        if subject == "no_locus_tag":
            subject = tabs[1].split("|")[6]
        if subject in geneclustergenes:
            subject = "h_" + subject
        if len(tabs[1].split("|")) > 6:
            locustag = tabs[1].split("|")[6]
        else:
            locustag = ""
        subject_genecluster = tabs[1].split("|")[0] + "_" + tabs[1].split("|")[1]
        subject_start = (tabs[1].split("|")[2]).split("-")[0]
        subject_end = (tabs[1].split("|")[2]).split("-")[1]
        subject_strand  = tabs[1].split("|")[3]
        subject_annotation = tabs[1].split("|")[5]
        perc_ident = int(float(tabs[2]) + 0.5)
        evalue = str(tabs[10])
        blastscore = int(float(tabs[11])+0.5)
        if query.split("|")[4] in seqlengths:
            perc_coverage = (float(tabs[3]) / seqlengths[query.split("|")[4]]) * 100
        else:
            feature_by_id = utils.get_feature_dict_protein_id(seq_record)
            seqlength = len(utils.get_aa_sequence(feature_by_id[query.split("|")[4]]))
            perc_coverage = (float(tabs[3]) / seqlength) * 100
        if firstquery == "y": #Only until the first blastline with good hit
            firstquery = "n"
            querylist.append(query)
            subjectlist = []
            querydict = {}
            subjectlist.append(subject)
            querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
            if subject_genecluster not in hitclusters:
                percid_per_cluster[subject_genecluster] = [perc_ident]
                hitclusters.append(subject_genecluster)
            last_query = query
        elif i == blastlines[-1]: #Only for the last blastline
            if query not in querylist:
                subjectlist = []
                querydict = {}
                subjectlist.append(subject)
                querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
                blastdict[query] = [subjectlist,querydict]
                querylist.append(query)
                if subject_genecluster not in hitclusters:
                    hitclusters.append(subject_genecluster)
                    percid_per_cluster[subject_genecluster] = [perc_ident]
                else:
                    percid_per_cluster[subject_genecluster].append(perc_ident)
            else:
                subjectlist.append(subject)
                querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
                blastdict[query] = [subjectlist,querydict]
                if subject_genecluster not in hitclusters:
                    hitclusters.append(subject_genecluster)
                    percid_per_cluster[subject_genecluster] = [perc_ident]
                else:
                    percid_per_cluster[subject_genecluster].append(perc_ident)
        else: #For all but the first and last blastlines
            if query not in querylist:
                blastdict[last_query] = [subjectlist,querydict]
                querylist.append(query)
                subjectlist = []
                querydict = {}
                subjectlist.append(subject)
                querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
                if subject_genecluster not in hitclusters:
                    hitclusters.append(subject_genecluster)
                    percid_per_cluster[subject_genecluster] = [perc_ident]
                else:
                    percid_per_cluster[subject_genecluster].append(perc_ident)
                last_query = query
            else:
                subjectlist.append(subject)
                querydict[subject] = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
                if subject_genecluster not in hitclusters:
                    hitclusters.append(subject_genecluster)
                    percid_per_cluster[subject_genecluster] = [perc_ident]
                else:
                    percid_per_cluster[subject_genecluster].append(perc_ident)
    #For plants, filter hitclusters to only keep those hits with at least one hit > 60% ID
    if options.taxon == "plants":
        hitclusters = [cluster for cluster in hitclusters if len([int(pid) for pid in percid_per_cluster[cluster] if int(pid) > 60]) > 0]
    return [blastdict,querylist,hitclusters]

def fastaseqlengths(seq_record):
    seqlengths = {}
    cdsfeatures = utils.get_cds_features(seq_record)
    for cds in cdsfeatures:
        seqlength = len(str(utils.get_aa_sequence(cds)))
        seqlengths[utils.get_gene_acc(cds)] = seqlength
    return seqlengths

def parse_blast(blastoutput, seq_record, minseqcoverage, minpercidentity):
    seqlengths = fastaseqlengths(seq_record)
    blastinfo = blastparse(blastoutput, minseqcoverage,
                            minpercidentity, seqlengths, seq_record)
    blastdict = blastinfo[0]
    querylist = blastinfo[1]
    hitclusters = blastinfo[2]

    return blastdict, querylist, hitclusters

def find_internal_orthologous_groups(internalhomologygroupsdict, iblastdict, iqueryclusternames, clusternumber):
    #find and store internal homologs
    groups = []
    for j in iqueryclusternames:
        if j in iblastdict:
            hits = iblastdict[j][0]
            group = []
            for k in hits:
                if k[:2] == "h_":
                    group.append(k[2:])
                elif k.count("|") > 4:
                    group.append(k.split("|")[4])
                else:
                    group.append(k)
            if j.split("|")[4] not in group:
                group.append(j.split("|")[4])
            x = 0
            for l in groups:
                for m in group:
                    if m in l:
                        del groups[x]
                        for n in l:
                            if n not in group:
                                group.append(n)
                        break
                x += 1
            group.sort()
            groups.append(group)
        else:
            groups.append([j.split("|")[4]])
    internalhomologygroupsdict[clusternumber] = groups
    return internalhomologygroupsdict

def internal_homology_blast(seq_record):
    options = config.get_config()
    #Run BLAST on gene cluster proteins of each cluster on itself to find internal homologs, store groups of homologs - including singles - in a dictionary as a list of lists accordingly
    with TemporaryDirectory(change=True):
        logging.info("Finding internal homologs in each gene cluster..")
        internalhomologygroupsdict = {}
        geneclusters = utils.get_sorted_cluster_features(seq_record)
        for genecluster in geneclusters:
            clusternumber = utils.get_cluster_number(genecluster)
            iqueryclusternames, iqueryclusterseqs, iqueryclusterprots = create_blast_inputs(genecluster, seq_record)
            utils.writefasta(iqueryclusternames, iqueryclusterseqs, "internal_input.fasta")
            blastoutput = run_internal_blastsearch()
            iblastdict, iquerylist, ihitclusters = parse_blast(blastoutput, seq_record, 25, 30)
            internalhomologygroupsdict = find_internal_orthologous_groups(internalhomologygroupsdict, iblastdict, iqueryclusternames, clusternumber)
    return internalhomologygroupsdict

def write_clusterblast_inputfiles(options, queryclusternames, queryclusterseqs):
    equalpartsizes = int(len(queryclusternames) / options.cpus)
    for i in range(options.cpus):
        if i == 0:
            setnames = queryclusternames[:equalpartsizes]
            setseqs = queryclusterseqs[:equalpartsizes]
        elif i == (options.cpus - 1):
            setnames = queryclusternames[(i*equalpartsizes):]
            setseqs = queryclusterseqs[(i*equalpartsizes):]
        else:
            setnames = queryclusternames[(i*equalpartsizes):((i+1)*equalpartsizes)]
            setseqs = queryclusterseqs[(i*equalpartsizes):((i+1)*equalpartsizes)]
        utils.writefasta(setnames, setseqs, "input" + str(i) + ".fasta")

def run_clusterblast_processes(options, searchtype="general"):

    processes = []
    for i in range(options.cpus):
        if searchtype == "general":
            blastdb = path.join(options.clusterblastdir, 'geneclusterprots.fasta')
        elif searchtype == "subclusters":
            blastdb = path.join(options.subclusterblastdir, 'subclusterprots.fasta')
        elif searchtype == "knownclusters":
            blastdb = path.join(options.knownclusterblastdir, 'knownclusterprots.fasta')
        if sys.platform == 'win32':
            blastpath, _, blastdb = blastdb.rpartition(os.sep)
            os.environ['BLASTDB'] = blastpath
        processes.append(Process(target=runblast, args=["input" + str(i) + ".fasta", blastdb]))
    for i in processes:
        i.start()
    time.sleep(10)
    while True:
        processrunning = "n"
        for i in processes:
            if i.is_alive():
                processrunning = "y"
        if processrunning == "y":
            time.sleep(5)
        else:
            break
    for i in processes:
        i.join()

def read_clusterblast_output(options):
    blastoutput = ""
    for i in range(options.cpus):
        fh = open("input" + str(i) + ".out","r")
        output = fh.read()
        fh.close()
        blastoutput = blastoutput + output
    return blastoutput

def write_raw_clusterblastoutput(outputfoldername, blastoutputs, searchtype="general"):
    if searchtype == "general":
        blastoutputfile = open(outputfoldername + os.sep + "clusterblastoutput.txt","a")
    elif searchtype == "subclusters":
        blastoutputfile = open(outputfoldername + os.sep + "subclusterblastoutput.txt","a")
    elif searchtype == "knownclusters":
        blastoutputfile = open(outputfoldername + os.sep + "knownclusterblastoutput.txt","a")
    for blastoutput in blastoutputs:
        blastoutputfile.write(blastoutput)
    blastoutputfile.close() # now all rows are written

def remove_queries_without_hits(querylist, blastdict):
    #Remove queries without hits
    querylist2 = []
    for i in querylist:
        if i in blastdict:
            querylist2.append(i)
        else:
            pass
    querylist = querylist2
    return querylist

def parse_clusterblast_dict(blastdict, querylist, clusters, hitclusternumber, hitclusterdata, allcoregenes):
    hitclusterdatalist = []
    nrhits = float(0)
    nrcoregenehits = float(0)
    cumblastscore = float(0)
    hitpositions = []
    hitposcorelist = []
    for j in querylist:
        querynrhits = 0
        querycumblastscore = float(0)
        nrhitsplus = "n"
        for k in blastdict[j][0]:
            if hitclusternumber == blastdict[j][1][k][0]:
                if blastdict[j][1][k][9] in clusters[hitclusternumber][0] and [querylist.index(j),clusters[hitclusternumber][0].index(blastdict[j][1][k][9])] not in hitpositions:
                    nrhitsplus = "y"
                    querynrhits += 1
                    blastscore = float(blastdict[j][1][k][6]) / 100000000
                    querycumblastscore = querycumblastscore + blastscore
                    hitclusterdatalist.append([j,k,blastdict[j][1][k][5],blastdict[j][1][k][6],blastdict[j][1][k][7],blastdict[j][1][k][8]])
                    hitclusterdata[hitclusternumber] = hitclusterdatalist
                    hitpositions.append([querylist.index(j),clusters[hitclusternumber][0].index(blastdict[j][1][k][9])])
        if nrhitsplus == "y":
            nrhits += 1
            if j.split("|")[4] in allcoregenes:
                nrcoregenehits += 0.01
                for hit in range(querynrhits):
                    hitposcorelist.append(1)
            else:
                for hit in range(querynrhits):
                    hitposcorelist.append(0)
            cumblastscore = cumblastscore + float(querycumblastscore)
    return hitclusterdata, nrhits, nrcoregenehits, cumblastscore, hitpositions, hitposcorelist

def find_clusterblast_hitsgroups(hitpositions):
    #Find groups of hits
    hitgroupsdict = {}
    for p in hitpositions:
        if p[0] not in hitgroupsdict:
            hitgroupsdict[p[0]] = [p[1]]
        else:
            hitgroupsdict[p[0]].append(p[1])
    return hitgroupsdict

def calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist, nrhits):
    query_givenscores_querydict = {}
    query_givenscores_hitdict = {}
    #Calculate synteny score; give score only if more than one hits (otherwise no synteny possible), and only once for every query gene and every hit gene
    synteny_score = 0
    z = 1
    if nrhits > 1:
        for p in hitpositions[:-1]:
            tandem = "n"
            #Check if a gene homologous to this gene has already been scored for synteny in the previous entry
            if p[1] in hitgroupsdict[hitpositions[z][0]]:
                tandem = "y"
            #Score entry
            if ((p[0] not in query_givenscores_querydict) or query_givenscores_querydict[p[0]] == 0) and ((p[1] not in query_givenscores_hitdict) or query_givenscores_hitdict[p[1]] == 0) and tandem == "n":
                q = hitpositions[z]
                if (abs(p[0] - q[0]) < 2) and abs(p[0]-q[0]) == abs(p[1]-q[1]):
                    synteny_score += 1
                    if hitposcorelist[z - 1] == 1 or hitposcorelist[z] == 1:
                        synteny_score += 1
                    query_givenscores_querydict[p[0]] = 1
                    query_givenscores_hitdict[p[1]] = 1
                else:
                    query_givenscores_querydict[p[0]] = 0
                    query_givenscores_hitdict[p[1]] = 0
            z += 1
    return synteny_score

def score_clusterblast_output(blastdict, querylist, hitclusters, clusters, allcoregenes):
    #Score BLAST output on all gene clusters
    #Rank gene cluster hits based on 1) number of protein hits covering >25% sequence length or at least 100aa alignment, with >30% identity and 2) cumulative blast score
    #Find number of protein hits and cumulative blast score for each gene cluster
    logging.info("   Scoring DIAMOND outputs on database of gene clusters...")
    hitclusterdict = {}
    hitclusterdata = {}
    for i in hitclusters:
        hitclusterdata, nrhits, nrcoregenehits, cumblastscore, hitpositions, hitposcorelist = parse_clusterblast_dict(blastdict, querylist, clusters, i, hitclusterdata, allcoregenes)
        hitgroupsdict = find_clusterblast_hitsgroups(hitpositions)
        synteny_score = calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist, nrhits)
        #Give bonus to gene clusters with >0 core gene hits
        if nrcoregenehits > 0:
            corebonus = 3
        else:
            corebonus = 0
        #sorting score is based on number of hits (discrete values) & cumulative blast score (behind comma values)
        sortingscore = nrhits + synteny_score + corebonus + nrcoregenehits + cumblastscore
        if len(set([pos[1] for pos in hitpositions])) > 1 and nrhits > 1:
            hitclusterdict[i] = sortingscore
    #Sort gene clusters
    rankedclusters = utils.sortdictkeysbyvaluesrev(hitclusterdict)
    rankedclustervalues = utils.sortdictkeysbyvaluesrevv(hitclusterdict)
    return rankedclusters, rankedclustervalues, hitclusterdict, hitclusterdata

def write_clusterblast_output(options, seq_record,clusterblastStorage, searchtype="general"):

    clusternumber = clusterblastStorage.clusternumber
    queryclusterprots = clusterblastStorage.queryclusterprots
    clusters = clusterblastStorage.clusters
    hitclusterdata = clusterblastStorage.hitclusterdata
    rankedclusters = clusterblastStorage.rankedclusters
    rankedclustervalues = clusterblastStorage.rankedclustervalues
    proteintags = clusterblastStorage.proteintags
    proteinlocations = clusterblastStorage.proteinlocations
    proteinannotations = clusterblastStorage.proteinannotations
    proteinstrands = clusterblastStorage.proteinstrands

    #Output for each hit: table of genes and locations of input cluster, table of genes and locations of hit cluster, table of hits between the clusters
    logging.info("   Writing output file...")
    currentdir = os.getcwd()
    if searchtype == "general":
        options.clusterblast_outputfolder = options.full_outputfolder_path + os.sep + "clusterblast"
        if not os.path.exists(options.clusterblast_outputfolder):
            os.mkdir(options.clusterblast_outputfolder)
        outputfolder = options.clusterblast_outputfolder
    elif searchtype == "subclusters":
        options.subclusterblast_outputfolder = options.full_outputfolder_path + os.sep + "subclusterblast"
        if not os.path.exists(options.subclusterblast_outputfolder):
            os.mkdir(options.subclusterblast_outputfolder)
        outputfolder = options.subclusterblast_outputfolder
    elif searchtype == "knownclusters":
        options.knownclusterblast_outputfolder = options.full_outputfolder_path + os.sep + "knownclusterblast"
        if not os.path.exists(options.knownclusterblast_outputfolder):
            os.mkdir(options.knownclusterblast_outputfolder)
        outputfolder = options.knownclusterblast_outputfolder
    os.chdir(outputfolder)
    out_file = open("cluster" + str(clusternumber) + ".txt","w")
    out_file.write("ClusterBlast scores for " + seq_record.id + "\n")
    out_file.write("\nTable of genes, locations, strands and annotations of query cluster:\n")
    feature_by_id = utils.get_feature_dict_protein_id(seq_record)
    for i in queryclusterprots:
        cds = feature_by_id[i]
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        out_file.write("\t".join([i, str(cds.location.start).replace(">","").replace("<",""), str(cds.location.end).replace(">","").replace("<",""), strand, utils.get_gene_annotation(cds)]) + "\t\n")
    out_file.write("\n\nSignificant hits: \n")
    z = 0
    for i in rankedclusters[:100]:
        out_file.write(str(z+1) + ". " + i + "\t" + clusters[i][1] + "\n")
        z += 1
    z = 0
    out_file.write("\n\nDetails:")
    for i in rankedclusters[:100]:
        value = "%.8f" % rankedclustervalues[z]
        nrhits = value.split(".")[0]
        if nrhits and int(nrhits) > 0:
            out_file.write("\n\n>>\n")
            cumblastscore = str(int(float(value.split(".")[1][2:])))
            out_file.write("\n".join([str(z+1) + ". " + i, "Source: " + clusters[i][1], "Type: " + clusters[i][2], "Number of proteins with BLAST hits to this cluster: " + nrhits,"Cumulative BLAST score: " + cumblastscore + "\n", "Table of genes, locations, strands and annotations of subject cluster:\n"]))
            clusterproteins = clusters[i][0]
            for j in clusterproteins:
                if j in proteinlocations and j in proteinannotations and j in proteinstrands:
                    if proteintags[j] == "no_locus_tag":
                        out_file.write(j + "\t")
                    else:
                        out_file.write(proteintags[j] + "\t")
                    out_file.write("\t".join([j, proteinlocations[j].split("-")[0], proteinlocations[j].split("-")[1], proteinstrands[j], proteinannotations[j]]) + "\n")
            out_file.write("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
            if i in list(hitclusterdata.keys()):
                tabledata = hitclusterdata[i]
                for x in tabledata:
                    w = 0
                    for y in x:
                        if w == 0:
                            out_file.write(str(y).split("|")[4] + "\t")
                            w += 1
                        else:
                            out_file.write(str(y) + "\t")
                    out_file.write("\n")
            else:
                out_file.write("data not found\n")
            out_file.write("\n")
            z += 1
    out_file.close()
    os.chdir(currentdir)

def perform_clusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    logging.info("Running DIAMOND gene cluster searches..")
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    with TemporaryDirectory(change=True) as tempdir:
        blastoutputs = []
        for genecluster in geneclusters:
            clusternumber = utils.get_cluster_number(genecluster)
            if options.debug and os.path.exists(options.dbgclusterblast + os.sep + "clusterblast" + os.sep + "cluster" + str(clusternumber) + ".txt"):
                logging.debug ("Skipping Clusterblast calculations, using results from %s instead" % options.dbgclusterblast + os.sep + "clusterblast"  + os.sep + "cluster" + str(clusternumber) + ".txt")
            else:

                logging.info("   Gene cluster " + str(clusternumber))
                queryclusternames, queryclusterseqs, queryclusterprots = create_blast_inputs(genecluster, seq_record)
                utils.writefasta(queryclusternames, queryclusterseqs, "input.fasta")
                if options.taxon == "plants":
                    out, err, retcode = run_diamond("input.fasta", path.join(options.clusterblastdir, "plantgeneclusterprots"), tempdir, options)
                else:
                    out, err, retcode = run_diamond("input.fasta", path.join(options.clusterblastdir, "geneclusterprots"), tempdir, options)
                if retcode != 0:
                    logging.error("Running diamond failed: returned %s, stderr: %r, stdout: %r", retcode, err, out)
                out, err, retcode = convert_to_tabular(tempdir)
                if retcode != 0:
                    logging.error("Converting daa failed: returned %s, stderr: %r, stdout: %r", retcode, err, out)

                with open("input.out", 'r') as fh:
                    blastoutput = fh.read()
                blastoutputs.append(blastoutput)
                minseqcoverage = 10
                minpercidentity = 30
                blastdict, querylist, hitclusters = parse_blast(blastoutput, seq_record, minseqcoverage, minpercidentity)
                querylist = remove_queries_without_hits(querylist, blastdict)
                allcoregenes = [utils.get_gene_acc(cds) for cds in utils.get_secmet_cds_features(seq_record)]
                rankedclusters, rankedclustervalues, hitclusterdict, hitclusterdata = score_clusterblast_output(blastdict, querylist, hitclusters, clusters, allcoregenes)

                # store all clusterblast related data in a utils.Storage object
                clusterblastStorage = utils.Storage()
                clusterblastStorage.clusternumber = clusternumber
                clusterblastStorage.queryclusterprots = queryclusterprots
                clusterblastStorage.clusters = clusters
                clusterblastStorage.hitclusterdata = hitclusterdata
                clusterblastStorage.rankedclusters = rankedclusters
                clusterblastStorage.rankedclustervalues = rankedclustervalues
                clusterblastStorage.proteintags = proteintags
                clusterblastStorage.proteinlocations = proteinlocations
                clusterblastStorage.proteinannotations = proteinannotations
                clusterblastStorage.proteinstrands = proteinstrands


                #write_clusterblast_output(options, seq_record, clusternumber, queryclusterprots, clusters, hitclusterdata, rankedclusters, rankedclustervalues, proteintags, proteinlocations, proteinannotations, proteinstrands)
                write_clusterblast_output(options, seq_record, clusterblastStorage)

        write_raw_clusterblastoutput(options.full_outputfolder_path, blastoutputs)
        logging.info("   DIAMOND search finished. Parsing results...")