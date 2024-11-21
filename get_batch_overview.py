#!/usr/bin/env python

import run_antismash
from antismash import config
from antismash import utils
from antismash.generic_modules import hmm_detection

import sys
import logging
import json
import multiprocessing
import numpy
from argparse import Namespace
from os import path

def get_mockup_config():
    options = Namespace()
    config.load_config(options)
    options.debug = True
    options.verbose = True
    options.fix_id_line = False
    options.input_type = "nucl"
    options.limit = -1
    options.gff3 = False
    options.hmmsearch_chunk = 10000
    options.cpus = multiprocessing.cpu_count()
    options.start = -1
    options.end = -1
    options.cutoff_multiplier = 1.00
    options.gene_num_cutoff = 0
    options.gene_num_cutoff_only = False
    options.cdh_use_binary = True
    options.genefinding = "glimmer"
    options.all_orfs = False
    options.glimmerhmm_train_folder = "arabidopsis"
    options.taxon = "plants"
    options.enabled_detection_models = ["plants"]
    options.enabled_cluster_types = hmm_detection.get_supported_cluster_types()
    options.outputfoldername = path.dirname(__file__)
    options.full_outputfolder_path = path.abspath(options.outputfoldername)
    options.full_hmmer = False
    options.ecpred = 'none'
    options.run_asf = False
    temp_ct = []
    for cl_type in options.enabled_cluster_types:
        if len(cl_type.split("/")) > 1:
            if cl_type.split("/")[0] in options.enabled_detection_models:
                temp_ct.append(cl_type)
        elif "default" in options.enabled_detection_models :
            temp_ct.append(cl_type)
    options.enabled_cluster_types = temp_ct
    run_antismash.apply_taxon_preset(options)
    return options

def main():
    multiprocessing.freeze_support()
    res_object = {}

    # get genome files
    files = []
    for line in open(sys.argv[1], 'r'):
        files.append(path.expanduser(line.replace("\n", "")))

    # mockup antismash run per files
    i = 1
    for fpath in files:
        res_object[fpath] = {}
        print "Processing %s... (%d/%d)" % (fpath, i, len(files));
        i += 1
        options = get_mockup_config()
        options.sequences = [fpath]
        config.set_config(options)
        run_antismash.setup_logging(options) #To-DO: get antismash logging to works!

        # load plugins
        plugins = run_antismash.load_detection_plugins()
        run_antismash.filter_plugins(plugins, options, options.enabled_cluster_types)

        # parse to seq_records
        seq_records = run_antismash.parse_input_sequences(options)
        options.next_clusternr = 1

        for seq_record in seq_records:
            if options.input_type == 'nucl':
                seq_records = [record for record in seq_records if len(record.seq) > 1000]
                if len(seq_records) < 1:
                    continue
            utils.sort_features(seq_record)
            run_antismash.strip_record(seq_record)
            utils.fix_record_name_id(seq_record, options)

            # fetch results_by_id
            feature_by_id = utils.get_feature_dict(seq_record)
            results = []
            results_by_id = {}
            for feature in utils.get_cds_features(seq_record):
                prefix = "%s:" % seq_record.id.replace(":", "_")
                gene_id = utils.get_gene_id(feature)
                if (prefix + gene_id) in options.hmm_results:
                    results_by_id[gene_id] = options.hmm_results[prefix + gene_id]
                    for res in results_by_id[gene_id]:
                        results.append(res)

            # ignore short aa's
            min_length_aa = 100
            short_cds_buffer = []
            for f in seq_record.features: # temporarily remove short aa
                if f.type == "CDS" and len(f.qualifiers['translation'][0]) < min_length_aa and not results_by_id.has_key(utils.get_gene_id(f)):
                    short_cds_buffer.append(f)
                    seq_record.features.remove(f)

            overlaps = utils.get_overlaps_table(seq_record)
            rulesdict = hmm_detection.create_rules_dict(options.enabled_cluster_types)
            # find total cdhit numbers in the chromosome
            total_cdhit = len(utils.get_cdhit_table(utils.get_cds_features(seq_record), options)[0])
            res_object[fpath][seq_record.id] = {"total_clusters" : 0, "total_genes" : len(overlaps[0]), "total_cdhit" : total_cdhit, "genes_with_hits" : 0, "largest_cdhit" : 0, "largest_domain_variations" : 0, "per_hits" : {}, "cluster_types" : {}}

            # filter overlap hits
            results, results_by_id = hmm_detection.filter_results(results, results_by_id, overlaps, feature_by_id)

            # count hits
            for gene_id in results_by_id:
                res_gene = results_by_id[gene_id]
                if len(res_gene) > 0:
                    res_object[fpath][seq_record.id]["genes_with_hits"] += 1
                for hsp in res_gene:
                    domain_name = hsp.query_id.replace("plants/", "")
                    if domain_name not in res_object[fpath][seq_record.id]["per_hits"]:
                        res_object[fpath][seq_record.id]["per_hits"][domain_name] = 0
                    res_object[fpath][seq_record.id]["per_hits"][domain_name] += 1

            # do cluster finding algorithm
            typedict = hmm_detection.apply_cluster_rules(results_by_id, feature_by_id, options.enabled_cluster_types, rulesdict, overlaps)
            hmm_detection.fix_hybrid_clusters_typedict(typedict)
            nseqdict = hmm_detection.get_nseq()
            for cds in results_by_id.keys():
                feature = feature_by_id[cds]
                if typedict[cds] != "none":
                    hmm_detection._update_sec_met_entry(feature, results_by_id[cds], typedict[cds], nseqdict)
            hmm_detection.find_clusters(seq_record, rulesdict, overlaps)
            seq_record.features.extend(short_cds_buffer)
            res_object[fpath][seq_record.id]["total_clusters"] += len(utils.get_cluster_features(seq_record))

            # do cluster specific and unspecific analysis
            if len(utils.get_cluster_features(seq_record)) > 0:
                run_antismash.cluster_specific_analysis(plugins, seq_record, options)
            run_antismash.unspecific_analysis(seq_record, options)

            #Rearrange hybrid clusters name alphabetically
            hmm_detection.fix_hybrid_clusters(seq_record)

            #before writing to output, remove all hmm_detection's subdir prefixes from clustertype
            for cluster in utils.get_cluster_features(seq_record):
                prod_names = []
                for prod in cluster.qualifiers['product']:
                    prod_name = []
                    for name in prod.split('-'):
                        prod_name.append(name.split('/')[-1])
                    prod_names.append("-".join(prod_name))
                cluster.qualifiers['product'] = prod_names
            for cds in utils.get_cds_features(seq_record):
                if 'sec_met' in cds.qualifiers:
                    temp_qual = []
                    for row in cds.qualifiers['sec_met']:
                        if row.startswith('Type: '):
                            clustertypes = [(ct.split('/')[-1]) for ct in row.split('Type: ')[-1].split('-')]
                            temp_qual.append('Type: ' + "-".join(clustertypes))
                        elif row.startswith('Domains detected: '):
                            cluster_results = []
                            for cluster_result in row.split('Domains detected: ')[-1].split(';'):
                                cluster_results.append(cluster_result.split(' (E-value')[0].split('/')[-1] + ' (E-value' + cluster_result.split(' (E-value')[-1])
                            temp_qual.append('Domains detected: ' + ";".join(cluster_results))
                        else:
                            temp_qual.append(row)
                    cds.qualifiers['sec_met'] = temp_qual

            #on plants, remove plant clustertype from hybrid types, and replace single
            #plant clustertype with "putative"
            for cluster in utils.get_cluster_features(seq_record):
                prod_names = []
                for prod in cluster.qualifiers['product']:
                    prod_name = list(set(prod.split('-')))
                    if (len(prod_name) > 1) and ("plant" in prod_name):
                        prod_name.remove("plant")
                    elif prod_name == ["plant"]:
                        prod_name = ["putative"]
                    prod_names.append("-".join(prod_name))
                cluster.qualifiers['product'] = prod_names
            for cds in utils.get_cds_features(seq_record):
                if 'sec_met' in cds.qualifiers:
                    temp_qual = []
                    for row in cds.qualifiers['sec_met']:
                        if row.startswith('Type: '):
                            clustertypes = list(set(row.split('Type: ')[-1].split('-')))
                            if (len(clustertypes) > 1) and ("plant" in clustertypes):
                                clustertypes.remove("plant")
                            elif clustertypes == ["plant"]:
                                clustertypes = ["putative"]
                            temp_qual.append('Type: ' + "-".join(clustertypes))
                        else:
                            temp_qual.append(row)
                    cds.qualifiers['sec_met'] = temp_qual

            # find largest cdhit number & largest domain diversity in a cluster
            res_object[fpath][seq_record.id]["average_cdhit"] = 0
            res_object[fpath][seq_record.id]["average_domain_variations"] = 0
            cdhit_numbers = []
            domain_numbers = []
            for cluster in utils.get_cluster_features(seq_record):
                cluster_type = utils.get_cluster_type(cluster)
                if cluster_type not in res_object[fpath][seq_record.id]["cluster_types"]:
                    res_object[fpath][seq_record.id]["cluster_types"][cluster_type] = 0
                res_object[fpath][seq_record.id]["cluster_types"][cluster_type] += 1
                num_cdhit = len(utils.get_cluster_cdhit_table(cluster, seq_record,options))
                num_domain = len(utils.get_cluster_domains(cluster, seq_record))
                cdhit_numbers.append(num_cdhit)
                domain_numbers.append(num_domain)
                if num_cdhit > res_object[fpath][seq_record.id]["largest_cdhit"]:
                    res_object[fpath][seq_record.id]["largest_cdhit"] = num_cdhit
                if num_domain > res_object[fpath][seq_record.id]["largest_domain_variations"]:
                    res_object[fpath][seq_record.id]["largest_domain_variations"] = num_domain
            if len(cdhit_numbers) > 0:
                res_object[fpath][seq_record.id]["average_cdhit"] = numpy.median(cdhit_numbers)
            if len(domain_numbers) > 0:
                res_object[fpath][seq_record.id]["average_domain_variations"] = numpy.median(domain_numbers)

        with open('result.js', 'w') as h:
            h.write('var result = %s;' % json.dumps(res_object, indent = 4))


if __name__ == "__main__":
    main()
