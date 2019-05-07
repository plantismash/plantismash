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

"""Excel and TXT output format module

"""

from pyExcelerator import *
from pyExcelerator.Workbook import *
import os
from antismash import utils
from os import path

name = "xls"
short_description = "xls_txt_output"
priority = 1

def write(seq_records, options):
    if options.input_type == 'prot':
        return
    #Open up TXT file and XLS record
    outfolder = options.full_outputfolder_path
    txtfile = open(path.join(outfolder, "geneclusters.txt"),"w")
    wb = Workbook()
    font1 = Font()
    style1 = XFStyle()
    style1.font = font1
    font1.bold = True
    ws0 = wb.add_sheet('0')
    ws0.write(0, 0, "Input accession number", style1)
    ws0.write(0, 1, "Input name", style1)
    ws0.write(0, 2, "Gene cluster type" ,style1)
    ws0.write(0, 3, "Gene cluster genes", style1)
    ws0.write(0, 4, "Gene cluster gene accessions", style1)
    if options.knownclusterblast:
      ws0.write(0, 5, "Compound with gene cluster of highest homology", style1)
    #For each gene cluster, write out info
    column = 1
    for seq_record in seq_records:
        clusters = utils.get_cluster_features(seq_record)
        for cluster in clusters:
            clustertype = utils.get_cluster_type(cluster)
            clusternr = utils.get_cluster_number(cluster)
            clustergenes = [utils.get_gene_id(cds) for cds in utils.get_cluster_cds_features(cluster, seq_record)]
            accessions = [utils.get_gene_acc(cds) for cds in utils.get_cluster_cds_features(cluster, seq_record)]
            ws0.write(column, 0, seq_record.id)
            try:
                ws0.write(column, 1, seq_record.description)
            except:
                ws0.write(column, 1, "Name to long to be contained in Excel cell; see txt file in downloadable zip archive.")
            ws0.write(column, 2, clustertype)
            try:
                ws0.write(column, 3, ";".join(clustergenes))
            except:
                ws0.write(column, 3, "Too many genes to be contained in Excel cell; see txt file in downloadable zip archive.")
            try:
                ws0.write(column, 4, ";".join(accessions))
            except:
                ws0.write(column, 4,"Too many genes to be contained in Excel cell; see txt file in downloadable zip archive.")
            if hasattr(seq_record, 'closestcompounddict') and \
               seq_record.closestcompounddict.has_key(clusternr):
                ws0.write(column, 5, seq_record.closestcompounddict[clusternr])
            column += 1
            txtfile.write("\t".join([seq_record.id, seq_record.description, clustertype, ";".join(clustergenes), ";".join(accessions)]) + "\n")
    wb.save(path.join(outfolder,  "%s.geneclusters.xls" % seq_record.id))
