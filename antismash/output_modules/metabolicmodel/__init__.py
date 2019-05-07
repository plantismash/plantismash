# vim: set fileencoding=utf-8 :
#

#
# Copyright (C) 2014 Tilmann Weber
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# Copyright (C) 2014 Hyun Uk Kim
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
write model sbml file and lists of enzymes/reactions
"""
import logging
import os
try:
    from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
    libImport = True
except ImportError, ImportWarning:
    libImport = False
try:
    import cPickle as pickle
except ImportError:
    import pickle

name = "write_metabolicmodel"
short_description = "write_metabolicmodel"
priority = 1000

def write(seq_records, options):
    logging.debug("start writing model file module")
    if ((options.modeling == "none" ) or (libImport == False)):
        logging.debug("No metabolic models available, skipping export")
        return
    elif "modeling_successful" in options:
        if not options.modeling_successful:
            logging.error("Error in modeling step, no outputfiles generated!")
            return
    else:
        logging.debug("Model selected: %s" % options.modeling )
    
    # Set up output directory (if not already present)
    model_dirname = "metabolicModel"
    basename = options.outputfoldername
    options.metabolicmodeldir = os.path.join(basename, model_dirname)
    logging.debug("Writing metabolic models to %r" % options.metabolicmodeldir)
    if not os.path.exists(options.metabolicmodeldir):
        os.mkdir(options.metabolicmodeldir)
        
    # Model file is always stored in seq_records[0]
    seq_record = seq_records[0]
    logging.debug("Trying to retrieve model data for %s" % seq_record.id)
    if options.extrarecord.has_key(seq_record.id):
        if options.extrarecord[seq_record.id].extradata.has_key("MetabolicModelDataObj"):
            target_model = pickle.loads(options.extrarecord[seq_record.id].extradata["MetabolicModelDataObj"])
            #Output files
	    #Model reloading and overwrtting are necessary for model consistency:
	    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
	    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
	    #Cobrapy IO module seems to have an error for adding new reactions
            logging.debug("Writing model SBML file...")
            write_cobra_model_to_sbml_file(target_model, options.metabolicmodeldir+os.sep+'antiSMASH_model_with_template_%s.xml' %(options.modeling))
            target_model = create_cobra_model_from_sbml_file(options.metabolicmodeldir+os.sep+'antiSMASH_model_with_template_%s.xml' %(options.modeling))
            write_cobra_model_to_sbml_file(target_model, options.metabolicmodeldir+os.sep+'antiSMASH_model_with_template_%s.xml' %(options.modeling))
            
            fp1 = open(options.metabolicmodeldir+os.sep+'antiSMASH_model_with_template_%s_reactions.txt' % options.modeling, "w")
            fp2 = open(options.metabolicmodeldir+os.sep+'antiSMASH_model_with_template_%s_metabolites.txt' % options.modeling, "w")
            fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Lower bound"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
            fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

            #Output on screen
            logging.debug("Number of genes in final model: %s", len(target_model.genes))
	    logging.debug("Number of reactions in final model: %s", len(target_model.reactions))
	    logging.debug("Number of metabolites in final model: %s", len(target_model.metabolites))

            for j in range(len(target_model.reactions)):
                rxn = target_model.reactions[j]
                print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.lower_bound, rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)
            
            for i in range(len(target_model.metabolites)):
                metab = target_model.metabolites[i]
                print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula, metab.compartment)
            
            fp1.close()
            fp2.close()
            
            options.metabolicmodel = model_dirname+os.sep+'antiSMASH_model_with_template_%s.xml' %(options.modeling)
            logging.debug("Model file %s generated" % options.metabolicmodel)
        else:
            logging.warning("options.extrarecord[%s].extradata does not have a key 'MetabolicModelObjData'\navailale keys: %s" % \
                             (seq_record.id,", ".join(options.extrarecord[seq_record.id].extradata.keys())))
        
