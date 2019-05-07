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
"""
Identify conserved active site residues in PFAM_Doman / aSDomain features
"""

import logging
import sys
import re
from os import path
from helperlibs.wrappers.io import TemporaryFile
# from tempfile import NamedTemporaryFile
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

# Ignore Biopython experimental warning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO

from antismash import utils

name = "ActiveSiteFinder"
short_description = "ActiveSiteFinder identifies conserved active sites in PFAM_Domain/aSDomain features"
priority = 20000

class active_site_finder(object):
    """Active site finder class; perfoms analysis of active sites
    identified from reference aa positions of HMM profiles"""

    def __init__(self, seq_record, options):
        "Initialize ASF object"

        # Set options
        if 'activeSiteFinderConf' not in options:
            options.activeSiteFinderConf = path.join(utils.get_full_path(__file__, ''), "config", "SignatureResources.xml")

        if 'activeSiteFinderHMMDir' not in options:
            options.activeSiteFinderHMMDir = path.join(utils.get_full_path(__file__, ''), "hmm")

        # Assign variables
        try:
            XMLtree = ET.parse(options.activeSiteFinderConf)
        except ET.ParseError:
            logging.exception("Could not load/parse ActiveSiteFinder configuration file %s.", options.activeSiteFinderConf)
            sys.exit(1)

        XMLroot = XMLtree.getroot()

        HmmProfilesFilenameObj = XMLroot.findall(".//Execute/database")

        self.seq_record = seq_record
        self.options = options
        self.XMLtree = XMLtree
        self.XMLroot = XMLroot
        self.HmmProfilesFilenameObj = HmmProfilesFilenameObj

        # Check prerequisites
        filecheck = self.check_prereqs()
        if len(filecheck) > 0:
            logging.exception("required files not found:\n%s", "\n".join(filecheck))
            sys.exit(1)

    def execute(self):
        """Execute complete analysis; includes XML parsing, export of feature aa sequences,
        external tool calling, evaluation of the results and addition of new annotation to the SeqFeature"""

        XMLroot = self.XMLroot
        seq_record = self.seq_record

        # get dictionary of all features matching primary_tag_type, tag and tag value obtained from XML file

        analysisResources = XMLroot.findall('./analysis')


        for analysisResource in analysisResources:

            # set variables
            analysisResourceName = analysisResource.attrib['name']
            analysisResourceType = analysisResource.attrib['type']
            logging.debug("**********************************************\nresource name %s\n*********************************", analysisResourceName)

            # Get information on features to search for from XML file
            primaryTagType = analysisResource.find('./Prerequisite/primary_tag_type').text
            tag = analysisResource.find('./Prerequisite/tag').text
            tagValue = analysisResource.find('./Prerequisite/tag_value').text

            # get scaffold / choice subtrees form XML file
            scaffoldXML = analysisResource.find('./Alignment/scaffold')
            predictionChoicesXML = analysisResource.findall('./Alignment/choice')

            # Get pre-annoated feature where we have rules defined in XML file
            SeqFeatureList = utils.get_all_features_of_type_with_query(seq_record, primaryTagType, tag, tagValue)

            # get dictionary with domain_id as key
            SeqFeature_byID = {}
            for SeqFeature in SeqFeatureList:
                # get domain ID from asDomain_id qualifier; each feature can only have one ID, thus we don't need to loop over the list
                domain_id = SeqFeature.qualifiers['asDomain_id'][0]

                SeqFeature_byID[domain_id] = SeqFeature

            # write multi-fasta tempfile with domain aa sequences and execute tool as defined in XML file
            # and parse with Biopython SearchIO
            results = self._run_external_tool(analysisResource, SeqFeatureList)
            logging.debug("found %s hsps in hmmer results", len(results))

            # Now cycle through all results for this analysis and associate the results
            # with the corresponding SeqFeature
            for result in results:
                logging.debug("found hit with %s", result.id)

                SeqFeature = SeqFeature_byID[result.id]

                if len(result.hsps) > 0:

                    if result.hsps[0].aln[0].id != result.id:
                        # This actually should never happen...
                        logging.exception("Result ID: %s", result.id)
                        logging.exception("Was looking for hit %s but got hit for %s instead", domain_id, result.id)
                        break

                    has_annotation = False
                    # identify scaffolds and annotate
                    ASF_string = self._get_scaffold_annotation(result, scaffoldXML)

                    if ASF_string:
                        logging.debug(ASF_string)
                        if self.options.QualifierTags.asf_scaffold not in SeqFeature.qualifiers:
                            SeqFeature.qualifiers[self.options.QualifierTags.asf_scaffold] = []
                        SeqFeature.qualifiers[self.options.QualifierTags.asf_scaffold].append(ASF_string)
                        has_annotation = True

                    # identify predictions / active sites and annotate
                    (ASF_stringList, choiceResultList) = self._get_prediction_annotation(result, predictionChoicesXML)

                    for ASF_string in ASF_stringList:
                        if self.options.QualifierTags.asf_choice not in SeqFeature.qualifiers:
                            SeqFeature.qualifiers[self.options.QualifierTags.asf_choice] = []
                        logging.debug("adding ASF choice info to %s %s..%s:", SeqFeature.type, SeqFeature.location.start, SeqFeature.location.end)
                        logging.debug(ASF_string)
                        SeqFeature.qualifiers[self.options.QualifierTags.asf_choice].append(ASF_string)
                        has_annotation = True

                    for choiceResult in choiceResultList:
                        if self.options.QualifierTags.asf_prediction not in SeqFeature.qualifiers:
                            SeqFeature.qualifiers[self.options.QualifierTags.asf_prediction] = []
                        logging.debug("adding ASF choiceResult info to %s %s..%s:", SeqFeature.type, SeqFeature.location.start, SeqFeature.location.end)
                        logging.debug(choiceResult)
                        SeqFeature.qualifiers[self.options.QualifierTags.asf_prediction].append(choiceResult)

                        # Also annotate choiceResult as sec_met qualifier to corresponding CDS feature

                        correspondingCDSFeatures = utils.get_all_features_of_type_with_query(seq_record, 'CDS', 'locus_tag', utils.get_gene_id(SeqFeature))
                        if len(correspondingCDSFeatures) == 0:
                            correspondingCDSFeatures = utils.get_all_features_of_type_with_query(seq_record, 'CDS', 'gene', utils.get_gene_id(SeqFeature))
                        if not len(correspondingCDSFeatures) == 1:
                            logging.warning("ASF: found %s entries for CDS with locus tag %s, skipping", len(correspondingCDSFeatures),
                                utils.get_gene_id(SeqFeature))
                            continue
                        correspondingCDSFeature = correspondingCDSFeatures[0]
                        if 'sec_met' not in correspondingCDSFeature.qualifiers:
                            correspondingCDSFeature.qualifiers['sec_met'] = []
                        sec_met_string = "ASF-prediction: "

                        # Calculate relative locations of domains

                        if SeqFeature.strand == 1:
                            start = ((SeqFeature.location.start - correspondingCDSFeature.location.start + 3) / 3) - 1
                            end = ((SeqFeature.location.end - correspondingCDSFeature.location.start + 3) / 3) - 1
                        else:
                            start = ((correspondingCDSFeature.location.end - SeqFeature.location.end + 3) / 3) - 1
                            end = ((correspondingCDSFeature.location.end - SeqFeature.location.start + 3) / 3) - 1
                        sec_met_string += SeqFeature.qualifiers['domain'][0] + " (" + str(start) + ".." + str(end) + "): "
                        sec_met_string += choiceResult
                        logging.debug("adding ASF-prediction data to sec_met qualifier of %s", utils.get_gene_id(correspondingCDSFeature))
                        correspondingCDSFeature.qualifiers['sec_met'].append(sec_met_string)
                    if has_annotation:
                        if self.options.QualifierTags.asf_note not in SeqFeature.qualifiers:
                            SeqFeature.qualifiers[self.options.QualifierTags.asf_note] = []
                        SeqFeature.qualifiers[self.options.QualifierTags.asf_note].append("ASF analyisis with definition %s (type %s)" % \
                                (analysisResourceName, analysisResourceType))
        return True

    def check_prereqs(self):
        "Check if all required files and applications are around"

        # Tuple is ( binary_name, optional)
        _required_binaries = [
            ('blastp', False),
            ('hmmpfam2', False),
            ('hmmscan', False)
        ]

        options = self.options
        failure_messages = []

        for binary_name, optional in _required_binaries:
            if utils.locate_executable(binary_name) is None and not optional:
                failure_messages.append("Failed to locate file: %r" % binary_name)

        # Get all HMM profile names from XML file

        for HMMProfile in self.HmmProfilesFilenameObj:
            if utils.locate_file(path.join(options.activeSiteFinderHMMDir, HMMProfile.text)) is None:
                failure_messages.append("Failed to locate file: %s" % HMMProfile.text)

        return failure_messages

    def _get_scaffold_annotation(self, result, scaffoldXML):
        "generate annotation from scaffold information"

        query_seq = result.hsps[0].aln[0].seq
        hmm_seq = result.hsps[0].aln[1].seq
        scaffoldPos = scaffoldXML.find('./scaffoldOffset').text
        scaffoldValue = scaffoldXML.find('./scaffoldValue').text
        if scaffoldXML.find('./scaffoldEmission'):
            scaffoldEmissionList = scaffoldXML.find('./scaffoldEmission').text.split(',')
        else:
            scaffoldEmissionList = []
        scaffoldPosList = scaffoldPos.split(',')
        scaffoldValueList = scaffoldValue.split(',')

        overallMatch = True

        # Calculate and print match overview line and print alignment and match line
        matchLineStr = []
        for i in range(0, len(hmm_seq)):
            offset = i + result.hsps[0].hit_start + 1
            if hmm_seq[i] == ".":
                matchLineStr.append(" ")
            if str(offset) in scaffoldPosList:
                matchLineStr.append("*")
            else:
                matchLineStr.append(" ")
        logging.debug("%s %s..%s", query_seq, result.hsps[0].query_start, result.hsps[0].query_end)
        logging.debug("%s %s..%s", hmm_seq, result.hsps[0].hit_start, result.hsps[0].hit_end)
        logging.debug("".join(matchLineStr))

        # Check scaffold matches

        extracted_aa_List = []
        match_List = []
        emission_List = []

        skip = False
        for i in range(0, len(scaffoldPosList)):

            scafPos = int(scaffoldPosList[i]) - 1
            scafValue = scaffoldValueList[i]

            # Check whether scafPos is within HSP coordinates
            if (result.hsps[0].hit_start > scafPos) or (result.hsps[0].hit_end < scafPos):
                logging.warn("scaffold coordinate %s outside hsp!", scafPos)
                overallMatch = False
                skip = True
                break
            # fix position for gap characters in hmm_hit
            try:
                fixedScafPos = self._fix_coordinates(scafPos - result.hsps[0].hit_start, hmm_seq)
            except ValueError:
                logging.error("gap-fixed scaffold coordinate %s outside hsp for original position: ", scafPos)
                overallMatch = False
                skip = True
                break

            # Check whether amino acid in hmm_seq[fixedScafPos] equals predifined aa from XML
            # (should always match) and thus may be removed when the module is thoroughly tested
            # This statement fails if the fixedScafPos is outside alignment due to gaps..., so I catch the exception...
            try:
                if hmm_seq[fixedScafPos].lower() != scafValue.lower():
                    logging.warn("ASF: aa extracted from hmm profile does not match predifined aa in XML config file!")
            except IndexError:
                logging.warn("gap-fixed scaffold coordinate %s outside hsp!", fixedScafPos)
                overallMatch = False
                skip = True
                break

            extracted_aa = query_seq[fixedScafPos]
            extracted_aa_List.append(extracted_aa)
            if len(scaffoldEmissionList) > 0:
                extracted_aa_Emission = int(scaffoldEmissionList[i])
            else:
                extracted_aa_Emission = "n.d."
            emission_List.append(extracted_aa_Emission)
            match = False

            # We have to use a RegEx here to allow negations and more complex queries; ignore case (?i)
            if re.match("(?i)" + scafValue, query_seq[fixedScafPos]):
                match = True
            else:
                overallMatch = False
            match_List.append(str(match))
            logging.debug("Scaffold coordinate %s; fixed scaffold coordinate %s, query aa %s; hmm aa %s; " \
                          "scaffold value %s; emission probability %s; match %s",
                          scafPos, fixedScafPos, extracted_aa, hmm_seq[fixedScafPos], scafValue, extracted_aa_Emission, match)

        logging.debug("Overall Scaffold Match: %s\n", str(overallMatch).upper())

        # Generate Feature qualifiers

        if not skip:
            ASF_string = "Scaffold coordinates: (%s); scaffold residues: (%s); expected: (%s); matchArray: (%s); " \
                         "emission probability array (%s); overall match: %s" % \
                         (",".join(scaffoldPosList), ",".join(extracted_aa_List), ",".join(scaffoldValueList), \
                         ",".join(match_List), ",".join(emission_List), str(overallMatch).upper())

            # logging.debug("adding ASF info to %s %s..%s:" % (SeqFeature.type, SeqFeature.location.start,
            # SeqFeature.location.end))
            return ASF_string
        return

    def _get_prediction_annotation(self, result, predictionChoicesXML):
        "gererate annotation from choices/prediciton information"

        choiceResultList = []
        ASF_stringList = []
        query_seq = result.hsps[0].aln[0].seq
        hmm_seq = result.hsps[0].aln[1].seq

        for choice in predictionChoicesXML:
            extracted_aa_List = []
            emission_List = []
            match_List = []
            skip = False
            predictionOffsetList = choice.find('./offset').text.split(',')
            predictionValueList = choice.find('./value').text.split(',')
            predictionResult = choice.attrib['result']
            predicitionComment = choice.find('./comment').text
            if choice.find('./valueEmission'):
                predictionValueEmissionList = choice.find('./valueEmission').text.split(',')
            else:
                predictionValueEmissionList = []

            choiceOverallMatch = True

            logging.debug("testing %s (%s):", predictionResult, predicitionComment)

            for i in range(0, len(predictionOffsetList)):
                choiceMatch = False
                predictionOffset = int(predictionOffsetList[i]) - 1
                predictionValue = predictionValueList[i]
                if len(predictionValueEmissionList) > 0:
                    predictionValueEmission = int(predictionValueEmissionList[i])
                else:
                    predictionValueEmission = "n.d."
                emission_List.append(predictionValueEmission)

                # Check whether coordinates are within HSP match of the HMM profile
                if (result.hsps[0].hit_start > predictionOffset) or (result.hsps[0].hit_end < predictionOffset):
                    logging.warn("choice/prediction coordinate %s outside hsp!", predictionOffset)
                    choiceOverallMatch = False
                    skip = True
                    break

                try:
                    fixed_predictionOffset = self._fix_coordinates(predictionOffset - result.hsps[0].hit_start, hmm_seq)
                except ValueError:
                    logging.warn("gap-fixed choice/prediction coordinate %s outside hsp! for result: %s, which has length: %s", \
                                 predictionOffset, result.id, len(query_seq))
                    choiceOverallMatch = False
                    skip = True
                    break

                # Check whether gap-fixed coordinate still is within returned query sequence
                if len(hmm_seq) < fixed_predictionOffset + 1:
                    logging.warn("gap-fixed choice/prediction coordinate %s outside hsp! for result: %s, which has length: %s", \
                                 fixed_predictionOffset, result.id, len(query_seq))
                    choiceOverallMatch = False
                    skip = True
                    break

                try:
                    extracted_aa = query_seq[fixed_predictionOffset]
                except IndexError:
                    # This error now actually should never be raised
                    logging.error("gap-fixed choice/prediction coordinate %s outside hsp! for result: %s, which has length: %s; this should never happen", \
                                 fixed_predictionOffset, result.id, len(query_seq))
                    choiceOverallMatch = False
                    skip = True
                    break

                extracted_aa_List.append(extracted_aa)

                if re.match("(?i)" + predictionValue, extracted_aa):
                    choiceMatch = True
                else:
                    choiceOverallMatch = False

                match_List.append(str(choiceMatch))
                logging.debug("Offset %s; fixed offset %s  Expected: %s; observed in query %s: observed in hmm %s; Emission %s; match %s", \
                              predictionOffset,
                              fixed_predictionOffset,
                              predictionValue,
                              extracted_aa,
                              hmm_seq[fixed_predictionOffset],
                              predictionValueEmission,
                              choiceMatch)

            logging.debug("Overall Match for prediction %s: %s", choice.attrib['result'], str(choiceOverallMatch).upper())
            logging.debug("================================")

            ASF_string = ""
            if not skip:
                ASF_string = "Description: %s, choice result: %s, choice coordinates: (%s); residues: (%s); " \
                             "expected for choice: (%s); matchArray: (%s); emission probability array (%s); overall match: %s" % \
                              (predicitionComment,
                               predictionResult,
                               ",".join(predictionOffsetList),
                               ",".join(extracted_aa_List),
                               ",".join(predictionValueList),
                               ",".join(match_List),
                               ",".join(emission_List),
                               str(choiceOverallMatch).upper())

                ASF_stringList.append(ASF_string)

            choice_string = ""
            if choiceOverallMatch:
                choice_string = "Full match for prediction: %s" % predictionResult
                choiceResultList.append(choice_string)
        return (ASF_stringList, choiceResultList)


    def _run_external_tool(self, analysisResource, SeqFeatureList):
        "Generate tempfile containing the extracted Feature sequences and run tool defined in XML file"

        # write fasta file for the features to tempfile
        # FastA header is aSDomain_id

        tempfile = ""
        fastafile = []
        for SeqFeature in SeqFeatureList:
            fastaHeader = SeqFeature.qualifiers['asDomain_id'][0]
            fastaSeq = SeqFeature.qualifiers['translation'][0]
            # Never write empty fasta entries
            if len(fastaSeq) == 0:
                logging.warn("No translation for %s, skipping", fastaHeader)
                continue
            fastafile.append(">%s\n" % fastaHeader)
            fastafile.append("%s\n" % fastaSeq)
        querydata = "".join(fastafile)
        # DEBUG

        # print "\nTempfile for %s:\n" % analysisResource.attrib['name']
        # print "".join(fastafile)

        UseSTDIN = "False"
        executeObj = analysisResource.find('./Execute')
        if 'UseSTDIN' in executeObj.attrib:
            UseSTDIN = executeObj.attrib['UseSTDIN']
        if UseSTDIN == "False":
            with TemporaryFile(prefix='antiSMASH_ASP') as tempfile:
                out_file = open(tempfile.name, "w")
                out_file.write(querydata)
                out_file.close()

        if len(fastafile) == 0:
            logging.warn("ASP: No features found containing feature/tag/value %s / %s / %s",
                         analysisResource.find('./Prerequisite/primary_tag_type').text,
                         analysisResource.find('./Prerequisite/tag').text,
                         analysisResource.find('./Prerequisite/tag_value').text)
            return []

        results = []
        if UseSTDIN == "True":
            results = self._execute_tool(analysisResource, stdin_data=querydata)
        else:
            results = self._execute_tool(analysisResource, fileName=tempfile.name)

        return results

    def _execute_tool(self, analysisResource, fileName=None, stdin_data=None):
        "Perform the external program execution"

        cmdlineList = []

        # Assemble commad line list

        # extract program name from XML
        executeObj = analysisResource.find('./Execute')
        cmdlineList.append(executeObj.attrib['program'])

        # Cycle through parameters in XML
        for parameter in list(analysisResource.findall('./Execute/parameters/parameter')):

            if 'prefix' in parameter.attrib:
                cmdlineList.append(parameter.attrib['prefix'])
            cmdlineList.append(parameter.text)

        # Get database name
        database = analysisResource.find('./Execute/database')
        if 'prefix' in database.attrib:
            cmdlineList.append(database.attrib['prefix'])
        # Add searchpath
        cmdlineList.append(utils.locate_file(path.join(self.options.activeSiteFinderHMMDir, database.text)))

        if fileName:
            # Get (optional) input file prefix (e.g. -query in blast)
            if 'inputfile_prefix' in executeObj.attrib:
                cmdlineList.append(executeObj.attrib['inputfile_prefix'])
            cmdlineList.append(fileName)

        if stdin_data:
            # Get (optional) prefix for stdin (e.g. "-" for hmmpfam / hmmscan
            if 'STDINprefix' in executeObj.attrib:
                cmdlineList.append(executeObj.attrib['STDINprefix'])

        logging.debug("ASF: %s; external program call:\n%s", analysisResource.attrib['name'], " ".join(cmdlineList))

        try:
            if fileName:
                logging.debug("Executing tool with file input")
                out, _, retcode = utils.execute(cmdlineList)
            else:
                logging.debug("Executing tools with STDIN input")
                out, _, retcode = utils.execute(cmdlineList, input=stdin_data)
        except OSError:
            logging.warn('OS error on execution of: %s', " ".join(cmdlineList))
            return []
        if retcode != 0:
            logging.warn('%s returned %s', cmdlineList[0], retcode)
            return []
        res_stream = StringIO(out)
        logging.debug('External program output: %s', res_stream)

        # Get Biopython parser information from XML
        biopython_parser = analysisResource.find('./Execute/BioPythonParser')
        try:
            results = list(SearchIO.parse(res_stream, biopython_parser.text))
        except Exception as e:
            logging.warn('Error parsing results for active site finder analysis: %s ; no hits will be reported', e)
            results = []

        return results

    def _fix_coordinates(self, coordinate, seq):
        "Fix coordinates in hmm reference seq if gaps are present"

        # Check if coordinate is within string range
        if len(seq) - seq.count('.') - 1 < coordinate:
            logging.error('tried to fix coordinate failed; coordinate exceeds sequence length. This should never happen!')
            raise ValueError('Coordinate not within sequence')


        numberOfGaps = seq[:coordinate].count('.')

        while (coordinate < len(seq)) and (seq[coordinate] == "."):
            coordinate += 1
            logging.debug("increase coordinate by 1")
        temp_coordinate = coordinate
        new_coordinate = coordinate + numberOfGaps

        # if seq[coordinate:new_coordinate] also contains gaps, we also have to add these;
        while seq[temp_coordinate:new_coordinate].count('.') > 0:
            t = new_coordinate
            new_coordinate += seq[temp_coordinate:new_coordinate].count('.')
            temp_coordinate = t

        return new_coordinate
