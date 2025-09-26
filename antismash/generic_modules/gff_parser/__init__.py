import logging
import sys
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from BCBio import GFF


def check_gff_suitability(options, sequences):
    if options.gff3:
        options.gff_ids = []
        # Some GFFs have a header, but some GFF parser functions break with it, so check for header and error out if
        # if exists.
        try:
            with open(options.gff3) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    else:
                        int(line.split('\t')[3])        # 4th column has to be a number (start)
                        int(line.split('\t')[4])        # 5th column has to be a number (end)
        except ValueError as e:
            logging.error('Parsing %r failed: %s', options.gff3, e)
            logging.error('It appears %r has a header. It should be removed or commented out for proper parsing.',
                          options.gff3)
            sys.exit(1)
        try:
            examiner=GFF.GFFExaminer()
            gff_data=examiner.available_limits(open(options.gff3))
            # Check if at least one GFF locus appears in sequence
            gff_ids=set([n[0] for n in gff_data['gff_id']])
            options.gff_ids = list(gff_ids)
            if len(gff_ids)==1 and len(options.all_record_ids)==1:
                # If both inputs only have one record, assume is the same, but first check coordinate compatibility
                logging.info("GFF3 and sequence have only one record. Assuming is the same as long as coordinates are"
                             "compatible.")
                limit_info=dict(gff_type=['CDS'])
                for record in GFF.parse(open(options.gff3),limit_info=limit_info):
                    break
                coord_max=max([n.location.end.real for n in record.features])
                if coord_max>len(sequences[0]):
                    logging.error('GFF3 record and sequence coordinates are not compatible.')
                    sys.exit(1)
                else:
                    options.single_entries=True
            elif len(gff_ids.intersection(set(options.all_record_ids)))==0:
                logging.error('No GFF3 record IDs match any sequence record IDs.')
                sys.exit(1)
            else:
                options.single_entries=False
            # Check GFF contains CDSs
            if not ('CDS',) in gff_data['gff_type']:
                logging.error('GFF3 does not contain any CDS.')
                sys.exit(1)
            # Check CDS are childless but not parentless
            if 'CDS' in set([n for key in examiner.parent_child_map(open(options.gff3)) for n in key]):
                logging.error('GFF3 structure is not suitable. CDS features must be childless but not parentless.')
                sys.exit(1)
        except AssertionError as e:
            logging.error('Parsing %r failed: %s', options.gff3, e)
            sys.exit(1)


def run(sequence, options):
    handle=open(options.gff3)
    # If there's only one sequence in both, read all, otherwise, read only appropriate part of GFF3.
    if options.single_entries:
        limit_info=False
    else:
        limit_info=dict(gff_id=[sequence.id])

    for record in GFF.parse(handle,limit_info=limit_info):
        for feature in record.features:
            if feature.type=='CDS':
                newFeature=feature
            else:
                newFeature=check_sub(feature,sequence,options)
                if not newFeature:
                    continue
                if not type(newFeature)==list:
                    newFeature=[newFeature]
                newFeature=[_f for _f in newFeature if _f]
                for n in newFeature:
                    if 'gene' in feature.qualifiers:
                        n.qualifiers['gene']=feature.qualifiers['gene']
                    elif 'name' in feature.qualifiers:
                        n.qualifiers['gene']=feature.qualifiers['name']
                    elif 'Name' in feature.qualifiers:
                        n.qualifiers['gene']=feature.qualifiers['Name']
                    else:
                        n.qualifiers['gene']=[feature.id]
                    sequence.features.append(n)


def check_sub(feature, sequence, options):
    # whether to use phase (codon start) to modify reported locations
    # Augustus, NCBI, and glimmerhmm report phase but have already adjusted the
    # locations and since they're the bulk of inputs, disable further modification
    if options.use_phase:
        MODIFY_LOCATIONS_BY_PHASE = True
    else:
        MODIFY_LOCATIONS_BY_PHASE = False

    newFeature=[]
    LocList=[]
    TransLocList=[]
    QualList={}
    topop=[]
    for sub in feature.sub_features:
        if sub.sub_features:                                # If there are sub_features, go deeper
            deepFeature=check_sub(sub,sequence,options)
            for n in deepFeature:
                newFeature.append(n)
        elif sub.type=='CDS':
            loc=[sub.location.start.real,sub.location.end.real]
            if MODIFY_LOCATIONS_BY_PHASE:
                if 'phase' in sub.qualifiers:
                    phase=int(sub.qualifiers['phase'][0])
                    if sub.strand==1:
                        loc[0]+=phase
                    else:
                        loc[1]-=phase
            LocList.append(FeatureLocation(loc[0],loc[1],strand=sub.strand))
            # Make sure CDSs lengths are multiple of three. Otherwise extend to next full codon.
            # This only applies for translation.
            modulus=(loc[1]-loc[0])%3
            if modulus==0:
                TransLocList.append(FeatureLocation(loc[0],loc[1],strand=sub.strand))
            else:
                if sub.strand==1:
                    TransLocList.append(FeatureLocation(loc[0],loc[1]+(3-modulus),strand=sub.strand))
                else:
                    TransLocList.append(FeatureLocation(loc[0]-(3-modulus),loc[1],strand=sub.strand))
            # For split features (CDSs), the final feature will have the same qualifiers as the children ONLY if
            # they're the same, i.e.: all children have the same "protein_ID" (key and value).
            for qual in list(sub.qualifiers.keys()):
                if not qual in QualList:
                    QualList[qual]=sub.qualifiers[qual]
                if qual in QualList and not QualList[qual]==sub.qualifiers[qual]:
                    topop.append(qual)

    for n in topop:                                            # Pop mismatching qualifiers over split features
        QualList.pop(n,None)
    QualList.pop('Parent',None)                                # Pop parent.

    newFeature.append(merge_cds(LocList, TransLocList, QualList, sequence))

    return newFeature


def merge_cds(LocList, TransLocList, QualList, sequence):
    if len(LocList) > 1:
        LocList = sorted(LocList, key=lambda x: x.start.real)
        TransLocList = sorted(TransLocList, key=lambda x: x.start.real)
        if LocList[0].strand == 1:
            newLoc = CompoundLocation(LocList)
        else:
            newLoc = CompoundLocation(list(reversed(LocList)))
            TransLocList = reversed(TransLocList)
    elif len(LocList) == 0:
        return None
    else:
        newLoc = LocList[0]

    cur_feature = SeqFeature(newLoc)
    cur_feature.qualifiers = QualList
    cur_feature.type = 'CDS'
    trans = ''.join([n.extract(sequence.seq).translate(stop_symbol='')._data for n in TransLocList])
    cur_feature.qualifiers['translation'] = [trans]

    return cur_feature
