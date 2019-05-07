#!/bin/bash
#
# This script extracts genebank accession number.version and the protein_ids and stores them in the asmash.gbacc_protin_mapping table
# NOTE: password to database hast to be provided via .pgpass file or PGUSER/PGPASSWORD environment variable
#
#
INPUT_FILE=$1
TMP_protid_FILE=$INPUT_FILE.protid.tmp
TMP_gid_FILE=$INPUT_FILE.gid.tmp

echo grepping GB_ACCESSION
GB_ACCESSION=$(grep "VERSION" $INPUT_FILE | perl -p -e 's/VERSION\s+(\S+).*/$1/')

echo grepping protein_id
PROTEIN_ID_LIST=$(grep "protein_id" $INPUT_FILE | perl -p -e 's/.*protein_id=\"(\S+)\".*/$1/')
echo grepping gi
GI_LIST=$(grep "db_xref=\"GI:" $INPUT_FILE | perl -p -e 's/.*db_xref=\"GI:(\S+)\".*/$1/')

#echo $GB_ACCESSION

if [ -e $TMP_protid_FILE ]; then
	rm $TMP_protid_FILE
fi

if [ -e $TMP_gid_FILE ]; then
	rm $TMP_gid_FILE
fi

echo writing protein_id_tempfile

for PROTEIN_ID in $PROTEIN_ID_LIST; do
	printf "%s\t%s\n" $GB_ACCESSION $PROTEIN_ID >> $TMP_protid_FILE
done

echo writing gi_tempfile
for GI in $GI_LIST; do
	printf "%s\t%s\n" $GB_ACCESSION $GI >> $TMP_gid_FILE
done

echo uploading content to asmash.gbacc_protid_mapping...
cat $TMP_protid_FILE | psql -h localhost -d biosql -c "\COPY asmash.gbacc_protid_mapping FROM STDIN DELIMITER E'\t'"

echo uploading content to asmash.gbacc_gi_mapping...
cat $TMP_gid_FILE | psql -h localhost -d biosql -c "\COPY asmash.gbacc_gi_mapping FROM STDIN DELIMITER E'\t'"