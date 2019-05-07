#!/bin/bash
set -x 
set -o pipefail

HMSEARCH_JNAME=$1
HMM_PROFILE_DIR=$2
HMMER_OUTPUT_DIR=$3"/"
mail=$4

HMSEARCH_FAILED=/megx/data/megx-portal/ufBGC/hmm-search/failed-jobs/

HMMER_PATH=/usr/local/software/hmmer-3.0b3-threads/bin/
HMMER_PROGRAM=hmmsearch
HMMER_EXEC=$HMMER_PATH$HMMER_PROGRAM


HMM_PROFILE=$(ls $HMM_PROFILE_DIR | tail -n+${SGE_TASK_ID} | head -n1)
HMM_PROFILE_SHORT=$(echo $HMM_PROFILE | sed -e "s/hmm$//g") 
HMMER_INPUT=$HMM_PROFILE_DIR$HMM_PROFILE
HMMER_DBFILE=/local/biodb/ufBGC/fasta/nr.fasta
HMMER_DOMTBLOUT=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"out"
HMMER_ERROR=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"error"
DB_ERROR=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"db.error"


HMMER_RESULT_TARGET=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"target"
HMMER_RESULT_ACCESSION=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"accession"
HMMER_RESULT_DESCRIPTION=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"description"
HMMER_RESULT_REST=$HMMER_OUTPUT_DIR$HMM_PROFILE_SHORT"rest"


cleanup() {
    [[ -f "$HMMER_RESULT_TARGET" ]] && rm -f "$HMMER_RESULT_TARGET"
    [[ -f "$HMMER_RESULT_ACCESSION" ]] && rm -f "$HMMER_RESULT_ACCESSION"
    [[ -f "$HMMER_RESULT_DESCRIPTION" ]] && rm -f "$HMMER_RESULT_DESCRIPTION"
    [[ -f "$HMMER_RESULT_REST" ]] && rm -f "$HMMER_RESULT_REST"
    [[ -f "$HMMER_RESULT_ACCESSION" ]] && rm -f "$HMMER_RESULT_ACCESSION"
}

trap cleanup EXIT

HMMER_RESULT=$({ $HMMER_EXEC --cpu $NSLOTS --domtblout $HMMER_DOMTBLOUT $HMMER_INPUT $HMMER_DBFILE ; } 2>&1)

if [ "$?" -ne "0" ]; then
    echo "ERROR PROCESSING PROFILE $HMM_PROFILE" > $HMMER_ERROR
    qdel $HMSEARCH_JNAME
    mail -s "HMM_search: job $JOB_ID task $SGE_TASK_ID failed" "$mail" <<EOF
ERROR PROCESSING PROFILE $HMM_PROFILE
Result: $HMMER_RESULT
EOF
    exit 2
fi 

grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f1 > $HMMER_RESULT_TARGET
grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f1 | awk -vFS='|' '{if (!$4) {print $5} else {print $4}}' > $HMMER_RESULT_ACCESSION
grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f3-22 > $HMMER_RESULT_REST
grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f23- | sed -e 's/\t/ /g' | tr -d "'" > $HMMER_RESULT_DESCRIPTION

DB_RESULT=$({ paste $HMMER_RESULT_TARGET $HMMER_RESULT_ACCESSION $HMMER_RESULT_REST $HMMER_RESULT_DESCRIPTION | psql -h antares -p 5434 -d megdb_r8 -c "\COPY ufbgc.hmmer_results FROM STDIN DELIMITER E'\t'"; } 2>&1 )

if [ "$?" -ne "0" ]; then
    echo "ERROR INSERTING $HMM_PROFILE SEARCH INTO DATABASE" > $DB_ERROR
    qdel $HMSEARCH_JNAME
    mail -s "HMM_search:job $JOB_ID task $SGE_TASK_ID failed" "$mail" <<EOF
ERROR INSERTING $HMM_PROFILE SEARCH INTO DATABASE
Result: $DB_RESULT
EOF
    exit 2
fi 


