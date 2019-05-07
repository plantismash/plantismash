#!/bin/bash
set -x 
set -o pipefail

HMSEARCH_JNAME=$1
HMSEARCH_OUTPUT_DIR=$2
mail=$3

HMSEARCH_FAILED=/megx/data/megx-portal/ufBGC/hmm-search/failed-jobs/
HMSEARCH_FINISHED=/megx/data/megx-portal/ufBGC/hmm-search/finished-jobs/

ERRORS=$(ls $HMSEARCH_OUTPUT_DIR | grep -c error)

if [ "$ERRORS" -ne "0" ]; then
    echo "ERROR FINISHING JOB $HMSEARCH_JNAME"
    mv $HMSEARCH_OUTPUT_DIR $HMSEARCH_FAILED
    exit 2
else
    echo "JOB $HMSEARCH_JNAME DONE"
    mv $HMSEARCH_OUTPUT_DIR $HMSEARCH_FINISHED
        mail -s "HMM_search:$JOB_ID done." "$mail" <<EOF
JOB $JOB_ID SUCCESSFULLY DONE.
EOF
    exit 0
fi 
