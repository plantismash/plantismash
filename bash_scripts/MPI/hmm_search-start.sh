#!/bin/bash
set -x 
set -o pipefail

NUM_THREADS=2

HMSEARCH_PATH=/megx/software/ufBGC/hmm_search/
HMM_PROFILE_DIR=$HMSEARCH_PATH"hmm_profiles/"

# Update profiles
cd $HMM_PROFILE_DIR
svn update

NUM_PROFILES=$(ls $HMM_PROFILE_DIR | wc -l)

# Prepare output folder
HMSEARCH_RUNNING=/megx/data/megx-portal/ufBGC/hmm-search/running-jobs/
HMSEARCH_JNAME=$(date -d "today" +"hmm_search-%Y%m%d%H%M")
HMSEARCH_OUTPUT_DIR=$HMSEARCH_RUNNING$HMSEARCH_JNAME

# Create output folder
[[ -d "$HMSEARCH_OUTPUT_DIR" ]] && rm -rf "$HMSEARCH_OUTPUT_DIR"
mkdir -p $HMSEARCH_OUTPUT_DIR

cd $HMSEARCH_OUTPUT_DIR

HMSEARCH_JNAME_FINISH=$HMSEARCH_JNAME"-finish"

mail="afernand@mpi-bremen.de"

echo "Submitting job array..."
qsub -t 1-$NUM_PROFILES -o $HMSEARCH_OUTPUT_DIR -q cluster.q -e $HMSEARCH_OUTPUT_DIR -l ga -j y -terse -R y -m s -M $mail -N $HMSEARCH_JNAME -pe threaded $NUM_THREADS $HMSEARCH_PATH"hmm_search.sh" $HMSEARCH_JNAME $HMM_PROFILE_DIR $HMSEARCH_OUTPUT_DIR $mail

echo "Finalizing job..."
qsub -hold_jid $HMSEARCH_JNAME -o $HMSEARCH_OUTPUT_DIR -q cluster.q -e $HMSEARCH_OUTPUT_DIR -l ga -j y -terse -R y -m s -M $mail -N $HMSEARCH_JNAME_FINISH $HMSEARCH_PATH"hmm_search-finish.sh" $HMSEARCH_JNAME $HMSEARCH_OUTPUT_DIR $mail
