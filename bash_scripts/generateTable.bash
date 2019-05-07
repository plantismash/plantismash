#!/bin/bash

HMMER_DOMTBLOUT=$1
TXTTBL="$HMMER_DOMTBLOUT.tab.gz"

if [ ! -f $HMMER_DOMTBLOUT ]; then
	echo "File $HMMER_DOMTBLOUT not found!"
	exit 1
fi

HMMER_RESULT_TARGET="$HMMER_DOMTBLOUT.target"
HMMER_RESULT_ACCESSION="$HMMER_DOMTBLOUT.accession"
HMMER_RESULT_REST="$HMMER_DOMTBLOUT.rest"
HMMER_RESULT_DESCRIPTION="$HMMER_DOMTBLOUT.desc"

cleanup() {
	    [[ -f "$HMMER_RESULT_TARGET" ]] && rm -f "$HMMER_RESULT_TARGET"
	    [[ -f "$HMMER_RESULT_ACCESSION" ]] && rm -f "$HMMER_RESULT_ACCESSION"
	    [[ -f "$HMMER_RESULT_DESCRIPTION" ]] && rm -f "$HMMER_RESULT_DESCRIPTION"
	    [[ -f "$HMMER_RESULT_REST" ]] && rm -f "$HMMER_RESULT_REST"
	    [[ -f "$HMMER_RESULT_ACCESSION" ]] && rm -f "$HMMER_RESULT_ACCESSION"
}

trap cleanup EXIT

grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f1 > $HMMER_RESULT_TARGET
grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f1 | awk -vFS='|' '{if (!$4) {print $5} else {print $4}}' > $HMMER_RESULT_ACCESSION
grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f3-22 > $HMMER_RESULT_REST
grep -v '^#' $HMMER_DOMTBLOUT | sed -e 's/\s\+/\t/g' | cut -f23- | sed -e 's/\t/ /g' | tr -d "'" > $HMMER_RESULT_DESCRIPTION

paste $HMMER_RESULT_TARGET $HMMER_RESULT_ACCESSION $HMMER_RESULT_REST $HMMER_RESULT_DESCRIPTION | gzip > $TXTTBL
