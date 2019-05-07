#!/bin/bash

GENERATEscript="/Users/tilmann/workspace/antiSMASH2/bash_scripts/generate_gb-acc_protein_id-mapping.pl"

if [ -f gi_out.txt ]; then
	rm gi_out.txt
fi

if [ -f protid_out.txt ]; then
	rm protid_out.txt
fi

if [ -f errors.txt ]; then
	rm errors.txt
fi

for fn in *.gbff; do
	cat $fn | $GENERATEscript
	if [ "$?" -ne "0" ]; then
		echo Fatal error porcessing $fn, quitting now!
		exit 1
	fi
done

