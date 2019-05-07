#!/bin/bash
# Don't forget to set PGPASSWORD env variable to allow database access without entering password!

PROCESS_CLUSTERBLAST='/Users/tilmann/workspace/antiSMASH2/bash_scripts/process_clusterblasttab_report.pl'

if [ ! -f $PROCESS_CLUSTERBLAST ]; then
	echo "Script $PROCESS_CLUSTERBLAST not found!"
	exit 1
fi

for fn in *.subclusterprot.blastp.out.gz; do
	if [ -f ${fn/.gz/.dberr.txt} ]; then
		rm ${fn/.gz/.dberr.txt}
	fi
	echo started loading $fn to database at $(date)
	DBRESULT=$({ gunzip -c $fn | $PROCESS_CLUSTERBLAST | psql -h localhost -d biosql -U biosql -c "\COPY asmash.subclusterblast_table FROM STDIN DELIMITER E'\t'"; } 2>&1 )
	if [ "$?" -ne "0" ]; then
    	echo "ERROR INSERTING $fn SEARCH INTO DATABASE" > ${fn/.gz/.dberr.txt}
    	echo $DBRESULT >> ${fn/.gz/.dberr.txt}
    fi
done