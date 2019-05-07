#!/bin/bash


for fn in *.gz; do
	echo started loading $fn to database at $(date)
	DBRESULT=$({ gunzip -c $fn | psql -h localhost -d biosql -U biosql -c "\COPY asmash.pfam_table FROM STDIN DELIMITER E'\t'"; } 2>&1 )
	if [ "$?" -ne "0" ]; then
    	echo "ERROR INSERTING $fn SEARCH INTO DATABASE" > ${fn/.gz/.dberr.txt}
    	echo $DBRESULT >> ${fn/.gz/.dberr.txt}
    fi
done