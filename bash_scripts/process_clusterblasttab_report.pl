#!/opt/local/bin/perl
#
# This script processes the blastp (-outfmt 6) outputfiles against the clusterblast.fasta database and converts them in a tabular format for storing in PostgreSQL
#
# Usage zcat blastpoutputfilename | process_clusterblasttab_report.pl |psql -h localhost -d biosql -U biosql -c "\COPY asmash.clusterblast_table FROM STDIN DELIMITER E'\t'"
# Output: STDOUT



use strict;

while (my $readln = <>) {
	chomp($readln);
	my ($query_string,
		$subject_string,
		$perc_id,
		$aln_length,
		$no_mism,
		$no_gaps,
		$query_start,
		$query_end,
		$subject_start,
		$subject_end,
		$evalue,
		$bitscore ) = split(/\t/, $readln);
		
	my (undef, $gi, undef, $prot_id) = split(/\|/, $query_string);
	
	print join("\t",($gi, $prot_id, $subject_string, $perc_id, $aln_length, $no_mism, $no_gaps,
                  $query_start, $query_end, $subject_start, $subject_end, $evalue,$bitscore)), "\n";
}