#!/opt/local/bin/perl
#
# This script extracts genebank accession number.version and the protein_ids and stores them in tab separated text files
# Input is read from STDIN and file
# recommended usage zcat <filename> | generate_gb_acc_protein_id_mapping.pl --gifile <gi_outfile> --protidfile <protid_outfile>
#
# If no output files are specified, defaults gi_out.txt and protid_out.txt are used


#TODO: Test with fungal sequence containing introns


use strict;
use Getopt::Long;


my $gi_outf = "gi_out.txt";
my $protid_outf = "protid_out.txt";
my $err_outf = "errors.txt";

GetOptions("gifile|g=s"		=> \$gi_outf,
		   "protidfile|p=s"	=> \$protid_outf,
		   "errfile|e=s" => \$err_outf );

open (GIFILE, ">>$gi_outf") || die "can't open $gi_outf for writing!\n";
open (PROTIDFILE, ">>$protid_outf") || die "cant't open $protid_outf for writing!\n";
open (ERRFILE, ">>$err_outf") || die "can't open $err_outf for writing \n";

my $genbank_id;
my $gi;
my $prot_id;
my %gi_hash;
my $start;
my $end;
my $strand;
my %protid_hash;


while (my $readln = <>) {
	chomp($readln);
	if ($readln =~ /^VERSION\s+(\S+).*/ ){
		$genbank_id = $1;
		$gi = undef;
		$prot_id = undef;
		$start = undef;
		$end = undef;
		next;
	}
	# Check for CDS coordinates split over two lines; script will fail if > 2lines...
	if (($readln =~/\s{5}CDS/) && ($readln !~/[\d\)]$/)) {
		my $nextline = <>;
		$nextline =~s/^\s+//;
		$readln .= $nextline;
	}

	if ($readln =~ /\s{5}CDS.*?(\d+).*[\.,]>?(\d+)\)*$/ ) {
		$start = $1;
		$end = $2;
		$strand = 1;
	}
	if ($readln =~/\s{5}CDS\s+complement/ ) {
		$strand = -1;
	}
	
	if ($readln =~ /.*\/db_xref=\"GI:(\S+)\".*/ ) {
		$gi = $1;
		
		if ((! defined $start) || (! defined $end)) {
			die "Start/end coordinates not assigned for $genbank_id, $gi"
		}
		if (not defined $gi_hash{$genbank_id.$gi}) {
			print GIFILE join("\t", $genbank_id, $gi, $start, $end, $strand), "\n";
			$gi_hash{$genbank_id.$gi} = defined;
		}		
		else {
			print ERRFILE "Duplicate:\t$genbank_id\t$gi\n";
		}
		
	}
	if ($readln =~ /.*\/protein_id=\"(\S+)\".*/ ) {
		$prot_id = $1;
		
		if ((! defined $start) || (! defined $end)) {
			die "Start/end coordinates not assigned for $genbank_id, $prot_id"
		}
		if (not defined $protid_hash{$genbank_id.$prot_id}) {
			print PROTIDFILE join("\t", $genbank_id, $prot_id, $start, $end, $strand), "\n";
			$protid_hash{$genbank_id.$prot_id} = defined;
		}
		else {
			print ERRFILE "Duplicate:\t$genbank_id\t$prot_id\n";
		}
	}
}