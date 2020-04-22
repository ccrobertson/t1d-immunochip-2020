#!/bin/perl
use strict;
use warnings;


my $phenofile = $ARGV[0];
my $exportmap = $ARGV[1];
my $phenofile_with_export = $ARGV[2]; 

open(MAP, "<$exportmap") or die("Cannot open $exportmap\n");
my %iidToExport;
while(<MAP>) {
	chomp;
	my ($fid, $iid, $export) = split(/\s+/, $_);
	$iidToExport{$iid} = $export;
}
close(MAP);


open(PHENOEXP,">$phenofile_with_export") or die("Cannot open $phenofile_with_export\n");
open(PHENO, "<$phenofile") or die("Cannot open $phenofile\n");
my %fidToExport;
while(<PHENO>) {
	chomp;
	my ($cohort, $fid, $iid, $fat, $mot, $sex, $t1d) = split(/\s+/, $_);
	
	my $export;
	if (exists($iidToExport{$iid})) {
		$export = $iidToExport{$iid};
		$fidToExport{$fid}{$export}=1;
	} else {
		$export = "0"
	}
	
	print PHENOEXP join("\t", $cohort, $fid, $iid, $fat, $mot, $sex, $t1d, $export),"\n";
	
}
close(PHENO);
close(PHENOEXP);

foreach my $fid (keys %fidToExport) {
	
	my @exports = keys %{$fidToExport{$fid}};
	if(scalar @exports > 1) {
		print join("\t", $fid, @exports),"\n";
	}
	
}