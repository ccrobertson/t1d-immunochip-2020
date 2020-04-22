#!/bin/perl
use strict;
use warnings;


my $kinfile = $ARGV[0];
my $pedfile = $ARGV[1];
my $phenofile = $ARGV[2];

#note, we must use the updated fam files from each iteration to get parent info !!!

my $iidToFatid = {};
my $iidToMotid = {};
open(PED, "<$pedfile") or die("Count not open $pedfile\n");
while(<PED>) {
	my ($fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
	$iidToFatid->{$iid} = $fat;
	$iidToMotid->{$iid} = $mot;
}
close(PED);


my $iidToCohort = {};
open(PHENO, "<$phenofile") or die("Count not open $phenofile\n");
while(<PHENO>) {
	my ($cohort, $fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
	$iidToCohort->{$iid} = $cohort;
}
close(PHENO);


#my $kinfile = $ARGV[0];
#my $pedfile = $ARGV[1];
#my $phenofile = $ARGV[2];
#
#my $iidToFatid = {};
#my $iidToMotid = {};
#
#open(PED, "<$pedfile") or die("Count not open $pedfile\n");
#while(<PED>) {
#	my ($fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
#	$iidToFatid->{$iid} = $fat;
#	$iidToMotid->{$iid} = $mot;
#}
#close(PED);


open(KIN, "<$kinfile") or die("Could not open $kinfile\n");
while(<KIN>) {
	
	chomp;
	if ($_ =~ /SNPREL/) {
		print join("\t", "FID1", "IID1", "FID2", "IID2", "Phi", "PEDREL", "SNPREL", "PEDCHILD", "PEDPARENT", "COHORT1", "COHORT2"),"\n";
		next;	
	}
	
	my ($fid1, $iid1, $fid2, $iid2, $phi, $snprel) = split(/\s+/, $_);
			
	if (!(exists $iidToCohort->{$iid1} && exists $iidToCohort->{$iid2})) {
		
		print join("\t", $fid1, $iid1, $fid2, $iid2, $phi, "UN", $snprel, "NA", "NA", "None", "None"),"\n";
		next;
	
	} else {
		
		my $cohort1 = $iidToCohort->{$iid1};
		my $cohort2 = $iidToCohort->{$iid2};
			
		#PO relationships
		#iid1 is child of iid2
		if ($iidToFatid->{$iid1} eq $iid2 || $iidToMotid->{$iid1} eq $iid2) {
			print join("\t", $fid1, $iid1, $fid2, $iid2, $phi, "PO", $snprel, $iid1, $iid2, $cohort1, $cohort2),"\n";
		}
		#iid2 is child of iid1
		elsif ($iidToFatid->{$iid2} eq $iid1 || $iidToMotid->{$iid2} eq $iid1) {
			print join("\t", $fid1, $iid1, $fid2, $iid2, $phi, "PO", $snprel, $iid2, $iid1, $cohort1, $cohort2),"\n";
		}
		
		#FS relationships (same father and mother)
		elsif ($iidToFatid->{$iid1} eq $iidToFatid->{$iid2} && $iidToMotid->{$iid1} eq $iidToMotid->{$iid2} && $iidToFatid->{$iid1} ne "0" && $iidToMotid->{$iid1} ne "0") {
			print join("\t", $fid1, $iid1, $fid2, $iid2, $phi, "FS", $snprel, "NA", "NA", $cohort1, $cohort2),"\n";
		}
		
		#Half siblings (only one parent is the same) 
		elsif ( ($iidToFatid->{$iid1} eq $iidToFatid->{$iid2} && $iidToFatid->{$iid1} ne "0") || ($iidToMotid->{$iid1} eq $iidToMotid->{$iid2} && $iidToMotid->{$iid1} ne "0") ) {
			print join("\t", $fid1, $iid1, $fid2, $iid2, $phi, "HS", $snprel, "NA", "NA", $cohort1, $cohort2),"\n";
		}
		
		#Unrelated (no 1st degree relation -- ignorning 2nd and 3rd degree relationships)
		else {
			print join("\t", $fid1, $iid1, $fid2, $iid2, $phi, "UN", $snprel, "NA", "NA", $cohort1, $cohort2),"\n";
		}
	
	}
}
close(KIN);


