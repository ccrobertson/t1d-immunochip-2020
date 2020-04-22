#!/bin/perl
use strict;
use warnings;


my $phenofile = $ARGV[0];
my $id_update = $ARGV[1];


my %phenoIDs;
open(OLD, "<$phenofile") or die("Cannot open $phenofile\n");
while(<OLD>) {
	chomp;
	if ($_=~/^Cohort/) {
		next;
	}
	my ($cohort, $fid, $iid, $fat, $mot, $sex, $t1d) = split(/\s+/, $_);
	$phenoIDs{"$fid:$iid"}=1;
}
close(OLD);


open(UPDATE, "<$id_update") or die("Cannot open $id_update\n");
my %oldToNew;
while(<UPDATE>) {
	chomp;
	my ($fid_old, $iid_old, $fid_new, $iid_new) = split(/\s+/, $_);
	if (exists($phenoIDs{"$fid_new:$iid_new"})) {
		next;
	} else {
		$oldToNew{"$fid_old:$iid_old"}{'FID'} = $fid_new;
		$oldToNew{"$fid_old:$iid_old"}{'IID'} = $iid_new;
	}
}
close(UPDATE);


open(OLD, "<$phenofile") or die("Cannot open $phenofile\n");
while(<OLD>) {
	chomp;
	if ($_=~/^Cohort/) {
		print $_, "\n";
		next;
	}
	
	my ($cohort, $fid, $iid, $fat, $mot, $sex, $t1d) = split(/\s+/, $_);	
		
	my $fid_new;
	my $iid_new;
	if (exists($oldToNew{"$fid:$iid"})) {
		$fid_new = $oldToNew{"$fid:$iid"}{'FID'};
		$iid_new = $oldToNew{"$fid:$iid"}{'IID'}; 
	} else {
		$fid_new = $fid;
		$iid_new = $iid;
	}
	print join("\t", $cohort, $fid_new, $iid_new, $fat, $mot, $sex, $t1d),"\n";
	
}
close(OLD);