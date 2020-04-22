#!/bin/perl
use strict;
use warnings;


my $phenofile = $ARGV[0];
my $parent_update = $ARGV[1];


open(UPDATE, "<$parent_update") or die("Cannot open $parent_update\n");
my %iidToParents;
while(<UPDATE>) {
	chomp;
	my ($fid, $iid, $fat, $mot) = split(/\s+/, $_);
	$iidToParents{"$fid:$iid"}{'FAT'} = $fat;
	$iidToParents{"$fid:$iid"}{'MOT'} = $mot;
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

		
	my $new_fat;
	my $new_mot;
	if (exists($iidToParents{"$fid:$iid"})) {
		$new_fat = $iidToParents{"$fid:$iid"}{'FAT'};
		$new_mot = $iidToParents{"$fid:$iid"}{'MOT'}; 
	} else {
		$new_fat = $fat;
		$new_mot = $mot;
	}
	print join("\t", $cohort, $fid, $iid, $new_fat, $new_mot, $sex, $t1d),"\n";
	
}
close(OLD);