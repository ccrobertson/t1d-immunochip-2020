#!/bin/perl
use strict;
use warnings;


my $phenofile = $ARGV[0];
my $sex_update = $ARGV[1];


open(UPDATE, "<$sex_update") or die("Cannot open $sex_update\n");
my %iidToSex;
while(<UPDATE>) {
	chomp;
	my ($fid, $iid, $sex) = split(/\s+/, $_);
	$iidToSex{"$fid:$iid"} = $sex;
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
		
	my $new_sex;
	if (exists($iidToSex{"$fid:$iid"})) {
		$new_sex = $iidToSex{"$fid:$iid"}; 
	} else {
		$new_sex = $sex;
	}
	print join("\t", $cohort, $fid, $iid, $fat, $mot, $new_sex, $t1d),"\n";
	
}
close(OLD);