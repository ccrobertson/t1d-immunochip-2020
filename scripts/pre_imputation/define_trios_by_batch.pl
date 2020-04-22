#!/bin/perl

my $group = $ARGV[0];
my $phenofile = $ARGV[1];
my $clusterfile = $ARGV[2];


my %groupIds;
open(CLUST,"<$clusterfile") or die("Cannot open $clusterfile\n");
while(<CLUST>) {
	chomp;
	my ($fid, $iid) = split(/\s+/, $_);

	#track iids in group
	$groupIds{$iid}=1;

}
close(CLUST);

my %iids;
my %iidToFat;
my %iidToMot;
my %iidToFid;
my %fidToCohort;
my %affectedPerson;
my %familySize;
open(IN, "<$phenofile") or die("Cannot open $phenofile");
while(<IN>) {
	chomp;
	my ($cohort, $fid, $iid, $fat, $mot, $sex, $t1d) = split(/\s+/,$_);

	if (exists($groupIds{$iid})) {
	#track iids
	$iids{$iid}=1;

	#track parents and fid
	$iidToFat{$iid}=$fat;
	$iidToMot{$iid}=$mot;
	$iidToFid{$iid}=$fid;

	#track cohort by family
	$fidToCohort{$fid}=$cohort;

	#track affected status
	if ($t1d eq "2") {
		$affectedPerson{$iid}=1;
	}

	#track how many people are in each family
	$familySize{$fid}++;

	}
}
close(IN);


#go through each affected individual
my %familiesToKeep;
my %affectedChildToKeep;
for my $iid (keys %affectedPerson) {

	#if their father and mother were also genotyped, keep the family in TDT analysis
	if (exists($iids{$iidToFat{$iid}}) && exists($iids{$iidToMot{$iid}})) {

		$familiesToKeep{$iidToFid{$iid}}=1;
		$affectedChildToKeep{$iid}=1;
	}


}

#get iids for all members of all families
for my $iid (keys %iids) {

	if (exists($familiesToKeep{$iidToFid{$iid}})) {
		print join("\t", $fidToCohort{$iidToFid{$iid}}.$iidToFid{$iid}, $iid),"\n";
	}

}
