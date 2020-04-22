#!/bin/perl
use strict;
use warnings;

########################################################################### 
#
# Usage: perl decide_dup_family_to_drop.pl {mkin} {duplicate_families} {phenofile}
#
##########################################################################


my $pedfile = $ARGV[0];
my $dupfams =  $ARGV[1];
my $phenofile = $ARGV[2];

#Get cohort ids for each family
my %iidToCohort;
my %iidToFid;
my %countCohortsPerFID;
open(PHENO, "<$phenofile") or die("Count not open $phenofile\n");
while(<PHENO>) {
	my ($cohort, $fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
	$iidToCohort{$iid} = $cohort;
	$iidToFid{$iid} = $fid;
	$countCohortsPerFID{$fid}{$cohort}=1;
}
close(PHENO);



#Count number of subjects in each family
my %countFamilyMembers;
open(PED, "<$pedfile") or die("Count not open $pedfile\n");
while(<PED>) {
	my ($fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
	if (exists $iidToCohort{$iid}) {
		my $uniq_fid = $iidToCohort{$iid}.$fid;
		$countFamilyMembers{$uniq_fid}{'MEMBERS'}{$iid}=1;
		$countFamilyMembers{$uniq_fid}{'MEMBERS'}{$iid}=1;
	}
}
close(PED);
foreach my $uniq_fid (keys %countFamilyMembers) {
	my @members = keys %{$countFamilyMembers{$uniq_fid}{'MEMBERS'}}; 
	$countFamilyMembers{$uniq_fid}{'COUNT'} = scalar @members;
	#print join("\t", $uniq_fid, $countFamilyMembers{$uniq_fid}{'COUNT'}),"\n";
}

#open(OUT, ">testing.txt") or die("Cannot open testing.txt\n");
#foreach my $fid (keys %countCohortsPerFID) {
#	print OUT join("\t", $fid, scalar keys %{$countCohortsPerFID{$fid}}, keys %{$countCohortsPerFID{$fid}}),"\n";
#}
#close(OUT);


#For each pair of duplicate families (at least one subject is in both families), choose family with most relateds
open(DUP,"<$dupfams") or die("Cannot open $dupfams\n");
while(<DUP>) {
	chomp;
	my ($fid1, $fid2) = split(/\s+/, $_);
	
	my @family1 = keys %{$countFamilyMembers{$fid1}{'MEMBERS'}};
	my @family2 = keys %{$countFamilyMembers{$fid2}{'MEMBERS'}};
	
	my $cohort1 = $iidToCohort{$family1[0]};
	my $cohort2 = $iidToCohort{$family2[0]};
	
	if ($cohort1=~/^T1DGC/ && $cohort2 =~/^T1DGC/) {
		if ($countFamilyMembers{$fid1}{'COUNT'} lt $countFamilyMembers{$fid2}{'COUNT'}) {
			foreach my $id (keys %{$countFamilyMembers{$fid1}{'MEMBERS'}}) {
				print join("\t", $iidToFid{$id}, $id, $cohort1),"\n";
			}
		} else {
			foreach my $id (keys %{$countFamilyMembers{$fid2}{'MEMBERS'}}) {
				print join("\t", $iidToFid{$id}, $id, $cohort2),"\n";
			}
		}
		 		
	} elsif ($cohort1=~/^T1DGC/ && $cohort2 !~/^T1DGC/) {
		foreach my $id (keys %{$countFamilyMembers{$fid2}{'MEMBERS'}}) {
			print join("\t", $iidToFid{$id}, $id, $cohort2),"\n";
		}
		
	} elsif ($cohort1!~/^T1DGC/ && $cohort2 =~/^T1DGC/) {
		foreach my $id (keys %{$countFamilyMembers{$fid1}{'MEMBERS'}}) {
			print join("\t", $iidToFid{$id}, $id, $cohort1),"\n";
		}
		
	} elsif ($cohort1!~/^T1DGC/ && $cohort2 !~/^T1DGC/) {
		if ($countFamilyMembers{$fid1}{'COUNT'} lt $countFamilyMembers{$fid2}{'COUNT'}) {
			foreach my $id (keys %{$countFamilyMembers{$fid1}{'MEMBERS'}}) {
				print join("\t", $iidToFid{$id}, $id, $cohort1),"\n";
			}
		} else {
			foreach my $id (keys %{$countFamilyMembers{$fid2}{'MEMBERS'}}) {
				print join("\t", $iidToFid{$id}, $id, $cohort2),"\n";
			}
		} 
		
	} else {
		die("ERROR\n");
	}
	
	
}
close(DUP);