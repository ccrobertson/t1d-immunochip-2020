#!/bin/perl
use strict;
use warnings;

########################################################################### 
#
# Usage: perl decide_dup_family_to_drop.pl {mkin} {duplicate_families}
#
##########################################################################


my $pedfile = $ARGV[0];
my $dupfams =  $ARGV[1];


#Count number of subjects in each family
my %countFamilyMembers;
open(PED, "<$pedfile") or die("Count not open $pedfile\n");
while(<PED>) {
	my ($fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
	$countFamilyMembers{$fid}{'MEMBERS'}{$iid}=1;
	$countFamilyMembers{$fid}{'MEMBERS'}{$iid}=1;
}
close(PED);
foreach my $fid (keys %countFamilyMembers) {
	my @members = keys %{$countFamilyMembers{$fid}{'MEMBERS'}}; 
	$countFamilyMembers{$fid}{'COUNT'} = scalar @members;
	#print join("\t", $fid, $countFamilyMembers{$fid}{'COUNT'}),"\n";
}

#For each pair of duplicate families (at least one subject is in both families), choose family with most relateds
open(DUP,"<$dupfams") or die("Cannot open $dupfams\n");
while(<DUP>) {
	chomp;
	my ($fid1, $fid2) = split(/\s+/, $_);
	#print join("\t", $fid1, $countFamilyMembers{$fid1}{'COUNT'}, $fid2, $countFamilyMembers{$fid2}{'COUNT'}),"\n";
	if ($countFamilyMembers{$fid1}{'COUNT'} lt $countFamilyMembers{$fid2}{'COUNT'}) {
		foreach my $id (keys %{$countFamilyMembers{$fid1}{'MEMBERS'}}) {
			print join("\t", $fid1, $id),"\n";
		}
	} else {
		foreach my $id (keys %{$countFamilyMembers{$fid2}{'MEMBERS'}}) {
			print join("\t", $fid2, $id),"\n";
		}
	} 
}
close(DUP);