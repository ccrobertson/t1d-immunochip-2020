#!/bin/perl
use strict;
use warnings;



my $mkin =  $ARGV[0]; 
my $expectedParentToChild = {};
my $inferredPO = {};
open(KIN1, "<$mkin") or die("Cannot open $mkin\n");
while(<KIN1>) {
	
	chomp;
	
	my ($fid1, $iid1, $fid2, $iid2, $phi, $pedrel, $snprel, $child,  $parent) = split("\t", $_);
	my $uniq1 = join("\t", $fid1, $iid1);
	my $uniq2 = join("\t", $fid2, $iid2);
			
	if ($pedrel eq "PO") {
		if ($child eq $iid1) {
			$expectedParentToChild -> {$uniq2} = $uniq1;	
		} elsif ($child eq $iid2) {
			$expectedParentToChild -> {$uniq1} = $uniq2;	
		} else {
			die("ERROR");
		}
	}
	
	if ($snprel eq "PO") {
		$inferredPO -> {$uniq1} = $uniq2;
		$inferredPO -> {$uniq2} = $uniq1;
	}


} 
close(KIN1);

open(KIN2, "<$mkin") or die("Cannot open $mkin\n");
while(<KIN2>) {
	
	chomp;
	
	my ($fid1, $iid1, $fid2, $iid2, $phi, $pedrel, $snprel, $child,  $parent) = split("\t", $_);
	my $uniq1 = join("\t", $fid1, $iid1);
	my $uniq2 = join("\t", $fid2, $iid2);
	
	my $uniq_child;
	if ($child eq $iid1) {
		$uniq_child = join("\t", $fid1, $iid1);
	} 
	elsif ($child eq $iid2) {
		$uniq_child = join("\t", $fid2, $iid2);
	} else {
		$uniq_child = "NA";
	}
	
	if (exists $inferredPO->{$uniq_child}) {
		my $infParent = $inferredPO->{$uniq_child};
		
		if (exists $expectedParentToChild-> {$infParent}){
			my $candidate = $expectedParentToChild->{$infParent};
			if ($uniq_child ne $candidate) {
				print join("\t", $uniq_child, $candidate),"\n";	
			}
		}		
	}
}
close(KIN2);



