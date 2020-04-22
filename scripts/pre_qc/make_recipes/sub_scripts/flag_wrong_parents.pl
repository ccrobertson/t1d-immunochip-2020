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
#print join("\t", "FIDCHILD", "IIDCHILD", "IIDPARENT", "FLAG"),"\n";
my %childToParentError;
while(<KIN2>) {
	
	chomp;
	
	my ($fid1, $iid1, $fid2, $iid2, $phi, $pedrel, $snprel, $child, $parent) = split("\t", $_);
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
	
	if ($pedrel eq "PO" && $snprel eq "UN") {
		
		if (exists $inferredPO->{$parent}) {
			#print join("\t", $uniq_child, $parent, "Likely_child_swap"),"\n";
			next; 
		} else {
			#print join("\t", $uniq_child, $parent, "Likely_wrong_parent"),"\n";
			print join("\t", $uniq_child, $parent),"\n";
			if (exists $childToParentError{$uniq_child}) {
				push @{$childToParentError{$uniq_child}}, $parent;
			} else {
				$childToParentError{$uniq_child}[0] = $parent;
			}
			
		}
		
	}
	
}
close(KIN2);



#my $fam = "relatedqc7.fam";
#open(FAM, "<$fam") or die("Cannot open $fam\n");
#while(<FAM>) {
#	chomp;
#	
#	my ($fid, $iid, $fat, $mot, $sex, $case) = split("\t", $_);
#	my $uniq = join("\t", $fid, $iid);
#	
#	if (exists $childToParentError{$uniq}) {
#		
#		my $parent1 = $childToParentError{$iid}[0];
#		#my $parent2 = $childToParentError{$iid}[1];
#
#		if ($parent1 eq $fat) {
#			print joint($fid, $iid, "0", $mot),"\n";	
#				
#		}
#			
#	} 
#		
#}
#close(FAM);














