#!/bin/perl
use strict;
use warnings;

my $fatupdate = $ARGV[0];
my $motupdate = $ARGV[1];
my $child_ids =  $ARGV[2]; ;

my $iidToFatUpdate = {};
open(FAT, "<$fatupdate") or die("Cannot open $fatupdate\n");
while(<FAT>) {
	chomp;
	
	my ($fid, $iid, $fat, $mot) = split(/\s+/, $_);
	
	$iidToFatUpdate -> {$iid} -> {'FATCOL'} = $fat;
	$iidToFatUpdate -> {$iid} -> {'MOTCOL'} = $mot;
 	
}
close(FAT);

my $iidToMotUpdate = {};
open(MOT, "<$motupdate") or die("Cannot open $motupdate\n");
while(<MOT>) {
	chomp;
	
	my ($fid, $iid, $fat, $mot) = split(/\s+/, $_);
	
	$iidToMotUpdate -> {$iid} -> {'FATCOL'} = $fat;
	$iidToMotUpdate -> {$iid} -> {'MOTCOL'} = $mot;
	
}
close(MOT);


open(CHILD, "<$child_ids") or die("Cannot open $child_ids\n");
while(<CHILD>) {
	chomp;
	my ($fid, $iid) = split(/\s+/, $_);
	
	if (exists $iidToFatUpdate -> {$iid} && exists $iidToMotUpdate -> {$iid}) {
		#print "BOTH","\n";
		print join("\t", $fid, $iid, $iidToFatUpdate -> {$iid}->{'FATCOL'}, $iidToMotUpdate -> {$iid}->{'MOTCOL'}),"\n";
	}
	
	elsif (exists $iidToFatUpdate -> {$iid} ) {
		#print "FAT only","\n";
		print join("\t", $fid, $iid, $iidToFatUpdate -> {$iid}->{'FATCOL'}, $iidToFatUpdate -> {$iid}->{'MOTCOL'}),"\n";
		
	}
	
	else {
		#print "MOT only","\n";
		print join("\t", $fid, $iid, $iidToMotUpdate -> {$iid}->{'FATCOL'}, $iidToMotUpdate -> {$iid}->{'MOTCOL'}),"\n";
		
	}
	
}



