#!/bin/perl
use strict;
use warnings;

#################################################################
# Usage: perl within_family_swaps.pl {mkin} {sex_map} {pedfile}
#################################################################


#INPUT FILES
my $mkin = $ARGV[0];
my $sexmap = $ARGV[1];
my $pedfile = $ARGV[2];

#OUTPUT FILES
my $parentswaps = "within_family_parent_swaps.txt";
my $siblingswaps = "within_family_sibling_swaps.txt";
my $pcswaps = "within_family_parentchild_swaps.txt";
my $sexupdate = "within_family_swaps_sexupdate.txt";
my $parentupdate = "within_family_swaps_parentupdate.txt";
my $categories = "family_categories.txt";
my $exclusions = "families_to_exclude.txt";


#Define expected pedigrees
my %iidToFID; 
my %ExpPedigree;
my %Mothers;
my %Fathers;
my %iidToFather;
my %iidToMother;
open(PED, "<$pedfile") or die("Count not open $pedfile\n");
while(<PED>) {
	my ($fid, $iid, $fat, $mot, $sex, $case) = split(/\s+/, $_);
	
	if (exists($iidToFID{$iid})) {die("ERROR duplicate individual id\n")}
	$iidToFID{$iid} = $fid;
	
	#Create hash of all subjects in each family
	$ExpPedigree{$fid}{'ALL'}{$iid}=1;
	#$ExpPedigree{$fid}{'ALL'}{$fat}=1;
	#$ExpPedigree{$fid}{'ALL'}{$mot}=1;

	$iidToFather{$iid}=$fat;
	$iidToMother{$iid}=$mot;
	
	#Create hash of all mothers and fathers 	
	if ($fat ne "0") {
		$Fathers{$fat}=1;
	}
	if ($mot ne "0") {
		$Mothers{$mot}=1;
	}
	
}
close(PED);
foreach my $fid (keys %ExpPedigree) {
	
	foreach my $id (keys %{$ExpPedigree{$fid}{'ALL'}}) {
		if (exists($Fathers{$id})) {
			$ExpPedigree{$fid}{'FAT'}{$id}=1;
			$ExpPedigree{$fid}{'PARENTS'}{$id}=1;
			next;
		} 
		if (exists($Mothers{$id})) {
			$ExpPedigree{$fid}{'MOT'}{$id}=1;
			$ExpPedigree{$fid}{'PARENTS'}{$id}=1;
			next;
		}
		if (!exists($Fathers{$id}) && !exists($Mothers{$id})) {
			$ExpPedigree{$fid}{'CHILDREN'}{$id}=1;
			next;
		}
		
	}
}

#Define SNP-based relationships
my %snpPO;
my %snpFS;
my %FamiliesWithRelErrors;
open(MKIN, "<$mkin") or die("Cannot open $mkin\n");
while(<MKIN>) {
	
	chomp;
	if ($_ =~/^FID1/) {
		next;
	}
	my ($fid1, $iid1, $fid2, $iid2, $phi, $pedrel, $snprel, $child,  $parent, $cohort1, $cohort2) = split("\t", $_);
	
	if ($fid1 ne $fid2) {
		die("MKIN contains between family relationships\n");
	}
		
	if ($snprel eq "PO") {
		$snpPO{$iid1}{$iid2}=1;
		$snpPO{$iid2}{$iid1}=1;
	} elsif ($snprel eq "FS") {
		$snpFS{$iid1}{$iid2}=1;
		$snpFS{$iid2}{$iid1}=1;
	}
	
	if ($pedrel ne $snprel) {
		$FamiliesWithRelErrors{$fid1}=1;
	}

}
close(MKIN);

#Define PED-based and SNP-based sex
my %sexIidToSex;
my %FamiliesWithSexErrors;
open(SEXMAP,"<$sexmap") or die("Cannot open $sexmap\n");
while(<SEXMAP>) {
	
	chomp;
	if ($_ =~/^FID/) {
		next;
	}
	my ($fid, $iid, $pedsex, $snpsex) = split(/\s+/, $_);
	
	#if sex is missing in pedigree, we assume there are no sex errors
	if ($pedsex eq 0) {
		$pedsex = $snpsex;
	}
	
	$sexIidToSex{$iid}{'PEDSEX'} = $pedsex;
	$sexIidToSex{$iid}{'SNPSEX'} = $snpsex;
	
	if ($pedsex ne $snpsex) {
		$FamiliesWithSexErrors{$fid}=1;
	}
	
}
close(SEXMAP);



#Go through each family and decipher most likely swap
open(CAT, ">$categories") or die("Cannot open $categories\n");
open(EXCLUDE, ">$exclusions") or die("Cannot open $exclusions\n");
my %P_swaps;
my %S_swaps;
my %PC_swaps;
my %FamiliesToExclude;
my @allfamilies = keys %ExpPedigree;
foreach my $fid (@allfamilies) {
	
	my $category = checkFamilyCategory($fid);	
	print CAT join("\t",$category, keys %{$ExpPedigree{$fid}{'ALL'}}),"\n";


	#identify families with ill-defined pedigrees
	if (exists($ExpPedigree{$fid}{'PARENTS'})) {
		#too many parents
		if (scalar(keys %{$ExpPedigree{$fid}{'PARENTS'}})>2){
			$FamiliesToExclude{$fid} = "More than two parents";
		} else {
			if (exists($ExpPedigree{$fid}{'FAT'})) {
				#too many fathers
				if (scalar(keys %{$ExpPedigree{$fid}{'FAT'}})>1) {
					$FamiliesToExclude{$fid} = "More than one father";
				#undefined father	
				} elsif (exists($ExpPedigree{$fid}{'CHILDREN'})) {
					my $fat = ( keys %{$ExpPedigree{$fid}{'FAT'}} )[0];
					my @children = keys %{$ExpPedigree{$fid}{'CHILDREN'}};	
					foreach my $child (@children) {
						if ($iidToFather{$child} ne $fat) {
							$FamiliesToExclude{$fid} = "Undefined parent";
						}
					}
				}
			} 
			if (exists($ExpPedigree{$fid}{'MOT'})) {
				#too many mothers
				if (scalar(keys %{$ExpPedigree{$fid}{'MOT'}})>1) {
					$FamiliesToExclude{$fid} = "More than one mother";
				#undefined mother	
				} elsif (exists($ExpPedigree{$fid}{'CHILDREN'})) {
					my $mot = ( keys %{$ExpPedigree{$fid}{'MOT'}} )[0];
					my @children = keys %{$ExpPedigree{$fid}{'CHILDREN'}};
					foreach my $child (@children) {
						if ($iidToMother{$child} ne $mot) {
							$FamiliesToExclude{$fid} = "Undefined parent";
						}
					}
				}
			}
		}	
	} #Also, maybe if there are siblings that have different parents??
	
	#skip family if pedigree is ill-defined
	if (exists($FamiliesToExclude{$fid})) {
		print EXCLUDE join("\t", $fid, $category, $FamiliesToExclude{$fid}),"\n";
		next;
	}
	
	######### TRIO #########
	if ($category eq "trio") {
		my $fat = ( keys %{$ExpPedigree{$fid}{'FAT'}} )[0] ;
		my $mot = ( keys %{$ExpPedigree{$fid}{'MOT'}} )[0] ;
		my $child = ( keys %{$ExpPedigree{$fid}{'CHILDREN'}} )[0];	
		
		my $relerrors = 0;
		if (!exists($snpPO{$fat}{$child})) {$relerrors++};
		if (!exists($snpPO{$mot}{$child})) {$relerrors++};
		if (exists($snpPO{$fat}{$mot})) {$relerrors++};
		if (exists($snpFS{$fat}{$mot})) {$relerrors++};

		#parent swaps
		if ($relerrors==0) {
			if ($sexIidToSex{$fat}{'SNPSEX'}==2 && $sexIidToSex{$mot}{'SNPSEX'}==1) {
				$P_swaps{$fat} = $mot;
				$P_swaps{$mot} = $fat;
				#print P_SWAPS join("\t", $category, $fid, $fat, $fid, $mot),"\n";
				#print P_SWAPS join("\t", $category, $fid, $mot, $fid, $fat),"\n";	
			}
		}
		
		#parent-child swaps
		if ($relerrors>0) {
			#mot<-->child
			if (!exists($snpPO{$fat}{$child}) && exists($snpPO{$fat}{$mot}) && exists($snpPO{$mot}{$child})) {
				if ($sexIidToSex{$mot}{'SNPSEX'}==1 && $sexIidToSex{$child}{'SNPSEX'}==2) {
					$PC_swaps{$mot} = $child;
					$PC_swaps{$child} = $mot;
					#print P_SWAPS join("\t", $category, $fid, $mot, $fid, $child),"\n";
					#print P_SWAPS join("\t", $category, $fid, $child, $fid, $mot),"\n";
				}	
			}	
			#fat<-->child
			if (!exists($snpPO{$mot}{$child}) && exists($snpPO{$mot}{$fat}) && exists($snpPO{$fat}{$child})) {
				if ($sexIidToSex{$fat}{'SNPSEX'}==2 && $sexIidToSex{$child}{'SNPSEX'}==1) {
					$PC_swaps{$fat} = $child;
					$PC_swaps{$child} = $fat;
					#print PC_SWAPS join("\t", $category, $fid, $fat, $fid, $child),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child, $fid, $fat),"\n";
				}	
			}
		}
		
		
	######### TRIOPLUS	#########
	} elsif ($category eq "trioplus") {
		my $fat = ( keys %{$ExpPedigree{$fid}{'FAT'}} )[0] ;
		my $mot = ( keys %{$ExpPedigree{$fid}{'MOT'}} )[0] ;
		my @children = keys %{$ExpPedigree{$fid}{'CHILDREN'}};
		my @childpairs = chooseTwo(@children);	

		my $relerrors = 0;
		if (!exists($snpPO{$fat}{$children[0]})) {$relerrors++};
		if (exists($snpFS{$fat}{$children[0]})) {$relerrors++};
		if (!exists($snpPO{$mot}{$children[0]})) {$relerrors++};
		if (exists($snpFS{$mot}{$children[0]})) {$relerrors++};
		foreach my $pair (@childpairs) {
			my ($child1, $child2) = split(":", $pair);		
			if (!exists($snpFS{$child1}{$child2})) {$relerrors++};
			if (exists($snpPO{$child1}{$child2})) {$relerrors++};
		}
		
		#parent swaps
		if ($relerrors==0) {
			if ($sexIidToSex{$fat}{'SNPSEX'}==2 && $sexIidToSex{$mot}{'SNPSEX'}==1) {
				$P_swaps{$fat} = $mot;
				$P_swaps{$mot} = $fat;
				#print P_SWAPS join("\t", $category, $fid, $fat, $fid, $mot),"\n";
				#print P_SWAPS join("\t", $category, $fid, $mot, $fid, $fat),"\n";	
			}
		}

		#sibling swaps
		if ($relerrors==0) {
			foreach my $pair (@childpairs) {
				my ($child1, $child2) = split(":", $pair);
				
				my $sexerrors = 0;
				if ($sexIidToSex{$child1}{'PEDSEX'} ne $sexIidToSex{$child1}{'SNPSEX'}) {$sexerrors++}; 
				if ($sexIidToSex{$child2}{'PEDSEX'} ne $sexIidToSex{$child2}{'SNPSEX'}) {$sexerrors++}; 
				if ($sexIidToSex{$child1}{'PEDSEX'} ne $sexIidToSex{$child2}{'SNPSEX'}) {$sexerrors++};
				if ($sexIidToSex{$child2}{'PEDSEX'} ne $sexIidToSex{$child1}{'SNPSEX'}) {$sexerrors++};
				
				if ($sexerrors==4) {
					$S_swaps{$child1} = $child2;
					$S_swaps{$child2} = $child1;
					#print S_SWAPS join("\t", $category, $fid, $child1, $fid, $child2),"\n";
					#print S_SWAPS join("\t", $category, $fid, $child2, $fid, $child1),"\n";
				} 
						
			}
		}
				
		#parent-child swaps
		if ($relerrors>0) {
			foreach my $pair (@childpairs) {
				my ($child1, $child2) = split(":", $pair);
				#mot<-->child1
				if (exists($snpPO{$mot}{$child1}) && !exists($snpPO{$fat}{$child1}) && exists($snpFS{$mot}{$child2})) {
					$PC_swaps{$mot} = $child1;
					$PC_swaps{$child1} = $mot;
					#print PC_SWAPS join("\t", $category, $fid, $mot, $fid, $child1),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child1, $fid, $mot),"\n";
				}	
				#fat<-->child1
				if (exists($snpPO{$fat}{$child1}) && !exists($snpPO{$mot}{$child1}) && exists($snpFS{$fat}{$child2})) {
					$PC_swaps{$fat} = $child1;
					$PC_swaps{$child1} = $fat;
					#print PC_SWAPS join("\t", $category, $fid, $fat, $fid, $child1),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child1, $fid, $fat),"\n";
				}
				#mot<-->child2
				if (exists($snpPO{$mot}{$child2}) && !exists($snpPO{$fat}{$child2}) && exists($snpFS{$mot}{$child1})) {
					$PC_swaps{$mot} = $child2;
					$PC_swaps{$child2} = $mot;
					print PC_SWAPS join("\t", $category, $fid, $mot, $fid, $child2),"\n";
					print PC_SWAPS join("\t", $category, $fid, $child2, $fid, $mot),"\n";
				}	
				#fat<-->child2
				if (exists($snpPO{$fat}{$child2}) && !exists($snpPO{$mot}{$child2}) && exists($snpFS{$fat}{$child1})) {
					$PC_swaps{$fat} = $child2;
					$PC_swaps{$child2} = $fat;
					#print PC_SWAPS join("\t", $category, $fid, $fat, $fid, $child2),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child2, $fid, $fat),"\n";
				}
			}
		}
		
		
	######### PO  #########
	} elsif ($category eq "po") {
		my $parent;
		if (exists($ExpPedigree{$fid}{'FAT'})) {
			$parent = ( keys %{$ExpPedigree{$fid}{'FAT'}} )[0] ;
		} elsif (exists($ExpPedigree{$fid}{'MOT'})) {
			$parent = ( keys %{$ExpPedigree{$fid}{'MOT'}} )[0] ;
		}
		my $child = ( keys %{$ExpPedigree{$fid}{'CHILDREN'}} )[0];

		my $relerrors = 0;
		if (!exists($snpPO{$parent}{$child})) {$relerrors++};
		if (exists($snpFS{$parent}{$child})) {$relerrors++};

		#parent-child swaps
		if ($relerrors==0) {
			
			my $sexerrors = 0;
			if ($sexIidToSex{$parent}{'PEDSEX'} ne $sexIidToSex{$parent}{'SNPSEX'}) {$sexerrors++};
			if ($sexIidToSex{$child}{'PEDSEX'} ne $sexIidToSex{$child}{'SNPSEX'}) {$sexerrors++};
			if ($sexerrors>0) {
				if ($sexIidToSex{$child}{'PEDSEX'} eq $sexIidToSex{$parent}{'SNPSEX'}) {
					$PC_swaps{$parent} = $child;
					$PC_swaps{$child} = $parent;
					#print PC_SWAPS join("\t", $category, $fid, $parent, $fid, $child),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child, $fid, $parent),"\n";
				}
			}
			 
		}
		
		
	######### POPLUS #########
	} elsif ($category eq "poplus") {
		my $parent;
		if (exists($ExpPedigree{$fid}{'FAT'})) {
			$parent = ( keys %{$ExpPedigree{$fid}{'FAT'}} )[0] ;
		} elsif (exists($ExpPedigree{$fid}{'MOT'})) {
			$parent = ( keys %{$ExpPedigree{$fid}{'MOT'}} )[0] ;
		}
		my @children = keys %{$ExpPedigree{$fid}{'CHILDREN'}};		
		my @childpairs = chooseTwo(@children);
		
		my $relerrors = 0;
		if (!exists($snpPO{$parent}{$children[0]})) {$relerrors++};
		if (exists($snpFS{$parent}{$children[0]})) {$relerrors++};
		foreach my $pair (@childpairs) {
			my ($child1, $child2) = split(":", $pair);		
			if (!exists($snpFS{$child1}{$child2})) {$relerrors++};
			if (exists($snpPO{$child1}{$child2})) {$relerrors++};
		}
		
		#parent-child swaps
		if ($relerrors>0) {
			foreach my $pair (@childpairs) {
				my ($child1, $child2) = split(":", $pair);
				#parent<-->child1
				if (exists($snpPO{$parent}{$child1}) && exists($snpPO{$child1}{$child2}) && exists($snpFS{$parent}{$child2})) {
					$PC_swaps{$parent} = $child1;
					$PC_swaps{$child1} = $parent;
					#print PC_SWAPS join("\t", $category, $fid, $parent, $fid, $child1),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child1, $fid, $parent),"\n";
				}	
				#parent<-->child2
				if (exists($snpPO{$parent}{$child2}) && exists($snpPO{$child1}{$child2}) && exists($snpFS{$parent}{$child1})) {
					$PC_swaps{$parent} = $child2;
					$PC_swaps{$child2} = $parent;
					#print PC_SWAPS join("\t", $category, $fid, $parent, $fid, $child2),"\n";
					#print PC_SWAPS join("\t", $category, $fid, $child2, $fid, $parent),"\n";
				}	
			}
		}

		#sibling swaps
		if ($relerrors==0) {
			foreach my $pair (@childpairs) {
				my ($child1, $child2) = split(":", $pair);
				
				my $sexerrors = 0;
				if ($sexIidToSex{$child1}{'PEDSEX'} ne $sexIidToSex{$child1}{'SNPSEX'}) {$sexerrors++}; 
				if ($sexIidToSex{$child2}{'PEDSEX'} ne $sexIidToSex{$child2}{'SNPSEX'}) {$sexerrors++}; 
				if ($sexIidToSex{$child1}{'PEDSEX'} ne $sexIidToSex{$child2}{'SNPSEX'}) {$sexerrors++};
				if ($sexIidToSex{$child2}{'PEDSEX'} ne $sexIidToSex{$child1}{'SNPSEX'}) {$sexerrors++};
								
				if ($sexerrors==4) {
					$S_swaps{$child1} = $child2;
					$S_swaps{$child2} = $child1;
					#print S_SWAPS join("\t", $category, $fid, $child1, $fid, $child2),"\n";
					#print S_SWAPS join("\t", $category, $fid, $child2, $fid, $child1),"\n";
				} 
			}		
		
		}
		
	######### SIBS	#########
	} elsif ($category eq "sibs") {
		my @children = keys %{$ExpPedigree{$fid}{'CHILDREN'}};	
		my @childpairs = chooseTwo(@children);
		
		my $relerrors = 0;
		foreach my $pair (@childpairs) {
			my ($child1, $child2) = split(":", $pair);		
			if (!exists($snpFS{$child1}{$child2})) {$relerrors++};
			if (exists($snpPO{$child1}{$child2})) {$relerrors++};
		}

		#sibling swaps
		if ($relerrors==0) {
			my @childpairs = chooseTwo(@children);
			foreach my $pair (@childpairs) {
				my ($child1, $child2) = split(":", $pair);
				
				my $sexerrors = 0;
				if ($sexIidToSex{$child1}{'PEDSEX'} ne $sexIidToSex{$child1}{'SNPSEX'}) {$sexerrors++}; 
				if ($sexIidToSex{$child2}{'PEDSEX'} ne $sexIidToSex{$child2}{'SNPSEX'}) {$sexerrors++}; 
				if ($sexIidToSex{$child1}{'PEDSEX'} ne $sexIidToSex{$child2}{'SNPSEX'}) {$sexerrors++};
				if ($sexIidToSex{$child2}{'PEDSEX'} ne $sexIidToSex{$child1}{'SNPSEX'}) {$sexerrors++};
				
				if ($sexerrors==4) {
					$S_swaps{$child1} = $child2;
					$S_swaps{$child2} = $child1;					
					#print S_SWAPS join("\t", $category, $fid, $child1, $fid, $child2),"\n";
					#print S_SWAPS join("\t", $category, $fid, $child2, $fid, $child1),"\n";
				} 
						
			}
		}
	
		
	} elsif ($category eq "unrel") {
		next;
		
	} else {
		print "ERROR - $fid was not assigned a category\n";
	}

}
close(EXCLUDE);
close(CAT);


#Print parent swaps to file 
open(P_SWAPS, ">$parentswaps") or die("Cannot open $parentswaps\n");
foreach my $id1 (keys %P_swaps) {
	my $id2 = $P_swaps{$id1};
	my $fid = $iidToFID{$id1};
	my $category = checkFamilyCategory($fid);
	print P_SWAPS join("\t", $category, $fid, $id1, $fid, $id2),"\n";
}
close(P_SWAPS);


#Print sibling swaps to file
open(S_SWAPS, ">$siblingswaps") or die("Cannot open $siblingswaps\n");
foreach my $id1 (keys %S_swaps) {
	my $id2 = $S_swaps{$id1};
	my $fid = $iidToFID{$id1};
	my $category = checkFamilyCategory($fid);
	print S_SWAPS join("\t", $category, $fid, $id1, $fid, $id2),"\n";
}
close(S_SWAPS);

#Print parent-child swaps to file
open(PC_SWAPS, ">$pcswaps") or die("Cannot open $pcswaps\n");
foreach my $id1 (keys %PC_swaps) {
	my $id2 = $PC_swaps{$id1};
	my $fid = $iidToFID{$id1};
	my $category = checkFamilyCategory($fid);
	print PC_SWAPS join("\t", $category, $fid, $id1, $fid, $id2),"\n";
}
close(PC_SWAPS);


#Create sex update file
open(my $SEXUPDATE,">$sexupdate") or die("Cannot open $sexupdate\n");
getSexUpdate($parentswaps, $SEXUPDATE);
getSexUpdate($pcswaps, $SEXUPDATE);
getSexUpdate($siblingswaps, $SEXUPDATE);
close($SEXUPDATE);

#Create parent update file
open(my $PARENTUPDATE,">$parentupdate") or die("Cannot open $parentupdate\n");
getParentUpdate($parentswaps, $PARENTUPDATE);
getParentUpdate($pcswaps, $PARENTUPDATE);
getParentUpdate($siblingswaps, $PARENTUPDATE);
close($PARENTUPDATE);


sub getSexUpdate {
	my $swapfile = $_[0];
	my $out = $_[1];
	open(IN, "<$swapfile") or die("Cannot open $swapfile\n");
	while(<IN>) {
		chomp;
		my ($category, $fid1, $iid1, $fid2, $iid2) = split(/\s+/, $_);
		print $out join("\t", $fid2, $iid2, $sexIidToSex{$iid2}{'PEDSEX'}),"\n";	
	}
	close(IN);
}

sub getParentUpdate {
	my $swapfile = $_[0];
	my $out = $_[1];
	open(IN, "<$swapfile") or die("Cannot open $swapfile\n");
	while(<IN>) {
		chomp;
		my ($category, $fid1, $iid1, $fid2, $iid2) = split(/\s+/, $_);
		print $out join("\t", $fid2, $iid2, $iidToFather{$iid2}, $iidToMother{$iid2}),"\n";	
	}
	close(IN);
}

sub checkFamilyCategory {
	my $fid = $_[0];
	my $category;
	if (exists($ExpPedigree{$fid}{'FAT'}) && exists($ExpPedigree{$fid}{'MOT'}) && exists($ExpPedigree{$fid}{'CHILDREN'})) {
		my $no_children = scalar(keys %{$ExpPedigree{$fid}{'CHILDREN'}});
		if ($no_children==1) {
			$category = "trio";
		} elsif ($no_children>1) {
			$category = "trioplus";
		}
		
	} elsif (exists($ExpPedigree{$fid}{'CHILDREN'}) && exists($ExpPedigree{$fid}{'PARENTS'})) {
		my $no_children = scalar(keys %{$ExpPedigree{$fid}{'CHILDREN'}});
		if ($no_children==1) {
			$category = "po";
		} elsif ($no_children>1) {
			$category = "poplus";
		}
		
	} elsif (exists($ExpPedigree{$fid}{'CHILDREN'})) {
		my $no_children = scalar(keys %{$ExpPedigree{$fid}{'CHILDREN'}});
		if ($no_children==1) {
			$category = "unrel";
		} elsif ($no_children>1) {
			$category = "sibs";
		}
	} else {
		$category = "ERROR";
	}
	return $category;
}

sub chooseTwo {
	my @list = @_;
	my %Pairs;
	foreach my $item1 (@list) {
		foreach my $item2 (@list) {
			if ($item1 ne $item2) {
				my $pair = "$item1:$item2";
				my $rpair = "$item2:$item1";
				if (!exists($Pairs{$pair}) && !exists($Pairs{$rpair})) {
					$Pairs{$pair}=1;
				}	
			}
		}
	}
	return keys %Pairs;	
}

sub printExpFam {
	my $fam = $_[0];
	print join("\t","PED-BASED FAMILY:"),"\n";
	print join("\t","FAMID", $fam),"\n";
	print join("\t", 'PARENTS', keys %{$ExpPedigree{$fam}{'PARENTS'}}),"\n";
	print join("\t", 'FATHER', keys %{$ExpPedigree{$fam}{'FAT'}}),"\n";
	print join("\t", 'MOTHER', keys %{$ExpPedigree{$fam}{'MOT'}}),"\n";
	print join("\t", 'CHILDREN', keys %{$ExpPedigree{$fam}{'CHILDREN'}}),"\n";
	print "\n";
}
#printExpFam('7038');
