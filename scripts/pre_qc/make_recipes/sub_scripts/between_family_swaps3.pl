#!/bin/perl
use strict;
use warnings;


################################################################
#	Usage: perl between_family_swaps2.pl {mkin} {family_swaps} {sexmap} {pedfile}

#INPUT FILES
my $mkin = $ARGV[0];
my $swaps = $ARGV[1];
my $sexmap = $ARGV[2];
my $pedfile = $ARGV[3];

#OUTPUT FILES
my $pcswaps = "between_family_parent_child_swaps.txt";
my $childswaps = "between_family_child_swaps.txt";
my $parentswaps = "between_family_parent_swaps.txt";
my $sexupdate = "between_family_swaps_sexupdate.txt";
my $parentupdate = "between_family_swaps_parentupdate.txt";
my $exclusions = "families_skipped.txt";

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
open(MKIN, "<$mkin") or die("Cannot open $mkin\n");
while(<MKIN>) {
	
	chomp;
	my ($fid1, $iid1, $fid2, $iid2, $phi, $pedrel, $snprel, $child,  $parent, $cohort1, $cohort2) = split("\t", $_);
		
	if ($snprel eq "PO") {
		$snpPO{$iid1}{$iid2}=1;
		$snpPO{$iid2}{$iid1}=1;
	} elsif ($snprel eq "FS") {
		$snpFS{$iid1}{$iid2}=1;
		$snpFS{$iid2}{$iid1}=1;
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


#Get pairs of families with unexpected relatedness
my %uniqSwapPairs;
open(SWAPS,"<$swaps") or die("Cannot open $swaps\n");
while(<SWAPS>) {
	chomp;
	my ($fid1, $fid2) = split(/\s+/, $_);
	if (!exists($uniqSwapPairs{$fid1."\t".$fid2}) && !exists($uniqSwapPairs{$fid2."\t".$fid1})) {
		$uniqSwapPairs{$fid1."\t".$fid2}=1;
		next;
	} 
}
close(SWAPS);

#Identify families with ill-defined pedigrees
my %FamiliesToExclude;
my @allfamilies = keys %ExpPedigree;
foreach my $fid (@allfamilies) {

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
	} 
}
#Also, maybe if there are siblings that have different parents??



#Go through each pair of families and define most likely swap
my %trackPairSwaps;
open(EXCLUDE, ">$exclusions") or die("Cannot open $exclusions\n");
my %P_swaps;
my %C_swaps;
my %PC_swaps;
my @pairs = keys %uniqSwapPairs;
foreach my $pair (@pairs) {
	
	my ($fid1, $fid2) = split("\t", $pair); 
	my $category1 = checkFamilyCategory($fid1);
	my $category2 = checkFamilyCategory($fid2);
	
	if (exists($FamiliesToExclude{$fid1}) || exists($FamiliesToExclude{$fid2})) {
		my $status1;
		if (exists($FamiliesToExclude{$fid1})) {
			$status1 = $FamiliesToExclude{$fid1};
		} else {
			$status1 = "OK";
		}
		my $status2;
		if (exists($FamiliesToExclude{$fid2})) {
			$status2 = $FamiliesToExclude{$fid2};
		} else {
			$status2 = "OK";
		}
		print EXCLUDE join("\t", $fid1, $category1, $status1, $fid2, $category2, $status2),"\n";
		next;
	}
	
	############# TWO TRIOS #############
	if (($category1 eq "trio" || $category1 eq "trioplus") && ($category2 eq "trio" || $category2 eq "trioplus")) {
				
		#DEFINE EXPECTED PEDIGREES FOR EACH FAMILY
		my $fat1 = ( keys %{$ExpPedigree{$fid1}{'FAT'}} )[0] ;
		my $mot1 = ( keys %{$ExpPedigree{$fid1}{'MOT'}} )[0] ;
		my @children1 = keys %{$ExpPedigree{$fid1}{'CHILDREN'}};	
		my $fat2 = ( keys %{$ExpPedigree{$fid2}{'FAT'}} )[0];
		my $mot2 = ( keys %{$ExpPedigree{$fid2}{'MOT'}} )[0];
		my @children2 = keys %{$ExpPedigree{$fid2}{'CHILDREN'}};
		my $no_children1 = scalar @children1;
		my $no_children2 = scalar @children2;
		print join("\t", $category1, $no_children1, $category2, $no_children2),"\n";
		
		
		#CHILD SWAPS
		#children1[i] <--> children2[j]
		foreach my $child1 (@children1) {
			my $relerrors1 = 0;
			if (!exists($snpPO{$fat1}{$child1})) {$relerrors1++};
			if (!exists($snpPO{$mot1}{$child1})) {$relerrors1++};
			if ($relerrors1==2) {
				foreach my $child2 (@children2) {
					my $relerrors2 = 0;
					if (!exists($snpPO{$fat2}{$child2})) {$relerrors2++};
					if (!exists($snpPO{$mot2}{$child2})) {$relerrors2++};
					if ($relerrors2==2) {
						my $swaperrors = 0;
						if (exists($snpPO{$fat2}{$child1})) {$swaperrors++};
						if (exists($snpPO{$mot2}{$child1})) {$swaperrors++};
						if (exists($snpPO{$fat1}{$child2})) {$swaperrors++};
						if (exists($snpPO{$mot1}{$child2})) {$swaperrors++};
						if ($swaperrors==4) {
							$C_swaps{$child1} = $child2;
							$C_swaps{$child2} = $child1;
							$trackPairSwaps{$pair} = "child_swap";
						}
					}
				}	
			}	
		}


		#PARENT SWAPS 
		#parent1 <--> parent2
		if (!exists($trackPairSwaps{$pair})) {
			my @parentpairs = ("$fat1:$fat2", "$mot1:$mot2", "$fat1:$mot2", "$mot1:$fat2");
			foreach my $pair (@parentpairs) {						
				my ($parent1, $parent2) = split(":", $pair);
				my $errors1 = 0;
				foreach my $child (@children1) {
					if (!exists($snpPO{$parent1}{$child})) {$errors1++}; 
				}
				my $errors2 = 0;
				foreach my $child (@children2) {
					if (!exists($snpPO{$parent2}{$child})) {$errors2++}; 
				} 
				 
				my $swap1 = 0;
				foreach my $child (@children2) {
					if (exists($snpPO{$parent1}{$child})) {$swap1++};
				}
				my $swap2 = 0;
				foreach my $child (@children1) {
					if (exists($snpPO{$parent2}{$child})) {$swap2++};
				}
				
				if ($errors1==$no_children1 && $errors2==$no_children2 && $swap1==$no_children2 && $swap2==$no_children1) {
					$P_swaps{$parent1} = $parent2;
					$P_swaps{$parent2} = $parent1;		
					$trackPairSwaps{$pair} = "parent_swap";
				}
			}
		}		
		
		#PARENT-CHILD SWAPS
		if (!exists($trackPairSwaps{$pair})) {
			#fat1 <--> child2 						
			if (exists($snpPO{$fat1}{$fat2}) && exists($snpPO{$fat1}{$mot2})) {	
				foreach my $child (@children2) {
					if (exists($snpPO{$child}{$children1[0]})) {
						$PC_swaps{$fat1} = $child;
						$PC_swaps{$child} = $fat1;
					}
				}
			}
			#mot1 <--> child2
			if (exists($snpPO{$mot1}{$fat2}) && exists($snpPO{$mot1}{$mot2})) {
				foreach my $child (@children2) {
					if (exists($snpPO{$child}{$children1[0]})) {
						$PC_swaps{$mot1} = $child;
						$PC_swaps{$child} = $mot1;
					}
				}
			}
			#fat2 <--> child1
			if (exists($snpPO{$fat2}{$fat1}) && exists($snpPO{$fat2}{$mot1})) {
				foreach my $child (@children1) {
					if (exists($snpPO{$child}{$children2[0]})) {
						$PC_swaps{$fat2} = $child;
						$PC_swaps{$child} = $fat2;
					}
				}
			}
			#mot2 <--> child1
			if (exists($snpPO{$mot2}{$fat1}) && exists($snpPO{$mot2}{$mot1})) {
				foreach my $child (@children1) {
					if (exists($snpPO{$child}{$children2[0]})) {
						$PC_swaps{$mot2} = $child;
						$PC_swaps{$child} = $mot2;
					}
				}
			}	
		}		
	}
} 
close(EXCLUDE);


#Print parent swaps to file 
open(P_SWAPS, ">$parentswaps") or die("Cannot open $parentswaps\n");
foreach my $id1 (keys %P_swaps) {
	my $id2 = $P_swaps{$id1};
	my $fid1 = $iidToFID{$id1};
	my $fid2 = $iidToFID{$id2};
	my $category1 = checkFamilyCategory($fid1);
	my $category2 = checkFamilyCategory($fid2);
	print P_SWAPS join("\t", $category1, $fid1, $id1, $category2, $fid2, $id2),"\n";
}
close(P_SWAPS);


#Print child swaps to file
open(C_SWAPS, ">$childswaps") or die("Cannot open $childswaps\n");
foreach my $id1 (keys %C_swaps) {
	my $id2 = $C_swaps{$id1};
	my $fid1 = $iidToFID{$id1};
	my $fid2 = $iidToFID{$id2};
	my $category1 = checkFamilyCategory($fid1);
	my $category2 = checkFamilyCategory($fid2);
	print C_SWAPS join("\t", $category1, $fid1, $id1, $category2, $fid2, $id2),"\n";
}
close(C_SWAPS);

#Print parent-child swaps to file
open(PC_SWAPS, ">$pcswaps") or die("Cannot open $pcswaps\n");
foreach my $id1 (keys %PC_swaps) {
	my $id2 = $PC_swaps{$id1};
	my $fid1 = $iidToFID{$id1};
	my $fid2 = $iidToFID{$id2};
	my $category1 = checkFamilyCategory($fid1);
	my $category2 = checkFamilyCategory($fid2);
	print PC_SWAPS join("\t", $category1, $fid1, $id1, $category2, $fid2, $id2),"\n";
}
close(PC_SWAPS);



#Create sex update file
open(my $SEXUPDATE,">$sexupdate") or die("Cannot open $sexupdate\n");
getSexUpdate($parentswaps, $SEXUPDATE);
getSexUpdate($pcswaps, $SEXUPDATE);
getSexUpdate($childswaps, $SEXUPDATE);
close($SEXUPDATE);

#Create parent update file
open(my $PARENTUPDATE,">$parentupdate") or die("Cannot open $parentupdate\n");
getParentUpdate($parentswaps, $PARENTUPDATE);
getParentUpdate($pcswaps, $PARENTUPDATE);
getParentUpdate($childswaps, $PARENTUPDATE);
close($PARENTUPDATE);



sub getSexUpdate {
	my $swapfile = $_[0];
	my $out = $_[1];
	open(IN, "<$swapfile") or die("Cannot open $swapfile\n");
	while(<IN>) {
		chomp;
		my ($category1, $fid1, $iid1, $category2, $fid2, $iid2) = split(/\s+/, $_);
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
		my ($category1, $fid1, $iid1, $category2, $fid2, $iid2) = split(/\s+/, $_);
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