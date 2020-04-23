#!/usr/bin/env python3

## count_hetalleles.py :: input raw bam + normal hetsites bed :: output hetsites bed with observed (simple) allele counts

    # from within het_counter.sh
    # usage: count_hetalleles.py --bam [normal/tumor].bam --hetsites_bed normal_hetsites.bed > subject_[normal/tumor]_hetcnts.bed 2> subject_[normal/tumor]_errcnts.bed

    # this script takes a bam file and a bed file with hetsites in the normal
    # exome sequences and extracts from the bam allele counts at each hetsite;
    # then writes allele counts at each hetsite to file

    # stderr reports locations with known heterozygotes < 80% of reads


import pysam
import sys
import argparse
import re


bases = ('A','C','G','T')
field_names = ('chrom','start','stop','info')
dp_thresh = 0


def get_known_alleles(info_str):
    '''
    info_str: AF=0.500;AD::A=100,C=50
    '''
    known_alleles = []
    for info in info_str.split(';'):
        if re.match('AD::',info):
            allele_cnts = info.split('::')[1].split(',')
            allele1 = allele_cnts[0].split('=')[0]
            allele2 = allele_cnts[1].split('=')[0]
            known_alleles.append(allele1)
            known_alleles.append(allele2)
    return known_alleles


def get_alleles(bamfile, chrom, zero_based_pos):
    '''
    Return a list of alleles observed in the reads aligned at a locus of
    interest defined by chrom and one_based_pos.

    The pileup call initiates an iterator over each position (i.e. column)
    in the set of positions (columns) composed of all positions in all reads
    that overlap the position of interest specified in the pileup call.

    i.e. (1) start with all reads that overlap a given position (2) take the
    union of all unique positions in that set of reads (3) iterate over each
    of those positions and capture all reads from step (1) that overlap
    '''
    alleles_at_pos = []
    for pileupcolumn in bamfile.pileup(chrom, zero_based_pos, zero_based_pos + 1):  ## iterator over each reference position (i.e. column)
        ## pysam includes all positions from all alignments that overlap the
        ## ref position of interest, so we exclude all positions different
        ## from the one of interest
        if pileupcolumn.reference_pos != zero_based_pos:
            continue
        for pileupread in pileupcolumn.pileups:  ## iterator over each read overlapping given position (i.e. column)
            ## query_position is None if is_del or is_refskip is set
            if not pileupread.is_del and not pileupread.is_refskip:
                aligned_read_allele = pileupread.alignment.query_sequence[pileupread.query_position]
                alleles_at_pos.append(aligned_read_allele)
    return alleles_at_pos


def allele_counts(alleles, known_alleles):
    '''
    Count alleles where alleles is a list of A,C,G,T
    '''
    counts = {x:0 for x in bases}  ## reset counter
    known_total = total = 0
    for allele in alleles:
        if allele in counts:
            counts[allele] = counts[allele] + 1
            total = total + 1
        if allele in known_alleles:
            known_total = known_total + 1
    return (counts, known_total, total)


parser = argparse.ArgumentParser(description='Count alleles at heterozygous sites')
parser.add_argument('--bam',help='target bam file',action='store',dest='bam_file')
parser.add_argument('--hetsites_bed',help='positions of hetsites',action='store',dest='hetsites_bed_file')
parser.add_argument('--min_count',help='minimum number of total reads at given position',action='store',dest='min_count',type=int,default=100)
parser.add_argument('--min_score',help='minimum score value for given position',action='store',dest='min_score',type=float,default=100) ## wrp originally had 200 here but 100 in find_hetsites
args = parser.parse_args()


## open bamfile for reading bases at location
bamfile = pysam.AlignmentFile(args.bam_file,'rb')


## scan bed file of heterozygous sites
with open(args.hetsites_bed_file,'r') as in_f:
    for in_line in in_f:
        in_line = in_line.strip('\n')
        ## get het site positions
        data = dict(list(zip(field_names,in_line.split('\t'))))
        ## get known alleles, func returns a list with the major, then minor allele, but not the counts
        known_alleles = get_known_alleles(data['info'])
        ## get observed alleles in bam, func returns a list of alleles at given chrom + position
        alleles = get_alleles(bamfile,data['chrom'],int(data['start']))
        ## get observed allele counts
        (counts, known_total, all_total) = allele_counts(alleles,known_alleles)
        ## check for large counts of unknown alleles or absent locations
        if all_total > dp_thresh and known_total/all_total > 0.8:
            allele_str = 'SD::' + ','.join(['%s=%d' % (x,counts[x]) for x in known_alleles])  ## SD is simple allele depth by counting aligned reads from bam
            print(';'.join((in_line,allele_str)))
        else:
            allele_str = 'BAD::' + ','.join(['%s=%d' % (x,counts[x]) for x in bases])
            sys.stderr.write('%s;%s ***\n' % (in_line,allele_str))


## close bamfile once alleles counted
bamfile.close()




### testing

# field_names = ('chrom','start','stop','type','score','strand','info')
# bamfile = pysam.AlignmentFile('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn.bam','rb')
# with open('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn_hetsites.bed') as in_f:
#     for in_line in in_f:
#         in_line = in_line.strip('\n')
#         ## get het site positions
#         data = dict(list(zip(field_names,in_line.split('\t'))))
#         ## get known alleles, func returns a list with the major, then minor allele, but not the counts
#         known_alleles = get_known_alleles(data['info'])  ## returns known_alleles, known_dict
#         ## get observed alleles in bam, func returns a list of alleles at given chrom + position
#         alleles = get_alleles(bamfile,data['chrom'],int(data['start']))  ## this takes a long time
#         ## get observed allele counts
#         (counts, known_total, all_total) = allele_counts(alleles,known_alleles)



## notes on pileup
## given a position in the setting of 101 bp reads, the union of all positions in
## all reads overlapping the given position will extend about 100 bases down- and
## upstream of the position, resulting in about 201 columns across all reads
## overlapping a given genomic position
##
## also, nsegments will give all reads overlapping a position, including lower
## quality (at given position) read mates that are filtered out by default in
## the pileup call (unless ignore_overlaps=False); so len(pileups) != nsegments
