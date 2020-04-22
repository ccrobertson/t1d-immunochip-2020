#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=largemem
#SBATCH --mem=100000

date

query=$1

~/bin/blat -fastMap -minIdentity=100 -out=psl ${GENOMES}/hg19/hg19.fa ${query}.fa ${query}.psl

date


