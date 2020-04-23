#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=50000

module load htslib
module load bcftools


group=$1
refpanel_nickname=$2
build=$3
chr=${SLURM_ARRAY_TASK_ID}

#Input file = unfiltered imputed data set
invcf=${release_branch}/${refpanel_nickname}/results/${group}/chr${chr}.dose.vcf.gz
outprefix=${impbench}/mega_${refpanel_nickname}_${group}_chr${chr}

if [ ! -f ${outprefix}_${build}.vcf.gz ]; then
	echo "Extracting subjects"
	bcftools view --samples-file ${impbench}/washu_batch1_mega_keeplist_${group}.txt --force-samples --output-type z --output-file ${outprefix}_${build}.vcf.gz ${invcf}
	echo "Indexing ${outprefix}_${build}.vcf.gz"
	tabix -f -p vcf ${outprefix}_${build}.vcf.gz
fi


if [ "${build}" != "gr38" ]; then

	module load picard/2.18.5

	echo "Liftover b37 to hg19"
	java -jar $EBROOTPICARD/picard.jar LiftoverVcf \
     I=${outprefix}_${build}.vcf.gz \
     O=${outprefix}_hg19.vcf.gz \
     CHAIN=${resources}/b37tohg19.chain \
     REJECT=${outprefix}_REJECT.vcf.gz \
     R=${resources}/ucsc.hg19.fasta \
     MAX_RECORDS_IN_RAM=10000


    echo "Liftover hg19 to b38"
  	java -jar $EBROOTPICARD/picard.jar LiftoverVcf \
     I=${outprefix}_hg19.vcf.gz \
     O=${outprefix}_gr38.vcf.gz \
     CHAIN=${resources}/hg19ToHg38.over.chain \
     REJECT=${outprefix}_hg19_REJECT.vcf.gz \
     R=${resources}/Homo_sapiens_assembly38.fasta \
     MAX_RECORDS_IN_RAM=10000

	echo "Indexing ${outprefix}_gr38.vcf.gz"
	tabix -f -p vcf ${outprefix}_gr38.vcf.gz

fi

echo "Extracting ichip regions"
tabix -h -R ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr38_chr${chr}.txt ${outprefix}_gr38.vcf.gz | bgzip > ${outprefix}_ichip_regions.vcf.gz

echo "Indexing ${outprefix}_ichip_regions.vcf.gz"
tabix -f -p vcf ${outprefix}_ichip_regions.vcf.gz

echo "DONE"
