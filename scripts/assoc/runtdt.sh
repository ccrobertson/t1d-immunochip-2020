#!/bin/bash

group=$1
dir=${fam_assoc}/${group}
cd ${dir}

#run in imputeded data
sbatch --array=1-22 ${scripts}/assoc/runtdt.slurm ${group}
