
set -euxo pipefail

comment "Setting up"

##############################################################################################
#
#	Build directory structure
#
##############################################################################################
echo "Building directory structure"

mkdir -p ${preqc_folder}
mkdir -p ${scripts} ${data_derived} ${data_raw} ${clean} ${prelim} ${logs}
mkdir -p ${qclists} ${sexqc} ${relqc} ${snpqc}


 ##############################################################################################
#
#	Retrieve raw files, copy to data
#
##############################################################################################
echo "Retrieving raw files"

#make directory writable
chmod -R u=rwx ${data_raw}/ 

#copy raw data
cp --recursive ${PROJECT_MEGA}/raw_data/* ${data_raw}/
#make directory writable
chmod -R u=rwx ${data_raw}/*

#copy other data
cp ${PROJECT_MEGA}/probe_alignment/mymapInfo ${data_raw}/snpqc_files/

