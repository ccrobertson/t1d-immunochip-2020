set -euxo pipefail

##############################################################################################
#
#    Validation of additions file (compare with mega_main file)
#
##############################################################################################
comment "Compare exports" 
bash ${scripts}/compare_exports.sh 



