#!/bin/usr/
setup idlutils; setup -r ~/primus/; setup -r ~/pro/
run=$1
Nransack=$2
Nrandom=$3

##########################################################################################
# Ransack, Random, and V_max,avail: 
##########################################################################################
#sh build_ransack.sh primus $Nransack $Nrandom
#idl -e "build_vmax_avail,"$run","$Nransack","$Nrandom",/primus"
##########################################################################################
# Build target and environment defining population:  
##########################################################################################
#idl -e "build_target_sample_mfdata_chh_test,"$run","$Nransack","$Nrandom",/literature,/primus"

#idl -e "build_environment_cylinder_mfdata_chh_test,"$run","$Nransack",/primus,/literature"

idl -e "build_smf_im_mf_vmax_mfdata_chh_test,"$run",/primus,/literature"
