#!/bin/usr/
setup idlutils; setup -r ~/primus/; setup -r ~/pro/ 

for run in 3 3 4 7 8; do 
    echo "run "$run
    idl -e "build_smf_im_mf_vmax,"$run",/twoenv"
done
