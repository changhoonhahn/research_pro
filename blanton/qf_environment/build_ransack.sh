#!/bin/usr/
Nransack=$2
Nrandom=$3
field_dir="/global/data/scr/chh327/primus/survey_regions/fields/"
ransack_dir="/global/data/scr/chh327/primus/data/ransack/"
idl -e ".r ransack_dat2fits.pro"

if [[ $1 == "primus" ]]; then
    ransack_fname="ransack_"
    random_fname="random_"
    fields=( "es1" "cosmos" "cfhtls_xmm" "cdfs" "xmm_swire" ) 
    ransack_command="ransack -r "$Nransack" "
    random_command="ransack -r "$Nrandom" "
    for i in $(seq 0 $((${#fields[@]}-1))); do 
        ransack_fname=$ransack_fname${fields[$i]}"_"
        random_fname=$random_fname${fields[$i]}"_"
#        polygon=$field_dir${fields[$i]}"_field_galex_window_2mask.ply"
        polygon=$field_dir"ransack_"${fields[$i]}".ply"
        ransack_command=$ransack_command$polygon" " 
        random_command=$random_command$polygon" " 
    done
    ransack_fname=$ransack_dir$ransack_fname$Nransack".dat"
    random_fname=$ransack_dir$random_fname$Nrandom".dat"
    ransack_command=$ransack_command$ransack_fname
    random_command=$random_command$random_fname
    echo $ransack_fname
    echo $ransack_command
    echo $random_fname
    echo $random_command
    if [[ ! -a $ransack_fname ]]; then  
        $ransack_command
        idl -e "ransack_dat2fits,'"$ransack_fname"'"
    fi
    if [[ ! -a $random_fname ]]; then  
        $random_command
        idl -e "ransack_dat2fits,'"$random_fname"'"
    fi 
elif [[ $1 == "sdss" ]]; then 
    field_dir='/global/data/scr/chh327/primus/science/mf/2165/'
    target_polygon=$field_dir"dr72bsafe0_galex_final.ply"
    edp_polygon=$field_dir"dr72bsafe0.ply"

    ransack_fname=$ransack_dir"ransack_sdss_edpfield_"$Nransack".dat"
    random_fname=$ransack_dir"random_sdss_targetfield_"$Nrandom".dat"
    ransack_command="ransack -r "$Nransack" "$edp_polygon" "$ransack_fname
    random_command="ransack -r "$Nrandom" "$target_polygon" "$random_fname
    echo $ransack_fname
    echo $ransack_command
    echo $random_fname
    echo $random_command
    $ransack_command
    idl -e "ransack_dat2fits,'"$ransack_fname"'"
    $random_command
    idl -e "ransack_dat2fits,'"$random_fname"'"
else 
    echo "survey must either be primus or sdss"
fi 

