#!/bin/bash

export DATE=$(date +%Y%m%d_%H%M%S)
section="$1"
file_path=data/"$DATE"_nap_amplification_DC_"$section"_stim_trunk.hdf5

for ttx in 0 1
do
    export ttx
    for ap5 in 0 1
    do
        export ap5
        python simulate_nap_EPSP_amplification.py --ap5=$ap5 --ttx=$ttx --section=$section --output-file-path=$file_path
    done
done
