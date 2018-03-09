#!/bin/bash

section=soma
export DATE=$(date +%Y%m%d_%H%M%S)
file_path=data/"$DATE"_nap_amplification_IO_DC_"$section"_stim_trunk.hdf5

for vrest in -65.0 -60.0
do
    export vrest
    for ttx in 0 1
    do
        export ttx
        python simulate_nap_EPSP_amplification_IO.py --ttx=$ttx --section=$section --output-file-path=$file_path \
            --vrest=$vrest
    done
done
