#!/bin/bash

export DATE=$(date +%Y%m%d_%H%M%S)
file_path=data/"$DATE"_nap_EPSC_amplification_soma.hdf5

for ttx in 0 1
do
    export ttx
    python simulate_nap_EPSC_amplification.py --ttx=$ttx --output-file-path=$file_path
    python simulate_nap_EPSC_amplification.py --ttx=$ttx --bar=1 --output-file-path=$file_path
    python simulate_nap_EPSC_amplification.py --ttx=$ttx --zd=1 --output-file-path=$file_path
done
