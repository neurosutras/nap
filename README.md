# nap
CA1 pyramidal neuron: Persistent Na current mediates steep synaptic amplification (Hsu et al., Neuron, 2018)

This repository contains the code used to simulate the full biophysically detailed CA1 pyramidal cell model described in the paper:
 
Hsu CL, Zhao X, Milstein AD, Spruston N. (2018) Persistent sodium current mediates the steep voltage dependence of spatial coding in 
hippocampal pyramidal neurons. Neuron 99, 1–16
 
Running these simulations requires:

1) NEURON - freely available from http://neuron.yale.edu

2) python2.7 - freely available from http://www.anaconda.com/download

3) open-source python modules: numpy, matplotlib, h5py, click, btmorph
 
For additional installation details, please refer to [install notes.txt](https://github.com/neurosutras/nap/blob/master/install%20notes.txt)
 
Usage:
 
Compile the .mod files (nrnivmodl on unix/linux/osx or mknrndll on mswin).

The following python scripts can be executed with input parameters via a command line interface, for example:
python simulate_nap_EPSC_amplification.py --section=soma --output-file-path=data/example_EPSC_amp.hdf5
python simulate_nap_EPSP_amplification.py --section=soma --output-file-path=data/example_EPSP_amp.hdf5
python simulate_nap_EPSP_amplification_IO.py --section=soma --vrest=-63 --output-file-path=data/example_EPSP_amp_IO.hdf5

or start an interactive ipython session and then run the scripts, for example:
run simulate_nap_EPSC_amplification --section=soma --output-file-path=data/example_EPSC_amp.hdf5


To export a series of simulations across a range of input parameters, execute the bash scripts:
./batch_nap_EPSC_amplification.sh
./batch_nap_EPSP_amplification.sh
./batch_nap_EPSP_amplification_IO.sh

This exports a set of simulation output data to an .hdf5 file. Then, a set of plot scripts can be executed to reproduce figures
from the paper, for example:

ipython
run plot_nap_EPSC_amplification  # generates plots similar to Figure S5
run plot_nap_EPSP_amplification  # generates plots similar to Figure 5C and S4 
run plot_nap_EPSP_amplification_IO  # generates plots similar to Figure 6B
