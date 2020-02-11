
# DPSNN simulation engine

## Build
The building environment needs an MPI-compliant C/C++ compiler.
To build the executable, invoke 'make'.

## Simulate
Prepare in the same folder the DPSNN executable and the DPSNN configuration files for the specific model to be simulated.
To run simulations, invoke:

source DPSNN_script <CFX> <CFY> <MPI_procs> <hosts_list> <MPI_BTL> <time> LIFCA

with:
- CFX: Number of columns in X
- CFY: Number of columns in Y
- MPI_procs: Number of MPI processes
- hosts_list: Hosts list, separated by commas
- MPI_BTL: MPI BTL parameters
- time: Total simulation time in seconds

A grid of X by Y cortical columns, distributed on the specified number of MPI processes, will be simulated for the required number of seconds on the servers from the hosts list.


------------------------

# NEST simulations

## Network model and custom configuration files for NEST

In order to simulate the model named modelScript.py with NEST 2.14.0, the following custom set-ups of the simulator are required:

1. The user must install the module here provided named "nestmldev". Throughout this module, a custom neuron model is automatically installed in the simulator. The custom neuron is named "iaf_lifca", it is a leaky-integrate-and-fire neuron with with spike-frequency adaptation (SFA) due to calcium- and sodium-dependent after-hyperpolarization (AHP) currents (see arXiv:1804.03441v1 for details on the neuron model). The model has been obtained using the NESTML tool: starting from an high-level description of the neuron behavior with the file "iaf_lifca.nestml", the tool automatically generate all the files required to create and compile the specific module for NEST.

2. The user must recompile the NEST simulator kernel after having included the two files topology_parameter.h and topologymodule.cpp here provided in substitution of the same files located in the nest source folder "topology". The two files contains two specific connection kernel rules, not available in the standard NEST simulator, required to create the appropriate delay distribution and connection rules used in the model. Specifically, "exponential2" is an exponential distribution not depending on the distance between source and target neurons, used in the generation of synaptic delays, while exponential3 is an exponential kernel distribution depending on the distance with a rounding at the third decimal digit.

The files required by 1 and 2 are respectively stored in the folders "customNeuronModel" and "customConnectivityKernelRules".

Please, refer to official NEST simulator documentation for the instructions about how to install custom modules (point 1) and how to recompile the NEST simulator (point 2).


------------------------

# Models configuration files

The set of 3 choosen configuration files must be moved in the folder containing the DPSNN executable and/or the NEST python model script in order to obtain the correct simulation.
The 4 proposed configurations allow simulations of a grid of cortical columns, each composed of 1250 leaky integrate-and-fire neurons with spike frequency adaptation, and exponentially decaying inter-modular synaptic connectivity with varying spatial decay constant.
Parameters are tuned to produce networks showing SW or AW dynamical states.
In detail:
- AW2.8Hz_lambda04: AW state with shorter intermodular connectivity (lambda = 0.4), mean firing rate @2.8Hz
- AW8.8Hz_lambda04: AW state with shorter intermodular connectivity (lambda = 0.4), mean firing rate @8.8Hz
- SW_lambda04: SW state with shorter intermodular connectivity (lambda = 0.4), mean firing rate @3.1Hz
- SW_lambda06: SW state with longer intermodular connectivity (lambda = 0.6), mean firing rate @3.1Hz

