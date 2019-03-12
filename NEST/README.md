# Network model and custom configuration files for NEST

In order to simulate the model named modelScript.py with NEST 2.14.0, the following custom set-ups of the simulator are required:

1. The user must install the module here provided named "nestmldev". Throughout this module, a custom neuron model is automatically installed in the simulator. The custom neuron is named "iaf_lifca", it is a leaky-integrate-and-fire neuron with with spike-frequency adaptation (SFA) due to calcium- and sodium-dependent after-hyperpolarization (AHP) currents (see arXiv:1804.03441v1 for details on the neuron model). The model has been obtained using the NESTML tool: starting from an high-level description of the neuron behavior with the file "iaf_lifca.nestml", the tool automatically generate all the files required to create and compile the specific module for NEST.

2. The user must recompile the NEST simulator kernel after having included the two files topology_parameter.h and topologymodule.cpp here provided in substitution of the same files located in the nest source folder "topology". The two files contains two specific connection kernel rules, not available in the standard NEST simulator, required to create the appropriate delay distribution and connection rules used in the model. Specifically, "exponential2" is an exponential distribution not depending on the distance between source and target neurons, used in the generation of synaptic delays, while exponential3 is an exponential kernel distribution depending on the distance with a rounding at the third decimal digit.

The files required by 1 and 2 are respectively stored in the folders "customNeuronModel" and "customConnectivityKernelRules".

Please, refer to official NEST simulator documentation for the instructions about how to install custom modules (point 1) and how to recompile the NEST simulator (point 2).
