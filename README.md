# Scaling of a large-scale simulation of synchronous slow-wave and asynchronous awake-like activity of a cortical model with long-range interconnections

The repository contains code and models used to run large scale simulation of spiking Neural Networks supporting a range of dynamic states, from synchronous slow wave activity (SWA), characteristic of deep sleep or anesthesia, to fluctuating, asynchronous activity during wakefulness (AW).
Scaling analyses have been performed on the proprietary mixed time and event driven DPSNN simulation engine (Distributed Simulator of Plastic Spiking Neural Networks), for both SW and AW configurations, also varying connectivity parameters from shorter to longer range, as well as varying the mean firing rate of network activity. This implies the evaluation of simulation performance and robustness under different loads of local computation and communication. Distributed simulations have been also executed on the NEST platform, for comparison and validation.

The folder **DPSNN** contains the source files required to build the simulator and the DPSNN_script necessary to launch simulations.

The folder **NEST** contains:
	- the model of the neural network for a 48by48 grid of columns
	- the files required to build the module implementing the custom neuron model used for simulations
	- the files for specific connection kernel rules, not available in the standard NEST simulator

The folder **config** contains four set of configurations files that can be used to run simulations, both with DPSNN and NEST.
