# Scaling of a large-scale simulation of synchronous slow-wave and asynchronous awake-like activity of a cortical model with long-range interconnections

Reference: Pastorelli E, Capone C, Simula F, Sanchez-Vives MV, Del Giudice P, Mattia M and Paolucci PS (2019) Scaling of a Large-Scale Simulation of Synchronous Slow-Wave and Asynchronous Awake-Like Activity of a Cortical Model With Long-Range Interconnections. Front. Syst. Neurosci. 13:33. doi: 10.3389/fnsys.2019.00033

Large-scale cortical model (several square centimeters of cortex) at biological neural and synaptic densities able to express both slow-waves and asynchronous awake states. The model is based on a bi-dimensional grid of neural populations, interconnected by long-range synaptic connectivity, which reflects the modular organization of the cortex.

Two implementations of the model are available. From one side, a NEST version of the model is distributed, in order to make it available to the large scientific community. At the same time, the model was also implemented using the proprietary mixed time and event driven DPSNN (Distributed Simulator of Plastic Spiking Neural Networks) simulator. In this natively distributed and parallel engine, the full neural system is represented by a network of C++ processes equipped with a communication infrastructure that can be easily interfaced with both Message Passing Interface (MPI) and custom software/hardware communication systems.

The simulations of different dynamical brain states constitute a challenge for parallel simulation of large multi-scale models of the brain. In DPSNN, efficient simulations resulted from the optimization of both the computational power and the inter-process communications. For this reason, speed-up measures as well as strong and weak scaling analysis was performed on the model simulated on DPSNN (Pastorelli et al, 2019).

We explored networks up to 192×192 modules, each composed of 1,250 integrate-and-fire neurons with spike-frequency adaptation, and exponentially decaying inter-modular synaptic connectivity with varying spatial decay constant. For the largest networks the total number of simulated synapses was over 70 billion.

# Content of the repository

The repository contains code and models used to run large scale simulation of spiking Neural Networks supporting a range of dynamic states, from synchronous slow wave activity (SWA), characteristic of deep sleep or anesthesia, to fluctuating, asynchronous activity during wakefulness (AW).
Scaling analyses have been performed on the proprietary mixed time and event driven DPSNN simulation engine (Distributed Simulator of Plastic Spiking Neural Networks), for both SW and AW configurations, also varying connectivity parameters from shorter to longer range, as well as varying the mean firing rate of network activity. This implies the evaluation of simulation performance and robustness under different loads of local computation and communication. Distributed simulations have been also executed on the NEST platform, for comparison and validation.

The folder **DPSNN** contains the source files required to build the simulator and the DPSNN_script necessary to launch simulations.

The folder **NEST** contains:
- the model of the neural network for a 48by48 grid of columns
- the files required to build the module implementing the custom neuron model used for simulations
- the files for specific connection kernel rules, not available in the standard NEST simulator

The folder **config** contains four set of configurations files that can be used to run simulations, both with DPSNN and NEST.
