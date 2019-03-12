# Configuration files

The set of 3 choosen configuration files must be moved in the folder containing the DPSNN executable and/or the NEST python model script in order to obtain the correct simulation.
The 4 proposed configurations allow simulations of a grid of cortical columns, each composed of 1250 leaky integrate-and-fire neurons with spike frequency adaptation, and exponentially decaying inter-modular synaptic connectivity with varying spatial decay constant.
Parameters are tuned to produce networks showing SW or AW dynamical states.
In detail:
- AW2.8Hz_lambda04: AW state with shorter intermodular connectivity (lambda = 0.4), mean firing rate @2.8Hz
- AW8.8Hz_lambda04: AW state with shorter intermodular connectivity (lambda = 0.4), mean firing rate @8.8Hz
- SW_lambda04: SW state with shorter intermodular connectivity (lambda = 0.4), mean firing rate @3.1Hz
- SW_lambda06: SW state with longer intermodular connectivity (lambda = 0.6), mean firing rate @3.1Hz
