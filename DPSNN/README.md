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
