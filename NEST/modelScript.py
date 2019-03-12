import nest
import nest.topology as tp
import numpy as np
import time
import socket
import matplotlib.pyplot as plt
import pylab
import sys
import math

nest.Install("nestmldev")

loc_h = nest.Rank()
if loc_h == 0:
    print ('MPI process %x on host' % nest.Rank(),socket.gethostname())
buildStart = time.time()
resolution = .1
resizeFactor = 1
debug = 0

#! Simulator Initialization
#! ==============================

nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads':6})
nest.SetKernelStatus({'overwrite_files': True})
nest.SetKernelStatus({'resolution': resolution, 'print_time': True})

msd = 1457427778               #int(time.time())              # Master seed...
N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]   # Number of virtual processes...
pyrngs = [np.random.RandomState(s) for s in range(msd, msd+N_vp)]  # Seeding the Python RNGs
nest.SetKernelStatus({'grng_seed' : msd+N_vp})                        # Seeding the global RNG
nest.SetKernelStatus({'rng_seeds' : range(msd+N_vp+1, msd+2*N_vp+1)}) # Seeding the per-process RNGs
totMPIproc = nest.NumProcesses()

if loc_h == 0:
    print ('Reading input files and setting parameters')
#! Setting Parameters
#! ==============================

# Network numbers
N_rows = 48
N_columns = 48
stencil = 4
ExtNeu = 400
ExtNeuF = float(ExtNeu)
simDuration = 10000.0
totConn = 1.5120
totConnPerMod = 0.9 
lambdaConn = 0.4

# Neuronal dynamic parameters
lambdaD = [0.037,0.037,0.14]

#! Read total number of pop from input column file
#! ==============================
with open('column.txt', 'r') as inputFile:
    for line in inputFile:
        if len(line.rstrip()) != 0:
            if line[0] != '#':
                B = line.split()
                N_pop_tot = (int(B[0])+1)
inputFile.close()

ID_pop = np.empty((N_pop_tot),dtype=np.int)
N_pop = np.empty((N_pop_tot),dtype=np.int)
JExt = np.empty((N_pop_tot))
DJExt = np.empty((N_pop_tot))
NuExt = np.empty((N_pop_tot))
Tau = np.empty((N_pop_tot))
Theta = np.empty((N_pop_tot))
H = np.empty((N_pop_tot))
Tarp = np.empty((N_pop_tot))
AlphaC = np.empty((N_pop_tot))
TauC = np.empty((N_pop_tot))
gC = np.empty((N_pop_tot))

sigmaScaleFactor = math.sqrt(ExtNeuF/400.)
if loc_h == 0:
    print ('sigmaScaleFactor '), sigmaScaleFactor

#! Read input column file
#! ==============================
with open('column.txt', 'r') as inputFile:
    for line in inputFile:
        if len(line.rstrip()) != 0:
            if line[0] != '#':
                B = line.split()
                item = (int(B[0]))
                ID_pop[item] = (int(B[0]))
                N_pop[item] = (int(B[1]))
                JExt[item] = (float(B[2])/resolution)
                DJExt[item] = (float(B[3])*sigmaScaleFactor)
                NuExt[item] = (float(B[4]))
                Tau[item] = (float(B[5]))
                Theta[item] = (float(B[6]))
                H[item] = (float(B[7]))
                Tarp[item] = (float(B[8]))
                AlphaC[item] = (float(B[9]))
                TauC[item] = (float(B[10]))
                gC[item] = (float(B[11]))
inputFile.close()

N_exf = N_pop[0]
N_exb = N_pop[1]
N_inh = N_pop[2]
N_tot = N_pop[0] +  N_pop[1] + N_pop[2]
ID_exf = ID_pop[0]
ID_exb = ID_pop[1]
ID_inh = ID_pop[2]

if debug:
    if loc_h == 0:
        print ('Total number of populations: '), N_pop_tot
        print ('Number of neurons per column: '), N_pop
        print ('Total number of neurons: '), N_tot
        print ('JExt: '), JExt
        print ('DJExt: '), DJExt
        print ('NuExt: '), NuExt
        print ('Tau: '), Tau
        print ('Theta: '), Theta
        print ('H: '), H
        print ('Tarp: '), Tarp
        print ('AlphaC: '), AlphaC
        print ('TauC: '), TauC
        print ('gC: '), gC

#! Read input bath file
#! ==============================
NoiseMtrx  = np.empty((stencil,stencil),dtype=np.int)
y = 0
with open('bathEfficacy.txt', 'r') as inputFile:
    for line in inputFile:
        if len(line.rstrip()) != 0:
            if line[0] != '#':
                B = line.split()
                for x in range (0,stencil):
                    NoiseMtrx[x][y] = (int(B[x]))
        y = y + 1

inputFile.close()

if debug:
    if loc_h == 0:
        print ("NoiseMtrx:")
        print (NoiseMtrx)

#! Read input connectivity file
#! ==============================
ConnectivityMatrix = np.empty((N_pop_tot,N_pop_tot,stencil,stencil))
J = np.empty((N_pop_tot,N_pop_tot))
DJ = np.empty((N_pop_tot,N_pop_tot))
Dmin = np.empty((N_pop_tot,N_pop_tot))
Dmax = np.empty((N_pop_tot,N_pop_tot))
with open('connectivity.txt', 'r') as inputFile:
    for line in inputFile:
        if len(line.rstrip()) != 0:
            if line[0] != '#':
                B = line.split()
                sourcePop = int(B[0])
                targetPop = int(B[1])
                J[sourcePop][targetPop] = float(B[2])/resolution
                DJ[sourcePop][targetPop] = (float(B[3]))
                Dmin[sourcePop][targetPop] = float(B[4])+1.
                Dmax[sourcePop][targetPop] = float(B[5])+1.
                for y in range(0,stencil):
                    for x in range(0,stencil):
                        index = 6 + y*stencil + x
                        ConnectivityMatrix[sourcePop][targetPop][y][x] = float(B[index])/resizeFactor
inputFile.close()

if debug:
    if loc_h == 0:
        print ("J:")
        print (J)
        print ("DJ:")
        print (DJ)
        print ("Dmin:")
        print (Dmin)
        print ("Dmax:")
        print (Dmax)
        print ("ConnectivityMatrix:")
        print (ConnectivityMatrix)

elementCreateStart = time.time()

#! Building neuron models
#! ==============================
# iaf_lifca neuron model parameters:
# V_m        double - Membrane potential in mV
# w          double - Adaptative state variable
# V_th       double - Spike threshold in mV
# V_reset    double - Reset potential of the membrane in mV
# E_L        double - Resting membrane potential in mV
# t_ref      double - Duration of refractory period in ms
# C_m        double - Capacitance of the membrane in pF
# Tau        double - Membrane time constant in ms
# tau_w      double - Adaptative time constant in ms
# q_w        double - Adaptative conductance in nS
# I_e        double - Constant input current in pA

if loc_h == 0:
    print ('Building neuron models')

exc_dict = {"Tau": Tau[ID_exf],"V_th": Theta[ID_exf],"V_reset": H[ID_exf],"t_ref": Tarp[ID_exf],"AlphaC": AlphaC[ID_exf],"tau_w": TauC[ID_exf],"g_w": gC[ID_exf]}
nest.CopyModel("iaf_lifca", "exf_iaf_neuron", params=exc_dict)
exc_dict = {"Tau": Tau[ID_exb],"V_th": Theta[ID_exb],"V_reset": H[ID_exb],"t_ref": Tarp[ID_exb],"AlphaC": AlphaC[ID_exb],"tau_w": TauC[ID_exb],"g_w": gC[ID_exb]}
nest.CopyModel("iaf_lifca", "exb_iaf_neuron", params=exc_dict)
inh_dict = {"Tau": Tau[ID_inh],"V_th": Theta[ID_inh],"V_reset": H[ID_inh],"t_ref": Tarp[ID_inh],"AlphaC": AlphaC[ID_inh],"tau_w": TauC[ID_inh],"g_w": gC[ID_inh]}
nest.CopyModel("iaf_lifca", "inh_iaf_neuron", params=inh_dict)

#! Poissonian noise
#! ==============================
if loc_h == 0:
    print ('Building Poissonian Noise rate generator')
def rateInit(pos, localNuExt):
    cfx = int(pos[0])
    cfy = int(pos[1])
    cfx = -1 * abs(cfx) + N_columns/2
    cfy = -1 * abs(cfy) + N_rows/2
    cfx = int(min(cfx,NoiseStencil_x))
    cfy = int(min(cfy,NoiseStencil_y))
    Cext = NoiseMtrx[cfx-1][cfy-1]
    return Cext*localNuExt
NoiseStencil_x = stencil
NoiseStencil_y = stencil

rate_exf = (400./ExtNeuF) * NuExt[0]
rate_exb = (400./ExtNeuF) * NuExt[1]
rate_inh = (400./ExtNeuF) * NuExt[2]
nest.SetDefaults('poisson_generator', {'rate': rate_exf})
nest.CopyModel('poisson_generator','noise_exf')
nest.SetDefaults('poisson_generator', {'rate': rate_exb})
nest.CopyModel('poisson_generator','noise_exb')
nest.SetDefaults('poisson_generator', {'rate': rate_inh})
nest.CopyModel('poisson_generator','noise_inh')

#! Spikes detector
#! ==============================
if loc_h == 0:
    print ('Building probes')
nest.SetDefaults('spike_detector', {'withtime': True,'withgid' : True,'to_file' : False})
nest.CopyModel('spike_detector','spikeTotal',{'to_file' : True})

elementCreateEnd = time.time()

#! Creating layers
#! ==============================
if loc_h == 0:
    print ('Building layers')

layerCreateStart = time.time()
layerBase = {'rows'     : N_rows, 
             'columns'  : N_columns,
             'extent'   : [(N_columns*1.0), (N_rows*1.0)]}

layerBase.update({'elements': ['exf_iaf_neuron',N_exf,
                               'exb_iaf_neuron',N_exb,
                               'inh_iaf_neuron',N_inh]})
ly = tp.CreateLayer(layerBase)

lyNoiseExf = tp.CreateLayer({'rows': 1, 'columns': 1, 
                               'elements': ['noise_exf',ExtNeu]})
lyNoiseExb = tp.CreateLayer({'rows': 1, 'columns': 1, 
                               'elements': ['noise_exb',ExtNeu]})
lyNoiseInh = tp.CreateLayer({'rows': 1, 'columns': 1, 
                               'elements': ['noise_inh',ExtNeu]})

lySpikeTotal = tp.CreateLayer({'rows': 1, 'columns': 1, 
                               'elements': ['spikeTotal']})

layerCreateEnd = time.time()

if loc_h == 0:
    print ('Setting rates of Poissonian noise')

#! Building connections
#! ==============================
if loc_h == 0:
    print ('Building connections')

networkConn = []

nest.CopyModel('static_synapse','dump')
localBaseExc = {'connection_type': 'convergent',
                'mask': {'circular': {'radius': 3.5}},
                'kernel': {'exponential3': {'a': totConnPerMod/totConn, 'tau': lambdaConn}},
                'synapse_model': 'dump',
                'allow_autapses': True,
                'allow_multapses': False}

for conn in[{'sources': {'model': 'exf_iaf_neuron'}, 
             'targets': {'model': 'exf_iaf_neuron'},
             'weights': {'normal':{'mean': J[0][0],'sigma': DJ[0][0]*J[0][0],'min': 0.0,'max': 2.0*J[0][0]}},
             'delays' : {'exponential2':{'tau': lambdaD[0],'min': Dmin[0][0],'max': Dmax[0][0]}}},
            {'sources': {'model': 'exf_iaf_neuron'}, 
             'targets': {'model': 'exb_iaf_neuron'},
             'weights': {'normal':{'mean': J[0][1],'sigma': DJ[0][1]*J[0][1],'min': 0.0,'max': 2.0*J[0][1]}},
             'delays' : {'exponential2':{'tau': lambdaD[0],'min': Dmin[0][1],'max': Dmax[0][1]}}},
            {'sources': {'model': 'exf_iaf_neuron'}, 
             'targets': {'model': 'inh_iaf_neuron'},
             'weights': {'normal':{'mean': J[0][2],'sigma': DJ[0][2]*J[0][2],'min': 0.0,'max': 2.0*J[0][2]}},
             'delays' : {'exponential2':{'tau': lambdaD[0],'min': Dmin[0][2],'max': Dmax[0][2]}}},
            {'sources': {'model': 'exb_iaf_neuron'}, 
             'targets': {'model': 'exf_iaf_neuron'},
             'weights': {'normal':{'mean': J[1][0],'sigma': DJ[1][0]*J[1][0],'min': 0.0,'max': 2.0*J[1][0]}},
             'delays' : {'exponential2':{'tau': lambdaD[1],'min': Dmin[1][0],'max': Dmax[1][0]}}},
            {'sources': {'model': 'exb_iaf_neuron'}, 
             'targets': {'model': 'exb_iaf_neuron'},
             'weights': {'normal':{'mean': J[1][1],'sigma': DJ[1][1]*J[1][1],'min': 0.0,'max': 2.0*J[1][1]}},
             'delays' : {'exponential2':{'tau': lambdaD[1],'min': Dmin[1][1],'max': Dmax[1][1]}}},
            {'sources': {'model': 'exb_iaf_neuron'}, 
             'targets': {'model': 'inh_iaf_neuron'},
             'weights': {'normal':{'mean': J[1][2],'sigma': DJ[1][2]*J[1][2],'min': 0.0,'max': 2.0*J[1][2]}},
             'delays' : {'exponential2':{'tau': lambdaD[1],'min': Dmin[1][2],'max': Dmax[1][2]}}}]:
    ndict = localBaseExc.copy()
    cdict = ndict.update(conn)
    networkConn.append(ndict)

localBaseInh = {'connection_type': 'convergent',
                'mask': {'circular': {'radius': .5}},
                'kernel': {'exponential3': {'a': totConnPerMod, 'tau': lambdaConn}},
                'synapse_model': 'dump',
                'allow_autapses': True,
                'allow_multapses': False}

for conn in[{'sources': {'model': 'inh_iaf_neuron'}, 
             'targets': {'model': 'exf_iaf_neuron'},
             'weights': {'normal':{'mean': J[2][0],'sigma': -DJ[2][0]*J[2][0],'min': 2.0*J[2][0],'max': 0.0}},
             'delays' : {'exponential2':{'tau': lambdaD[2],'min': Dmin[2][0],'max': Dmax[2][0]}}},
            {'sources': {'model': 'inh_iaf_neuron'}, 
             'targets': {'model': 'exb_iaf_neuron'},
             'weights': {'normal':{'mean': J[2][1],'sigma': -DJ[2][1]*J[2][1],'min': 2.0*J[2][1],'max': 0.0}},
             'delays' : {'exponential2':{'tau': lambdaD[2],'min': Dmin[2][1],'max': Dmax[2][1]}}},
            {'sources': {'model': 'inh_iaf_neuron'}, 
             'targets': {'model': 'inh_iaf_neuron'},
             'weights': {'normal':{'mean': J[2][2],'sigma': -DJ[2][2]*J[2][2],'min': 2.0*J[2][2],'max': 0.0}},
             'delays' : {'exponential2':{'tau': lambdaD[2],'min': Dmin[2][2],'max': Dmax[2][2]}}}]:
    ndict = localBaseInh.copy()
    cdict = ndict.update(conn)
    networkConn.append(ndict)

#! Connecting Network
#! ==============================
if loc_h == 0:
    print ('Connecting Network')
recurrentStart = time.time()
[tp.ConnectLayers(ly,ly,conn) for conn in networkConn]
recurrentEnd = time.time()

#! Connecting Poissonian noise
#! ==============================
if loc_h == 0:
    print ('Connecting Poissonian noise')

noiseBase = {'connection_type': 'divergent',
             'allow_autapses': False,
             'allow_multapses': False}

exfNoise = {'sources': {'model': 'noise_exf'}, 
            'targets': {'model': 'exf_iaf_neuron'},
            'weights': {'normal':{'mean': JExt[0],'sigma': DJExt[0]*JExt[0],
                                  'min': 0.0,'max': 2.0*JExt[0]}}}
exbNoise = {'sources': {'model': 'noise_exb'}, 
            'targets': {'model': 'exb_iaf_neuron'},
            'weights': {'normal':{'mean': JExt[1],'sigma': DJExt[1]*JExt[1],
                                  'min': 0.0,'max': 2.0*JExt[1]}}}

inhNoise = {'sources': {'model': 'noise_inh'}, 
            'targets': {'model': 'inh_iaf_neuron'},
            'weights': {'normal':{'mean': JExt[2],'sigma': DJExt[2]*JExt[2],
                                  'min': 0.0,'max': 2.0*JExt[2]}}}

poissonStart = time.time()
ndict = noiseBase.copy()
cdict_exfNoise = ndict.update(exfNoise)
tp.ConnectLayers(lyNoiseExf,ly,ndict)
cdict_exfNoise = ndict.update(exbNoise)
tp.ConnectLayers(lyNoiseExb,ly,ndict)
cdict_exfNoise = ndict.update(inhNoise)
tp.ConnectLayers(lyNoiseInh,ly,ndict)
poissonEnd = time.time()

#! Connecting probes
#! ==============================
if loc_h == 0:
    print ("Connecting probes")

probeStart = time.time()
cdict_spikeTotal_exf = {'connection_type': 'convergent',
                        'sources': {'model': 'exf_iaf_neuron'},
                        'targets': {'model': 'spikeTotal'},
                        'allow_autapses': False,
                        'allow_multapses': False}
tp.ConnectLayers(ly,lySpikeTotal,cdict_spikeTotal_exf)

cdict_spikeTotal_exb = {'connection_type': 'convergent',
                        'sources': {'model': 'exb_iaf_neuron'},
                        'targets': {'model': 'spikeTotal'},
                        'allow_autapses': False,
                        'allow_multapses': False}
tp.ConnectLayers(ly,lySpikeTotal,cdict_spikeTotal_exb)

cdict_spikeTotal_inh = {'connection_type': 'convergent',
                        'sources': {'model': 'inh_iaf_neuron'},
                        'targets': {'model': 'spikeTotal'},
                        'allow_autapses': False,
                        'allow_multapses': False}
tp.ConnectLayers(ly,lySpikeTotal,cdict_spikeTotal_inh)
probeEnd = time.time()
buildEnd = time.time()

#! Simulation
#! ==============================
if loc_h == 0:
    print ('End of building phase')
    print ('Building time:                      %.2f s') % (buildEnd-buildStart)
    print ('Elements creation:                  %.2f s') % (elementCreateEnd-elementCreateStart)
    print ('Layers creation:                    %.2f s') % (layerCreateEnd-layerCreateStart)
    print ('Building recurrent connection time: %.2f s') % (recurrentEnd-recurrentStart)
    print ('Building poisson connection time:   %.2f s') % (poissonEnd-poissonStart)
    print ('Building probe connection time:     %.2f s') % (probeEnd-probeStart)
    print ('sigmaScaleFactor:                   %.7f  ') % sigmaScaleFactor
    print ('Simulating')
startSimulate = time.time()
nest.Simulate(simDuration)
endSimulate = time.time()

#! Print information
#! ==============================
elapsed = endSimulate - startSimulate
if loc_h == 0:
    if elapsed > 0.0:
        print ('Simulation time: %f sec' % (elapsed))
        print ('Realtime factor: %1.4f' % (float(simDuration/1000.)/elapsed))
    
kernelStatus = nest.GetKernelStatus()

if loc_h == 0:
    print (' ')
    print ('Simulator Kernel Status:')
    print ('Number of MPI processes = '), kernelStatus['num_processes'] #equal to: totMPIproc
    print ('Local number of threads = '), kernelStatus['local_num_threads']
    print ('Number of virtual processes = '), kernelStatus['total_num_virtual_procs'] #equal to: N_vp
    print ('MPI Process rank = '), loc_h
    print (' ')

    print ('Neurons / column  :'), N_tot
    print ('       ExitatoryF :'), N_pop[0]
    print ('       ExitatoryB :'), N_pop[1]
    print ('       Inhibitory :'), N_pop[2]
    print ('Virtual Processes :'), N_vp
    print ('SimTime Resolution:'), kernelStatus['resolution']
    print ('Master Seed       :'), msd
    print ('Simulated time    : %.2f s') % (simDuration/1000.)
    print ('Building time:                      %.2f s') % (buildEnd-buildStart)
    print ('Elements creation:                  %.2f s') % (elementCreateEnd-elementCreateStart)
    print ('Layers creation:                    %.2f s') % (layerCreateEnd-layerCreateStart)
    print ('Building recurrent connection time: %.2f s') % (recurrentEnd-recurrentStart)
    print ('Building poisson connection time:   %.2f s') % (poissonEnd-poissonStart)
    print ('Building probe connection time:     %.2f s') % (probeEnd-probeStart)
    print ('Simulation time:                    %.2f s') % (endSimulate-startSimulate)


if loc_h == 0:
    orig_stdout = sys.stdout
    f = open('logNEST_%sx%s_H%s_T%s.log' % (N_rows,N_columns,kernelStatus['num_processes'],kernelStatus['local_num_threads']), 'w')
    sys.stdout = f
    print ('Simulation log for a grid problem size %dx%d') % (N_rows,N_columns)
    print ('execution date: ') + (time.strftime("%d/%m/%Y"))
    print (' ')
    print ('Simulator Kernel Status:')
    print ('Number of MPI processes = '), kernelStatus['num_processes'] #equal to: totMPIproc
    print ('Local number of threads = '), kernelStatus['local_num_threads']
    print ('Number of virtual processes = '), kernelStatus['total_num_virtual_procs'] #equal to: N_vp
    print ('MPI Process rank = '), loc_h
    print (' ')

    print ('Neurons / column  :'), N_tot
    print ('       ExitatoryF :'), N_pop[0]
    print ('       ExitatoryB :'), N_pop[1]
    print ('       Inhibitory :'), N_pop[2]
    print ('Virtual Processes :'), N_vp
    print ('SimTime Resolution:'), kernelStatus['resolution']
    print ('Master Seed       :'), msd
    print ('Simulated time    : %.2f s') % (simDuration/1000.)
    print ('Building time:                      %.2f s') % (buildEnd-buildStart)
    print ('Elements creation:                  %.2f s') % (elementCreateEnd-elementCreateStart)
    print ('Layers creation:                    %.2f s') % (layerCreateEnd-layerCreateStart)
    print ('Building recurrent connection time: %.2f s') % (recurrentEnd-recurrentStart)
    print ('Building poisson connection time:   %.2f s') % (poissonEnd-poissonStart)
    print ('Building probe connection time:     %.2f s') % (probeEnd-probeStart)
    print ('Simulation time:                    %.2f s') % (endSimulate-startSimulate)
    sys.stdout = orig_stdout
    f.close()

    fpar = open('param.txt', 'w')
    sys.stdout = fpar
    print ('# Parameters for the NEST utilities')
    print ('# the file can be used with spikeConverter and rateGenerator')
    print ('')
    print ('  CFX = %d	  # Grid size in x direction') % (N_rows)
    print ('  CFY = %d	  # Grid size in y direction') % (N_columns)
    print ('  H = %d	  # Number of processes') %  kernelStatus['total_num_virtual_procs']
    print ('  pop0 = %d	  # Neurons in pop0') % (N_pop[0])
    print ('  pop1 = %d	  # Neurons in pop1') % (N_pop[1])
    print ('  pop2 = %d	  # Neurons in pop2') % (N_pop[2])
    print ('  offset = 2	  # Offset of the first NEST neuron ID') 
    print ('  bin = 5  	  # Binning in ms, used to generate rates from spikes')
    sys.stdout = orig_stdout
    fpar.close()

plt.show()
