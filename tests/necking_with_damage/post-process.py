#!/usr/bin/env python

# Calculate the macroscopic force-displacement plot

# Undeformed cross-section area of the specimen
undeformed_area = 10

import os, fnmatch
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


mesh_name = 'Meshless.0.geo' # file containing the initial positions
d_name = 'Displacement' # base name for displacement files
f_name = 'ForceDensity' # base name for force density files
a_name = 'Area' # base name for nodal area files
datapath = './data/' # directory containing the data files
dirpath = './Output/' # directory containing the output files
plotsdir = './Plots/' # directory where the plots will be saved

f = open(dirpath+mesh_name, 'r')

# skip the first few lines
for i in range(5):
    f.readline()

# next line is the number of nodes
num_nodes = int(f.readline())

# read the initial position values
positions = np.zeros((num_nodes,3))

for idx in range(num_nodes):
    line = f.readline()
    positions[idx,0] = np.fromstring(line, dtype=float, sep=' ')[1] # x-coord
    positions[idx,1] = np.fromstring(line, dtype=float, sep=' ')[2] # y-coord
    positions[idx,2] = np.fromstring(line, dtype=float, sep=' ')[3] # z-coord


# find the number of time steps in the calculation
num_time_steps = len(fnmatch.filter(os.listdir(dirpath), f_name+'.*.res'))

# read data
displacement = np.zeros((num_time_steps,num_nodes,3))
forceDensity = np.zeros((num_time_steps,num_nodes,3))
nodal_area = np.zeros((num_time_steps,num_nodes))


# iterate over all time steps
for ts in range(num_time_steps):

    # read displacement values (vector)
    fname = dirpath+d_name+'.'+str(ts)+'.res'
    f = open(fname) # open file
    lines = f.readlines() # read lines

    for idx, line in enumerate(lines[1:]):
        displacement[ts,idx*2:(idx+1)*2,0] = np.fromstring(line, dtype=float, sep=' ')[0::3] # x-coord
        displacement[ts,idx*2:(idx+1)*2,1] = np.fromstring(line, dtype=float, sep=' ')[1::3] # y-coord
        displacement[ts,idx*2:(idx+1)*2,2] = np.fromstring(line, dtype=float, sep=' ')[2::3] # z-coord

    # read force density values (vector)
    fname = dirpath+f_name+'.'+str(ts)+'.res'
    f = open(fname) # open file
    lines = f.readlines() # read lines

    for idx, line in enumerate(lines[1:]):
        forceDensity[ts,idx*2:(idx+1)*2,0] = np.fromstring(line, dtype=float, sep=' ')[0::3] # x-coord
        forceDensity[ts,idx*2:(idx+1)*2,1] = np.fromstring(line, dtype=float, sep=' ')[1::3] # y-coord
        forceDensity[ts,idx*2:(idx+1)*2,2] = np.fromstring(line, dtype=float, sep=' ')[2::3] # z-coord

    # read nodal area values (scalar)
    fname = dirpath+a_name+'.'+str(ts)+'.res'
    f = open(fname) # open file
    lines = f.readlines() # read lines

    for idx, line in enumerate(lines[1:]):
        nodal_area[ts,idx*6:(idx+1)*6] = np.fromstring(line, dtype=float, sep=' ')

# Essential BCs were applied to top row of nodes
ymax = positions[:,1].max() # y-coordinate of the top row
top_nodes = positions[:,1] > ymax-1.0e-10

# Convert force density per unit area to force
nodal_force = np.einsum('ijk,ij->ijk', forceDensity, nodal_area)

# Compute the reaction force for top row of nodes
Force = - nodal_force[:,top_nodes,1].sum(axis=1) # y-coord

# Compute engineering stress
stress = Force/undeformed_area

Unorm = np.sqrt( np.einsum('ijk->i',displacement**2) / num_nodes)

PD_data = np.loadtxt(datapath+'PD_data.dat')
PD_Unorm = PD_data[:,0]
PD_Force = PD_data[:,1]
PD_stress = PD_Force / undeformed_area

IGA_data = np.loadtxt(datapath+'iga_data.dat')
IGA_Force = IGA_data[:,3]
IGA_Unorm = IGA_data[:,4]
IGA_stress = IGA_Force / undeformed_area

ambati_KL_data = np.loadtxt(datapath+'ambati_KL.dat')
ambati_KL_Unorm = ambati_KL_data[:,0]
ambati_KL_Force = ambati_KL_data[:,1]

ambati_solid_data = np.loadtxt(datapath+'ambati_solid.dat')
ambati_solid_Unorm = ambati_solid_data[:,0]
ambati_solid_Force = ambati_solid_data[:,1]

plt.figure(figsize=(10,5))

plt.plot(IGA_Unorm, IGA_Force*1e-3, 'r--', label='Alaydin et al., 2021 - IGA KL Shell') 
plt.plot(ambati_KL_Unorm, ambati_KL_Force, 'g--', label='Ambati et al., 2018 - IGA KL Shell') 
plt.plot(ambati_solid_Unorm, ambati_solid_Force, 'b--', label='Ambati et al., 2018 - IGA 3D Solid') 
plt.plot(PD_Unorm, PD_Force, 'k--', label='Peridynamic KL Shell without damage') 
plt.plot(Unorm, Force*1e-3, 'k-', lw=2, label='Peridynamic KL Shell with damage') 

plt.xlabel(r'${\rm U}_{\rm Norm}$ (mm)', fontsize=22)
plt.ylabel(r'Force (kN)', fontsize=22)

plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)

plt.legend(fontsize=14, loc=8)

plt.grid()
plt.xlim(0.0, 6.5)
plt.savefig(plotsdir+'necking_plate_load_displacement.png', format='png', dpi=200, bbox_inches='tight')
