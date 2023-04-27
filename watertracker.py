import MDAnalysis as mda
from MDAnalysisTests.datafiles import PDB, DCD
import numpy as np
import datetime

""" Script for tracking the internal waters from the inner hollow of the protein into the bulk solvent, used for identifying vaccuum space """


dcd = ['']
pdb = ''
u = mda.Universe(pdb,dcd)
print (u.dimensions)
allwater = u.select_atoms('name OH2') #This is the atomgroup of all the water oxygens in the system
COG = (allwater.center_of_geometry(pbc=True, compound='group', unwrap=False)) #System dependent - centre of geometry
print(COG)
print ("There are", len(allwater), 'waters')

H2O_in=[]
H2O_out=[]
frames=u.trajectory

#System dependent - assumes apoferritin is central in the system so the COG should be the hollow centre of 
#protein filled with water
print('Starting at', datetime.datetime.now())
for t in frames:
    waters_internal=u.select_atoms("name OH2 and point 100 100 100 35")
    waters_external=allwater.subtract(waters_internal)
    H2O_in.append(len(waters_internal))
    H2O_out.append(len(waters_external))

print('Finished tracking at', datetime.datetime.now())
H2O_in_array=np.array(H2O_in)
np.savetxt('internalwater.dat', H2O_in_array)
H2O_out_array=np.array(H2O_out)
np.savetxt('externalwater.dat', H2O_out_array)

import matplotlib.pyplot as plt
y1=H2O_in_array
y2=H2O_out_array
xmax=(len(H2O_in_array))
min_in=np.amin(H2O_in_array)
max_in=np.amax(H2O_in_array)
min_out=np.amin(H2O_out_array)
max_out=np.amax(H2O_out_array)

x=np.arange(0,xmax,1)

f, (ax2, ax) = plt.subplots(2, 1, sharex=True)
ax2.plot(y2, 'teal', label='298K')
ax2.legend(loc="lower left")
ax.plot(y1, 'teal')
ax.set_ylim((min_in-500), (max_in+500))
ax2.set_ylim((min_out-500), (max_out+500))
ax.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.xaxis.tick_top()
ax2.set_title('Water Transport')
plt.xlabel('Trajectory Frame No.')
plt.savefig('watertracker_full.png', dpi=300)
