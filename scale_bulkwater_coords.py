import MDAnalysis as mda
from MDAnalysisTests.datafiles import PDB, DCD
import numpy as np
import os
import argparse

path = os.getcwd() + "/"

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--pdb", type=str)
parser.add_argument("-dcd", "--traj", help="trajectory file", default="full_27_out.dcd")
args = parser.parse_args()

u = mda.Universe(args.pdb, args.traj)
dimensions = (u.dimensions)
print('System dimensions are:', u.dimensions)

# Define central coordinate x y z #

minx = np.min(u.atoms.positions[:,0])
miny = np.min(u.atoms.positions[:,1])
minz = np.min(u.atoms.positions[:,2])
maxx = np.max(u.atoms.positions[:,0])
maxy = np.max(u.atoms.positions[:,1])
maxz = np.max(u.atoms.positions[:,2])

c_x = (minx + maxx) / 2
c_y = (miny + maxy) / 2
c_z = (minz + maxz) / 2

COM = u.atoms.center_of_mass()
print (COM)

############################################################
#                      ATOMGROUPS                          #
############################################################

protein = u.select_atoms("protein")
ions_all = u.select_atoms("name SOD CLA")
ions_int = ions_all.select_atoms(f'name SOD CLA and point {c_x} {c_y} {c_z} 65')
bulk_ions = ions_all - ions_int
all_water = u.select_atoms("resname TIP3")
all_water_ox = all_water.select_atoms("resname TIP3 and name OH2")
oxy_internal = u.select_atoms(f"name OH2 and point {c_x} {c_y} {c_z} 65")
oxy_resids = oxy_internal.residues.resids

selection_str = 'resname TIP3 and resid '
for resi in oxy_resids:
        selection_str += f'{resi} '

internal_waters = u.select_atoms(selection_str)
bulk_water = all_water - internal_waters
internal_protein_system = protein + internal_waters + ions_int #EDITED HERE - FORGOT INTERNAL IONS

#internal_protein_system.atoms.write('internal_protein_system.pdb')

bulk_water_ox = bulk_water.select_atoms("name OH2")
bulk_water_hyd = bulk_water.select_atoms("name H1 H2")


print('Central coordinate is by Calculation:', c_x, c_y, c_z)

############################################################
#                      SYSTEM DENSITY                      #
############################################################

volume_A=mda.lib.mdamath.box_volume(u.dimensions)
print('System volume is:', volume_A, 'cubic Å')
volume_cm = ((volume_A)*10**-24)
n_Av = (6.022*(10**23))
mr_h2o = (18.015)
mass_system = ((mr_h2o*(len(all_water_ox)))/(n_Av))

density=(mass_system/volume_cm)
print('Original system density is:', density, 'g/cm^3')

################### LDA ice density is approx. 0.92 g/cm^3. sf  variable should be density / 0.92 #####################
################### Because we're talking xyz dimensions, actually need the cube root of 0.92, which is 0.972 #########

sf=(density/0.972)
if sf < 1:
    raise ValueError("Invalid scale factor. Coordinates must be expanded not compressed.")
print ('System volume needs scaling by a scale factor of:', sf)

print ('Scaling beginning...')

#################################################################
#                      BULK WATER SCALING                       #
#################################################################

# Oxygen

bulk_water_ox_old = bulk_water_ox.positions
bulk_water_ox.positions = (bulk_water_ox.positions)*(sf)
bulk_water_ox_new = bulk_water_ox.positions

print ('Oxygen scaling complete')

# Hydrogen

change_ox = bulk_water_ox_new - bulk_water_ox_old
change_ox=(np.repeat(np.repeat(change_ox,2,axis=0),1,axis=1))
bulk_water_hyd.positions = bulk_water_hyd.positions+change_ox

print ('Hydrogen translation complete')

bulkwater_new = mda.Merge(bulk_water_ox.atoms, bulk_water_hyd.atoms)

# Ion Translation

bulk_ions.positions = bulk_ions.positions*sf

print ('Ions translation complete')

###############################################################
#            INTERNAL PROTEIN SYSTEM TRANSLATION              #
###############################################################

external_system_scaled = mda.Merge(bulkwater_new.atoms, bulk_ions.atoms)

# New Centre of Geometry

scaled_COM = external_system_scaled.atoms.center_of_mass()
COM_trans = tuple(map(lambda i, j: i - j, COM, scaled_COM))

external_system_translated = external_system_scaled.atoms.translate(COM_trans)
u_final = mda.Merge(external_system_translated.atoms, internal_protein_system.atoms)


###############################################################
#                SAVE AS PDB AND VISUALISE                    #
###############################################################

print ('Writing scaled PDB...')
u_final.atoms.write('6z6u_SCALED_FINAL_03.23_NPT.pdb')
print ('PDB write complete! NOW VISUALLY CHECK!')

quit()
