import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rdf
import argparse
import datetime
import seaborn as sns
import pandas as pd

sns.set_style('darkgrid')

parser = argparse.ArgumentParser()
parser.add_argument("-psf", "--psf", type=str, help='NAMD topology PSF file')
#parser.add_argument("-dcd", "--dcd", type=str, help='NAMD DCD trajectory')
parser.add_argument('-w', "--workdir", type=str, help='Current working directory')
args = parser.parse_args()

u = mda.Universe(args.psf)
#print("Number of Frames: ", len(u.trajectory))

#u.trajectory[124]
#water_oxy = u.select_atoms("name OH2 and point 100 100 100 50")
#water_hyd = u.select_atoms("name H1 H2 and point 100 100 100 50")

water_oxy = u.select_atoms("name OH2")
water_hyd = u.select_atoms("name H1 H2")

print ("Coordinates loaded at ", datetime.datetime.now(),". RDF beginning now ...") 

rdf_oo = rdf.InterRDF(water_oxy, water_oxy, range=(0.0, 12.0), nbins=200, exclusion_block=(1,1))
print ('RDF OO starting at ', datetime.datetime.now())
rdf_oo.run()
print ('RDF OO finished at ', datetime.datetime.now())

rdf_oh = rdf.InterRDF(water_oxy, water_hyd, range=(0.0, 12.0), nbins=200, exclusion_block=(1,1))
print ('RDF OH starting at ', datetime.datetime.now())
rdf_oh.run()
print ('RDF OH finished at ', datetime.datetime.now())

rdf_hh = rdf.InterRDF(water_hyd, water_hyd, range=(0.0, 12.0), nbins=200, exclusion_block=(1,1))
print ('RDF HH starting at ', datetime.datetime.now())
rdf_hh.run()
print ('RDF HH finished at ', datetime.datetime.now())

#df = pd.DataFrame({'Bin_OO': rdf_oo.bins, 'RDF_OO': rdf_oo.rdf, 'Bin_OH': rdf_oh.bins, 'RDF_OH': rdf_oh.rdf, 'Bin_HH': rdf_hh.bins, 'RDF_HH': rdf_hh.rdf})
#df.to_csv(index=False)

#plt.axhline(1)
plt.plot(rdf_oo.bins, rdf_oo.rdf, label='O-O', color='mediumblue')
plt.plot(rdf_oh.bins, rdf_oh.rdf, label='O-H', color='forestgreen')
plt.plot(rdf_hh.bins, rdf_hh.rdf, label='H-H', color='lightseagreen')
plt.xlabel('Distance, $\AA$', fontsize=16)
plt.ylabel('Radial Distribution g(r)', fontsize=16)
plt.legend(fancybox=True, shadow=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(f'{args.workdir}/<filename>.png', bbox_inches='tight', dpi=1000)
