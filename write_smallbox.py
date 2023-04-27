import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rdf
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-psf", "--psf", type=str, help='NAMD topology PSF file')
parser.add_argument("-dcd", "--dcd", type=str, help='NAMD trajectory')
parser.add_argument('-w', "--workdir", type=str, help='Current working directory')
args = parser.parse_args()

u = mda.Universe(args.psf, args.dcd)

u.trajectory[-1]
water_oxy = u.select_atoms("(prop abs x <= 40.0 and prop abs y <= 40.0 and prop abs z <= 40.0)", updating=True)
length = max(water_oxy.atoms.positions[:, 0]) - min(water_oxy.atoms.positions[:, 0])
print(length)
u.dimensions = [length, length, length, 90, 90, 90]



water_oxy.write(f"{args.workdir}<filename>.pdb")
