# icemodelling
Sample input scripts and analysis scripts for use alongside the models on Zenodo (https://doi.org/10.5281/zenodo.4415835)

Forcefields used: CHARMM 36 (toppar_c36_jul17)

Files included: 

  Input files:
  
      protein_input.namd - This is the namd run script with all the MD parameters included, as used for the apoferritin containing system, edit as required
      
      submit_namd.sh - This is an example submit script used for submission to ARCHER2 with the required resources to run MD of the apoferritin system
  
  Analysis files:
  
    coor2pdb.tcl - This is the tcl script used for converting the VMD binary output coor file into a pdb file. Intended for use through vmd without the GUI
    
    watertracker.py - This is the python script used for tracking the internal waters in the hollow centre of the apoferritin molecule. MDAnalysis required
    
    Energy_Time.ipynb - This jupyter notebook reads in a text file of total energies and writes the plot 
    
    Density_Time.ipynb - This jupyter notebook reads in a text file of system volumes, converts to density, and writes the plot
    
    write_smallbox.py - This is the python script for selecting a smaller cube of water from the final frame of the larger 645 Angstrom trajectory and writing it to pdb format
    
    rdf.py - This python script plots the RDF for water using the pdb written out from write_smallbox.py
    
