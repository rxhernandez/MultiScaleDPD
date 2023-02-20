#!/usr/local/opt/python@3.8/bin/python3
import lib

# read structure pdb and MARTINI force field files
martini_pdb = open("step3_charmm2gmx.pdb").readlines()
protein_top = open("PROA_P.itp").read()
# initialize data
data = lib.Data(75, 75, 75)
# read MARTINI protein to data
data.add_molecule(martini_pdb[3:40], protein_top, "protein", "protein")
# fill water into the box
data.fill_water(ion_concentration=0.15)
#dump input
data.dump('protein.data')
