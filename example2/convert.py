#!/usr/local/opt/python@3.8/bin/python3
import lib

# read structure pdb and MARTINI force field files
martini_pdb = open("martini.pdb").readlines()
protein_top = open("PROA_P.itp").read()
heme_top = open("heme-CG_corrected.itp").read()

# initialize data
data = lib.Data(75, 75, 75)
# read MARTINI protein to data
data.add_molecule(martini_pdb[:236], protein_top, "protein", "protein")
# read MARTINI heme to data
data.add_molecule(martini_pdb[236:255], heme_top, "protein", "heme")
# add bond that connect heme to backbone
data.add_bond(40, 255, lib.k_constraints, 0.239)
data.add_bond(182, 255, lib.k_constraints, 0.239)
# angles for bond connecting heme and backbone
data.add_angle(40, 255, 237, 500, 90)
data.add_angle(40, 255, 239, 500, 90)
data.add_angle(40, 255, 241, 500, 90)
data.add_angle(40, 255, 243, 500, 90)
data.add_angle(182, 255, 237, 500, 90)
data.add_angle(182, 255, 239, 500, 90)
data.add_angle(182, 255, 241, 500, 90)
data.add_angle(182, 255, 243, 500, 90)
# fill water into the box
data.fill_water(ion_concentration=0.15)

data.dump('protein.data', "neigh_modify exclude group heme heme\n")
