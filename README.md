# MultiScaleDPD

This repository contains the source code accompanying the paper: "[Construction of Multiscale Dissipative Particle Dynamics (DPD) Models from Other Coarse-Grained Models](https://pubs.acs.org/doi/full/10.1021/acsomega.4c01868)".

And the (PENDING) paper: "Coarse-Grained Dissipative Particle Dynamics (DPD) Simulation of Cytochrome \textit{c} Facilitated Binding Between Lipid Bilayers and Citrate-Coated Gold Nanoparticles"


## Overview

MultiScaleDPD facilitates the construction of multiscale DPD models derived from existing coarse-grained models, enabling simulations in LAMMPS.

## Features

- **Model Conversion**: Transform MARTINI models into DPD-compatible formats.
- **LAMMPS Integration**: Generate input files for seamless execution in LAMMPS.
- **Customization**: Support for advanced configurations, including the incorporation of additional molecular components like heme groups.

## Getting Started

### Prerequisites

- **Software**:
	- [CHARMM-GUI](www.charmm-gui.org)
	- [LAMMPS](https://github.com/lammps/lammps)
	- NumPy

### Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/rxhernandez/MultiScaleDPD.git
   ```



<hr>

Documentation
----------------

Process of building a DPD protein to run in LAMMPS:
1. build a MARTINI model of your protein using 
CHARMM-gui (https://www.charmm-gui.org/?doc=input/martini.solution) with martini22 force field. You need two files from the result, charmm-gui-xxxxxxxxx/gromacs/step3_charmm2gmx.pdb is the structure file, and charmm-gui-xxxxxxxxx/gromacs/PROA_P.itp is the MARTINI force field file.
2. Put these two files under the same dir as lib.py and convert.py
3. Open convert.py and step3_charmm2gmx.pdb with a text editor. Modify the line numbers in 
convert.py according to the lines in step3_charmm2gmx.pdb
4. run convert.py (numpy is required)
5. equilibration.lmps and production.lmps are the input files for LAMMPS. Run equilibration first, 
then run production.
6. Advanced useage example avaiable in example2. We add the heme to the pdb generated by CHARMM-gui 
(step3_charmm2gmx.pdb) to generate martini.pdb. Then adding customize bonds and angles that connecting heme and protein backbone as defined in convert.py.
*Note: make sure the LAMMPS has bug fixes for July 2022

<hr>

Authors
----------------

The MultiScaleDPD codes and databaess were developed by Yinhan Wang

Contributors can be found [here](https://github.com/rxhernandez/MultiScaleDPD/graphs/contributors).

<hr>

Citing
----------------
If you use database or codes, please cite the paper:

> Wang, Y; Hernandez R.; _ACS Omega_ **9**, 17667 (2024). [https://doi.org/10.1021/acsomega.4c01868]

Acknowledgment
----------------

This work was supported by 
the National Science Foundation
under Grant No. CHE-2001611, the NSF 884 Center for Sustainable Nanotechnology (CSN). The CSN is 885 part of the Centers for Chemical Innovation Program. 

<hr>

License
----------------

MultiScaleDPD code and databases are distributed under terms of the [MIT License](https://github.com/rxhernandez/MultiScaleDPD/blob/main/LICENSE).
