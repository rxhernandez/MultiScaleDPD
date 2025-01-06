# MultiScaleDPD

This repository contains the source code accompanying the paper: "[Construction of Multiscale Dissipative Particle Dynamics (DPD) Models from Other Coarse-Grained Models](https://pubs.acs.org/doi/full/10.1021/acsomega.4c01868)".

And the (PENDING) paper: "Coarse-Grained Dissipative Particle Dynamics (DPD) Simulation of Cytochrome \textit{c} Facilitated Binding Between Lipid Bilayers and Citrate-Coated Gold Nanoparticles"


## Overview

MultiScaleDPD facilitates the construction of multiscale DPD models derived from existing coarse-grained models, enabling simulations in LAMMPS.

## Features

- **Model Conversion**: Transform MARTINI models into DPD-compatible formats.
- **LAMMPS Integration**: Generate input files for seamless execution in LAMMPS.
- **Customization**: Support for advanced configurations, including the incorporation of additional molecular components like heme groups.
- **Improved DPD Interactions**: Incorporates enhanced pair interaction potentials for DPD simulations using `pair_dpd.cpp`.

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
2. **Navigate to the Directory**:
   ```bash
   cd MultiScaleDPD
   ```

3. **Install Dependencies**:
   ```bash
   pip install numpy
   ```

### Building LAMMPS

To use the improved `pair_dpd.cpp` file, you need to build LAMMPS with the appropriate modifications:

1. **Obtain LAMMPS Source Code**:
	- Clone the LAMMPS repository:
   ```bash
   git clone https://github.com/lammps/lammps.git
   cd lammps
   ```

2. **Replace `pair_dpd.cpp`**:
	- Copy the provided `pair_dpd.cpp` file into the `src/DPD-BASIC` directory of the LAMMPS source code:
   ```bash
   cp /path/to/MultiScaleDPD/pair_dpd.cpp src/DPD-BASIC/
   ```

3. **Build LAMMPS**:
	- Follow the official instructions to compile LAMMPS with your desired settings:
	[Build LAMMPS](https://docs.lammps.org/Build.html)

### Usage

1. **Build a MARTINI Model**:
   - Utilize [CHARMM-GUI](https://www.charmm-gui.org/?doc=input/martini.solution) with the MARTINI22 force field to construct your protein model.
   - Obtain the following files from the output:
     - `step3_charmm2gmx.pdb` (structure file)
     - `PROA_P.itp` (MARTINI force field file)

2. **Prepare Files**:
   - Place `step3_charmm2gmx.pdb` and `PROA_P.itp` in the same directory as `lib.py` and `convert.py`.

3. **Modify `convert.py`**:
   - Open `convert.py` and `step3_charmm2gmx.pdb` in a text editor.
   - Adjust line numbers in `convert.py` to correspond with `step3_charmm2gmx.pdb`.

4. **Run Conversion Script**:
   ```bash
   python convert.py
   ```

5. **Execute LAMMPS Simulations**:
	- Utilize `equilibration.lmps` and `production.lmps` as input files for LAMMPS.
	- Run the equilibration phase first, followed by the production phase.

6. **Advanced Usage**:
	- Refer to `example2` for incorporating additional components, such as heme groups.
	- Ensure LAMMPS includes bug fixes up to July 2022.

## Citation

If you utilize this codebase or database, please cite:

> Wang, Y; Hernandez R.; _ACS Omega_ **9**, 17667 (2024). [https://doi.org/10.1021/acsomega.4c01868]


## Acknowledgment

This work was supported by 
the National Science Foundation
under Grant No. CHE-2001611, the NSF 884 Center for Sustainable Nanotechnology (CSN). The CSN is 885 part of the Centers for Chemical Innovation Program. 

## License

MultiScaleDPD code and databases are distributed under terms of the [MIT License](https://github.com/rxhernandez/MultiScaleDPD/blob/main/LICENSE).
