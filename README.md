# Building-Lattices

## Introduction

The function of this repository is to write lattice configuration files for crystalline solids. There are two methods by which this can be done: bulk building, in which material unit cells are arranged into a rectangular supercell; and grid building, in which an empty 3D grid is created and filled with unit cells according to specification (allowing for substructuring within the system). A list with short descriptions of each program is given below, followed by sections with more detailed information and instructions.

### Programs
- bulkBuild.py
-- Program used to carry out bulk building
- gridBuild.py
-- Program used to carry out grid building
- surfaceSi.py
-- A grid building program for surface reconstructed Silicon
- writeLattice.py
-- Functions for writing various types of configuration files

In each of the programs listed above, excluding writeLattice.py, there is a section for choosing parameters and a section for carrying out the protocol (clearly marked). Unless you are attempting to extend the function of these programs for your own use, it is not recommended that you make any changes to the latter section.

All programs in this repository require numpy.

## bulkBuild.py

To perform build building, open this file and change the following variables in the parameters section:
- **name:** A name for your system (string)
- **lc:** A lattice constant (float)
- **a1, a2, a3:** Unit cell generating vectors (lists)
- **bs:** Position vectors for each atom in your basis (nested list)
- **ms:** Masses for each atom in your basis (list)
- **dims:** Supercell dimensions as multiples of generating vectors (list)
- **fType:** The file type you'd like to save (string, "xyz" or "dat")

Once these parameters are set, save and run the program:

$: python3 buildLattices.py

## gridBuild.py

## surfaceSi.py

## writeLattice.py
