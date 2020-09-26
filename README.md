# Building-Lattices

## Introduction

The function of this repository is to write lattice configuration files for crystalline solids. There are two methods by which this can be done: bulk building, in which material unit cells are arranged into a rectangular supercell; and grid building, in which an empty 3D grid is created and filled with unit cells according to specification (allowing for substructuring within the system). A list with short descriptions of each program is given below, followed by sections with more detailed information and instructions.

### Programs

- bulkBuild.py
  - Program used to carry out bulk building
- gridBuild.py
  - Program used to carry out grid building
- surfaceSi.py
  - A grid building program for surface reconstructed Silicon
- writeLattice.py
  - Functions for writing various types of configuration files

In each of the programs listed above, excluding writeLattice.py, there is a section for choosing parameters and a section for carrying out the protocol (clearly marked). Unless you are attempting to extend the function of these programs for your own use, it is not recommended that you make any changes to the latter section.

All programs in this repository require numpy.

## bulkBuild.py

During bulk building, a perfect crystal unit cell is repeated a specified number of times in three dimensions in order to produce a rectangular crystal supercell. To perform bulk building, open buildBulk.py and change the following variables in the parameters section:
- **name:** A name for your system (string)
- **lc:** A lattice constant (float)
- **a1, a2, a3:** Unit cell generating vectors (lists)
- **bs:** Position vectors for each atom in your basis (nested list)
- **ms:** Masses for each atom in your basis (list)
- **dims:** Supercell dimensions as multiples of generating vectors (list)
- **fType:** The file type you'd like to save (string, "xyz" or "dat")

Once all parameters have been set, save and run the program:

$: python3 bulkBuild.py

## gridBuild.py

During grid building, an empty unit cell grid of specified dimensions is filled with perfect crystal unit cells according to specified conditions. Regions of the grid may be left empty, allowing for the creation of substructured systems (e.g. membranes, wires, and substrates with surface features). To perform grid building, open gridBuild.py and change the following variables in the parameters section:
- **name:** A name for your system (string)
- **lc:** A lattice constant (float)
- **a1, a2, a3:** Unit cell generating vectors (lists)
- **bs:** Position vectors for each atom in your basis (nested list)
- **ms:** Masses for each atom in your basis (list)
- **dims:** Grid dimensions as multiples of generating vectors (list)
- **fill:** List of grid elements that should be filled (list)
- **fType:** The file type you'd like to save (string, "xyz" or "dat")

The variable "fill" can become very cumbersome for large systems when defined element by element. In most cases, it is much quicker and easier to define "fill" using logical conditions on the x, y, and z coordinates of the grid elements. This can be done as follows (starting from grid definition):

```
# Define grid dimensions
dims = [8,8,8]

# Assign filled elements of grid
fill = [[[(z>=2)*(z<6) \
    for z in range(dims[2])] for y in range(dims[1])] for x in range(dims[0])]
```

In this example, a cubic 8x8x8 unit cell grid has been defined, but only the elements with z>=2 and z<6 are filled, creating an 8x8x4 membrane of unit cells with 2 unit cell vacant regions above and below. Using this flexible syntax, many different kinds of structures (such as the ones mentioned above can be produced.

Once all parameters have been set, save and run the program:

$: python3 gridBuild.py

## surfaceSi.py

During grid building, an empty unit cell grid of specified dimensions is filled with perfect crystal unit cells according to specified conditions. Regions of the grid may be left empty, allowing for the creation of substructured systems (e.g. membranes, wires, and substrates with surface features). To perform grid building, open gridBuild.py and change the following variables in the parameters section:
- **name:** A name for your system (string)
- **lc:** A lattice constant (float)
- **a1, a2, a3:** Unit cell generating vectors (lists)
- **bs:** Position vectors for each atom in your basis (nested list)
- **ms:** Masses for each atom in your basis (list)
- **dims:** Grid dimensions as multiples of generating vectors (list)
- **fill:** List of grid elements that should be filled (list)
- **fType:** The file type you'd like to save (string, "xyz" or "dat")

The variable "fill" can become very cumbersome for large systems when defined element by element. In most cases, it is much quicker and easier to define "fill" using logical conditions on the x, y, and z coordinates of the grid elements. This can be done as follows (starting from grid definition):

```
# Define grid dimensions
dims = [8,8,8]

# Assign filled elements of grid
fill = [[[(z>=2)*(z<6) \
    for z in range(dims[2])] for y in range(dims[1])] for x in range(dims[0])]
```

In this example, a cubic 8x8x8 unit cell grid has been defined, but only the elements with z>=2 and z<6 are filled, creating an 8x8x4 membrane of unit cells with 2 unit cell vacant regions above and below. Using this flexible syntax, many different kinds of structures (such as the ones mentioned above can be produced.

Once all parameters have been set, save and run the program:

$: python3 gridBuild.py

## writeLattice.py

The writeLattice.py program contains several functions for writing configuration data calculated by other programs in this repository to files that are readable by various simulation software libraries (e.g. LAMMPS .dat files). These functions are called from within the active part of the other scripts, and in most cases do not require the user's attention. If a configuration is needed for a specific simulation software that will not take .xyz or .dat formats, a tailored function will need to be added to this program, and the active part of the building program will need to be edited to call the new function.
