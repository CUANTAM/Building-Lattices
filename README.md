# Building-Lattices

The file buildLattice.py is a simple program that writes a
configuration file in either a .dat or .xyz format. To use
the program, edit the section with the header:

''' This is the part where you set parameters, edit away! '''

The parameters that must be provided are:
- A name for your system
- A lattice constant
- Unit cell vectors (a1,a2,a3)
- Basis vectors (bs=b1,...,bn)
- Masses for each basis atom
- Types (1,2,...N) for each basis atom
- Supercell dimensions as a multiple of unit cell vectors
- The file type (xyz or dat) you'd like to save

Once the parameters are set, simply save and run the program:

$: python3 buildLattices.py

and a new file, "(name).(fType)", should appear in your
directory. The parameters currently set in the file are
for a 2x2x2 conventional unit cell germanium lattice.
