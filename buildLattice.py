### buildLattice.py
### Morgan Henderson, September 17, 2020
### Produces a lattice configurations

''' This is the part where you set parameters, edit away! '''

# Define a name for the lattice
name = "Ge_2x2x2-CC"

# Define a lattice constant for your basis
lc = 5.658#5.432009155

# Define the bounding box for your basis
a1 = [1,0,0]
a2 = [0,1,0]
a3 = [0,0,1]

# Define atomic positions within your basis
bs = [[0,0,0],\
    [.25,.25,.25],\
    [.5,0,.5],\
    [.75,.25,.75],\
    [.5,.5,0],\
    [.75,.75,.25],\
    [0,.5,.5],\
    [.25,.75,.75]]

masses = [72.63 for i in range(len(bs))]

# Define the type of each atoms in your basis
types = [1 for i in range(len(bs))]

# Define the dimensions of the supercell to be used
dims = [2,2,2]

# Choose which type of file to write (dat or xyz)
fType = "dat"

''' This is the part that runs things, don't edit probably. '''

# Import useful modules
from numpy import zeros, array, unique

# Define arrays for parameters
base = array([a1,a2,a3])*lc
cell = array(bs)*lc
ls = array(dims)

# Determine system information from inputs
nCell = array(cell).shape[0]
nAtoms = int(nCell*ls.prod())
lims = dims[0]*base[0]+dims[1]*base[1]+dims[2]*base[2]
mset = unique(array(masses))

# Determine atomic positions
pos = zeros((nAtoms,3))
anum = 0
for i in range(dims[0]):
    for j in range(dims[1]):
        for k in range(dims[2]):
            for a in range(nCell):
                cellPos = i*base[0]+j*base[1]+k*base[2]
                pos[anum,:] = cell[a]+cellPos
                anum += 1

# Save configuraiton in requested format
if fType == "dat":

    # Write a .dat configuration file
    preamble = "Start File for LAMMPS\n\n"\
        "%d atoms\n\n"%anum+\
        "%d atom types\n\n"%mset.size+\
        "%5f %5f xlo xhi\n"%(0,lims[0])+\
        "%5f %5f ylo yhi\n"%(0,lims[1])+\
        "%5f %5f zlo zhi\n\n"%(0,lims[2])+\
        "Masses\n\n"
    for m in range(len(mset)):
        preamble += "%d %s\n"%(m+1,mset[m])
    preamble += "\nAtoms\n\n"
    with open(name+".dat","w") as datFile:
        datFile.write(preamble)
        for a in range(nAtoms):
            datFile.write("%d\t%d\t%5f\t%5f\t%5f\n"\
                %(a+1,types[a%nCell],pos[a,0],pos[a,1],pos[a,2]))

elif fType == "xyz":

    # Write a .dat configuration file
    with open(name+".xyz","w") as xyzFile:
        xyzFile.write("id mass x y z")
        for a in range(nAtoms):
            xyzFile.write("%d\t%.3f\t%.4f\t%.4f\t%.4f\n"\
                %(a+1,mset[types[a%nCell]],pos[a,0],pos[a,1],pos[a,2]))

else:
    print("ERROR: fType not understood.")
    quit()
