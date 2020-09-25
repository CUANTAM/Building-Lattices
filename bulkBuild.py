### bulkBuild.py
### Morgan Henderson, September 23, 2020
### Builds a bulk lattice configuration

''' This is the part where you set parameters, edit away! '''

# Define a name for the lattice
name = "Si_CC"

# Define lattice constant(s) (scalar or 3 element list)
lc = 5.432009155

# Define unit cell unit vectors
a1 = [1,0,0]
a2 = [0,1,0]
a3 = [0,0,1]

# Define atomic positions within basis (multiples of lc)
bs = [[0.0,0.0,0.0],\
    [.25,.25,.25],\
    [.5,0.0,.5],\
    [.75,.25,.75],\
    [.5,.5,0.0],\
    [.75,.75,.25],\
    [0.0,.5,.5],\
    [.25,.75,.75]]

# Define mass of each atom in basis
masses = [28.085 for i in range(len(bs))]

# Define supercell dimensions
dims = [1,1,1]

# Choose output file type (dat or xyz)
fType = "dat"

''' This is the part that runs things, don't edit probably. '''

# Import useful modules
from writeLattice import writeDat, writeXYZ
from numpy import zeros, array

# Define arrays for parameters
a = lc if type(lc)!=list else array(lc)
rl = array([a1,a2,a3])*a
rb = array(bs)*a
li = array(dims)

# Determine total number of cells and atoms
nb = int(rb.shape[0])
nl = int(li.prod())
na = int(nl*nb)

# Determine box 
R = li[0]*rl[0]+li[1]*rl[1]+li[2]*rl[2]
m = nl*masses

# Determine atomic positions
r = zeros((na,3))
for i in range(li[0]):
    for j in range(li[1]):
        for k in range(li[2]):
            for b in range(nb):
                ID = (i*li[1]*li[2]+j*li[2]+k)*nb+b
                r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]

# Save configuraiton in requested format
if fType == "dat": writeDat(name,m,r,R)
elif fType == "xyz": writeXYZ(name,m,r)
else:
    print("ERROR: requested file type not understood.")
    quit()