### gridBuild.py
### Morgan Henderson, September 23, 2020
### Builds a lattice configuration using a grid

''' This is the part where you set parameters, edit away! '''

# Define a name for the lattice
name = "SiFin_M-12x8x6_F-4x8x2"

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

# Define grid dimensions
dims = [12,8,12]

# Assign filled elements of grid (change first line)
fill = [[[(z>=3)*(z<9) or (x>=4)*(x<8)*(z>=1)*(z<11) \
    for z in range(dims[2])] for y in range(dims[1])] for x in range(dims[0])]

# Choose output file type (dat or xyz)
fType = "dat"

''' This is the part that runs things, don't edit probably. '''

# Import useful modules
from writeLattice import writeDat, writeXYZ
from numpy import zeros, array

# Define arrays for parameters
a = array(lc) if len(lc)>1 else lc
rl = array([a1,a2,a3])*a
rb = array(bs)*a
li = array(dims)
fl = array(fill)

# Determine total number of cells and atoms
nb = int(rb.shape[0])
nl = int(fl.sum())
na = int(nl*nb)

# Determine box 
R = li[0]*rl[0]+li[1]*rl[1]+li[2]*rl[2]
m = nl*masses

# Determine atomic positions
r = zeros((na,3))
ID = 0
for i in range(li[0]):
    for j in range(li[1]):
        for k in range(li[2]):
            if fl[i,j,k]:
                for b in range(nb):
                    r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]
                    ID += 1

# Save configuraiton in requested format
if fType == "dat": writeDat(name,m,r,R)
elif fType == "xyz": writeXYZ(name,m,r)
else:
    print("ERROR: requested file type not understood.")
    quit()