### surfaceSi.py
### Morgan Henderson, September 23, 2020
### Builds a surface reconstructed Si lattice configuration using a grid

''' This is the part where you set parameters, edit away! '''

# Define a name for the lattice
name = "SiFin-SR_M-12x8x6_F-4x8x2"

# Set Si lattice constant (0 for default)
aSi = 0

# Set Si mass (0 for default)
mSi = 0

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

# Set Si lattice constant
lc = aSi if aSi!=0 else 5.431

# Determine shifting factor
dNN = lc/2**.5
dBond = lc*3**.5/4
shift = (dNN-dBond)/2**1.5

# Define unit 
rl = array([[1,0,0],[0,1,0],[0,0,1]])*lc
rb = array([[0.0,0.0,0.0],[.25,.25,.25],[.5,0.0,.5],\
    [.75,.25,.75],[.5,.5,0.0],[.75,.75,.25],\
    [0.0,.5,.5],[.25,.75,.75]])*lc
li = array(dims)
fl = array(fill)

# Determine total number of cells and atoms
nb = int(rb.shape[0])
nl = int(fl.sum())
na = int(nl*nb)
m = array([mSi if mSi!=0 else 28.085 for i in range(8)]*nl)

# Determine box boundaries
R = li[0]*rl[0]+li[1]*rl[1]+li[2]*rl[2]

# Calculate atomic positions, looping over grid
ID = 0
r = zeros((na,3))
for i in range(li[0]):
    for j in range(li[1]):
        for k in range(li[2]):
            if fl[i,j,k]:

                # Determine whether this is a surface cell
                top = not fl[i,j,k+1]
                bot = not fl[i,j,k-1]

                # Arrange non-surface cell atoms normally
                if not (top or bot):
                    for b in range(nb):
                        r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]
                        ID += 1

                # Arrange surface cell atoms appropriately
                else:
                    for b in range(nb):
                        rc= (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]

                        # Shift positions for surface atoms
                        if bot and rb[b,2]==rb[:,2].min():
                            right = (round((rc[1]-rc[0])/lc))%2
                            if right: rc += array([-shift,shift,0])
                            else: rc += array([shift,-shift,0])
                            
                        if top and rb[b,2]==rb[:,2].max():
                            right = round(abs(rc[1]-rc[0])/lc-.5)%2
                            if (rc[1]>rc[0]):
                                if right: rc += array([shift,-shift,0])
                                else: rc += array([-shift,shift,0])
                            else:
                                if right: rc += array([-shift,shift,0])
                                else: rc += array([shift,-shift,0])

                        r[ID,:] = rc
                        ID += 1

# Save configuraiton in requested format
if fType == "dat": writeDat(name,m,r,R)
elif fType == "xyz": writeXYZ(name,m,r)
else:
    print("ERROR: requested file type not understood.")
    quit()