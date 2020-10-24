### surfaceSi.py
### Morgan Henderson
### October 22, 2020
### Builds Si lattice configurations with surface modifcations using a grid

''' This is the part where you set parameters, edit away! '''

# Define a name for the lattice
name = "SiMemP_1x6x4-CC"

# Indicate whether to reconstruct (100) surfaces with 2x1 dimers
reconstruct = False

# Indicate whether to passivate (100) surfaces with H
passivate = True

# Set Si lattice constant (0 for default)
aSi = 0

# Set Si mass (0 for default)
mSi = 0

# Define grid dimensions
dims = [1,6,6]

# Assign filled elements of grid (change first line)
fill = [[[(z>=1)*(z<5)  \
    for z in range(dims[2])] for y in range(dims[1])] for x in range(dims[0])]

# Choose output file type (dat or xyz)
fType = "xyz"

# Choose whether to check surface bond lengths
check = True

''' This is the part that runs things, don't edit probably. '''

# Handle any errors based on provided inputs
if 0 in dims:
    print("\nERROR: System dimensions cannot be 0.\n")
if reconstruct and (dims[0]==1 or dims[2]==1):
    print("\nERROR: Cannot properly reconstruct systems with unitary x- or y-dimension.\n")
    quit()

# Import useful modules
from writeLattice import writeDat, writeXYZ
from numpy import pi, sin, cos, zeros, array, rint

# Set Si lattice constant
aSi = aSi if aSi!=0 else 5.431
mSi = mSi if mSi!=0 else 28.0855

# Define constants for 2x1 reconstruction
if reconstruct:
    dNN = aSi/2**.5 # Si nearest-neighbor distance
    dDimer = aSi*3**.5/4+.05*passivate # Dimer bond length
    sSi = (dNN-dDimer)/2**1.5 # Shifting factor for surface Si atoms

# Define constants for H passivation
if passivate:
    mH = 1.00784 #Hydrogen mass
    theta = (69.8 if reconstruct else 39.0)*(pi/180) # Hydrogen bond angle (from surface)
    dH = 1.54 if reconstruct else 1.50 # Hydrogen bond length, angstroms
    zH = dH*sin(theta) # Hydrogen z-position relative to surface, angstroms
    sH = dH*cos(theta)/2**.5 # Hydrogen xy-position relative to surface, angstroms

# Define unit cell vectors and basis atom positions
rl = array([[1,0,0],[0,1,0],[0,0,1]])*aSi
rb = array([[0.0,0.0,0.0],[.25,.25,.25],[.5,0.0,.5],\
    [.75,.25,.75],[.5,.5,0.0],[.75,.75,.25],\
    [0.0,.5,.5],[.25,.75,.75]])*aSi

# Determine total number of Si cells and atoms
nb = int(rb.shape[0])
fl = array(fill)
nl = int(fl.sum())
nSi = int(nl*nb)

# Detect [0 0 1] surfaces in the system, determine number of surface Si atoms
top = array([[[fl[i,j,k] and not fl[i,j,k+1] \
    for k in range(dims[2])] for j in range(dims[1])] for i in range(dims[0])])
bot = array([[[fl[i,j,k] and not fl[i,j,k-1] \
    for k in range(dims[2])] for j in range(dims[1])] for i in range(dims[0])])
nSurf = 2*(top.sum()+bot.sum())
if passivate: nH = nSurf*(1+(not reconstruct))
nAtoms = nSi if not passivate else nSi+nH

# Define variables for construction loop
ID = 0
r = zeros((nAtoms,3))
m = zeros(nAtoms)
surfaceSi = zeros(nAtoms,dtype=int)

# Determine positions for a reconstructed & passivated system
if (reconstruct and passivate):
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                if fl[i,j,k]:

                    # Arrange non-surface atoms normally
                    if not (top[i,j,k] or bot[i,j,k]):
                        for b in range(nb):
                            r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]
                            m[ID] = mSi
                            ID += 1

                    # Arrange surface atoms according to Northrup (1991)
                    else:
                        for b in range(nb):
                            rSi = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]

                            # Shift bottom surface Si atoms, add H atoms for passivation
                            if bot[i,j,k] and rb[b,2]==rb[:,2].min():
                                surfaceSi[ID] = 1
                                right = (round((rSi[1]-rSi[0])/aSi))%2
                                rSi += array([sSi-2*right*sSi,-sSi+2*right*sSi,0])
                                rH = rSi + array([-sH+2*right*sH,sH-2*right*sH,-zH])
                                r[ID:ID+2] = [rSi,rH]
                                m[ID:ID+2] = [mSi,mH]
                                ID += 2
                                
                            # Shift top surface Si atoms, add H atoms for passivation
                            elif top[i,j,k] and rb[b,2]==rb[:,2].max():
                                surfaceSi[ID] = 1
                                right = round(abs(rSi[1]-rSi[0])/aSi-.5)%2

                                # Determine whether dimer row is odd/even (top only)
                                if (rSi[1]>rSi[0]):
                                    rSi += array([-sSi+2*right*sSi,sSi-2*right*sSi,0])
                                    rH = rSi + array([sH-2*right*sH,-sH+2*right*sH,zH])
                                else:
                                    rSi += array([sSi-2*right*sSi,-sSi+2*right*sSi,0])
                                    rH = rSi + array([-sH+2*right*sH,sH-2*right*sH,zH])
                                r[ID:ID+2] = [rSi,rH]
                                m[ID:ID+2] = [mSi,mH]
                                ID += 2

                            else:
                                r[ID,:] = rSi
                                m[ID] = mSi
                                ID += 1

# Determine positions for a reconstructed system
elif reconstruct and not passivate:
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                if fl[i,j,k]:

                    # Arrange non-surface atoms normally
                    if not (top[i,j,k] or bot[i,j,k]):
                        for b in range(nb):
                            r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]
                            m[ID] = mSi
                            ID += 1

                    # Arrange surface atoms according to Appelbaum (1976)
                    else:
                        for b in range(nb):
                            rSi = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]

                            # Shift bottom surface Si atoms, add H atoms for passivation
                            if bot[i,j,k] and rb[b,2]==rb[:,2].min():
                                right = (round((rSi[1]-rSi[0])/aSi))%2
                                rSi += array([sSi-2*right*sSi,-sSi+2*right*sSi,0])
                                surfaceSi[ID] = 1
                                r[ID,:] = rSi
                                m[ID] = mSi
                                ID += 1
                                
                            # Shift top surface Si atoms, add H atoms for passivation
                            elif top[i,j,k] and rb[b,2]==rb[:,2].max():
                                right = round(abs(rSi[1]-rSi[0])/aSi-.5)%2

                                # Determine whether dimer row is odd/even (top only)
                                if (rSi[1]>rSi[0]):
                                    rSi += array([-sSi+2*right*sSi,sSi-2*right*sSi,0])
                                else:
                                    rSi += array([sSi-2*right*sSi,-sSi+2*right*sSi,0])
                                surfaceSi[ID] = 1
                                r[ID,:] = rSi
                                m[ID] = mSi
                                ID += 1

                            else:
                                r[ID,:] = rSi
                                m[ID] = mSi
                                ID += 1

# Determine positions for a passivated system
elif passivate and not reconstruct:
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                if fl[i,j,k]:

                    # Arrange non-surface atoms normally
                    if not (top[i,j,k] or bot[i,j,k]):
                        for b in range(nb):
                            r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]
                            m[ID] = mSi
                            ID += 1

                    # Arrange surface atoms according to Northrup (1991)
                    else:
                        for b in range(nb):
                            rSi = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]

                            # Shift bottom surface Si atoms, add H atoms for passivation
                            if bot[i,j,k] and rb[b,2]==rb[:,2].min():
                                surfaceSi[ID] = 1
                                rH1 = rSi + array([sH,-sH,-zH])
                                rH2 = rSi + array([-sH,sH,-zH])
                                r[ID:ID+3] = [rSi,rH1,rH2]
                                m[ID:ID+3] = [mSi,mH,mH]
                                ID += 3
                                
                            # Shift top surface Si atoms, add H atoms for passivation
                            elif top[i,j,k] and rb[b,2]==rb[:,2].max():
                                surfaceSi[ID] = 1
                                rH1 = rSi + array([sH,-sH,zH])
                                rH2 = rSi + array([-sH,sH,zH])
                                r[ID:ID+3] = [rSi,rH1,rH2]
                                m[ID:ID+3] = [mSi,mH,mH]
                                ID += 3

                            else:
                                r[ID,:] = rSi
                                m[ID] = mSi
                                ID += 1

# Determine positions for ordinary system
elif not (reconstruct or passivate):
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                if fl[i,j,k]:

                    # Arrange all atoms normally
                    for b in range(nb):
                        r[ID,:] = (i*rl[0]+j*rl[1]+k*rl[2])+rb[b]
                        m[ID] = mSi
                        ID += 1

# Determine box boundaries
R = dims[0]*rl[0]+dims[1]*rl[1]+dims[2]*rl[2]

# If desired, check bond lengths against experimental values
if check and (reconstruct or passivate):

    # Identify surface Si atoms
    rSurf = array([r[i] for i in range(nAtoms) if surfaceSi[i]])

    if reconstruct:

        # Define counting variable for properly formed dimers
        nDimer = 0

    if passivate:

        # Identify H atoms, define counting variable for properly formed hydrogen bonds
        rH = array([r[i] for i in range(nAtoms) if m[i]==mH])
        nHBond = 0

    # Loop over surface Si atoms
    for i in range(nSurf):

        if reconstruct:

            # Loop over pairs of surface Si atoms
            for j in range(i):
                rij = rSurf[i]-rSurf[j]

                # Apply periodic boundary conditions to separation vector
                rPBC = rij/R
                rPBC -= rint(rPBC)
                rij = rPBC*R

                # Determine dimer bond lengths, check against Appelbaum (1976)/Northrup (1991)
                dij = (rij.dot(rij))**.5
                if round(dij,2)==round(dDimer,2):
                    nDimer += 1

        if passivate:

            # Loop over surface Si/H atoms pairs
            for k in range(nH):
                rik = rSurf[i]-rH[k]

                # Apply periodic boundary conditions to separation vector
                rPBC = rik/R
                rPBC -= rint(rPBC)
                rik = rPBC*R

                # Determine hydrogen bond lengths, check against Northrup (1991)
                dik = (rik.dot(rik))**.5
                if round(dik,2)==round(dH,2):
                    nHBond += 1

    # Indicate surface bond formation status and bond lengths
    if reconstruct and passivate: status = "properly" if (nDimer==0.5*nSurf and nHBond==nH) else "improperly"
    elif reconstruct and not passivate: status = "properly" if nDimer==0.5*nSurf else "improperly"
    elif passivate and not reconstruct: status = "properly" if nHBond==nH else "improperly"        
    print("\nSurface bonds formed %s:\n\t"%status+\
            "Number of Surface Si Atoms = %d;"%nSurf)
    if reconstruct:
        print("\tNumber of Surface Dimer Bonds = %d;\n\t"%nDimer+\
            "Surface Dimer Bond Lengths = %.2f Angstroms\n"%round(dDimer,2))
    if passivate:
        print("\tNumber of Surface Hydrogen Bonds = %d\n\t"%nHBond+\
            "Surface Hydrogen Bond lengths = %.2f Angstroms\n"%round(dH,2))

elif check and not (reconstruct or passivate):
    print("WARNING: No surface modifications to check.")

# Save configuraiton in requested format
if fType == "dat": writeDat(name,m,r,R)
elif fType == "xyz": 
    if not passivate:
        t = ["Si" for i in range(nAtoms)]
        writeXYZ(name,r,t)
    else:
        t = ["Si" if m[i]==mSi else "H" for i in range(nAtoms)]
        writeXYZ(name,r,t)

else:
    print("ERROR: requested file type not understood.")
    quit()