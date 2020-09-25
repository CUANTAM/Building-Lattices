### writeLattice.py
### Morgan Henderson, September 23, 2020
### Functions for writing several types of lattice configuration files

# Write a LAMMPS (.dat) configuration file
def writeDat(fName,m,r,R):

    # Get information directly from inputs
    na = len(r)
    mset = list(dict.fromkeys(m))
    t = [mset.index(mi)+1 for mi in m]
    
    # Construct preamble string
    preamble = "Start File for LAMMPS\n\n"\
        "%d atoms\n\n"%na+\
        "%d atom types\n\n"%len(mset)+\
        "%5f %5f xlo xhi\n"%(0,R[0])+\
        "%5f %5f ylo yhi\n"%(0,R[1])+\
        "%5f %5f zlo zhi\n\n"%(0,R[2])+\
        "Masses\n\n"
    for m in range(len(mset)):
        preamble += "%d %s\n"%(m+1,mset[m])
    preamble += "\nAtoms\n\n"

    # Write preamble and atomic information to file
    with open(fName+".dat","w") as datFile:
        datFile.write(preamble)
        for a in range(na):
            datFile.write("%d\t%d\t%5f\t%5f\t%5f\n"\
                %(a+1,t[a],r[a][0],r[a][1],r[a][2]))

# Write a Cartesian (.xyz) configuration file with IDs and masses
def writeXYZ(fName,m,r):
    N = len(m)
    with open(fName+".xyz","w") as xyzFile:
        xyzFile.write("id mass x y z")
        for a in range(N):
            xyzFile.write("%d\t%.3f\t%.4f\t%.4f\t%.4f\n"\
                %(a+1,m[a],r[a][0],r[a][1],r[a][2]))