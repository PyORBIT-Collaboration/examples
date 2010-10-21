#------------------------------------------------
#ORBIT_MPI and pyORBIT benchmark without spacecharge
#------------------------------------------------

import sys
import math

import orbit_mpi

from bunch import Bunch

from spacecharge import SpaceChargeCalc2p5D
from spacecharge import Boundary2D

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

print "Start."
#------------------------------
#Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
b.readBunch("pyorbit_bunch_input.dat")
b.mass(0.93827231)
b.macroSize(1.0)

energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

lattice_length = 2.0

# dump initial bunch for ORBIT_MPI input
bunch_pyorbit_to_orbit(lattice_length, b, "orbit_mpi_bunch_input.dat")

#make 2.5D space charge calculator
sizeX = 32
sizeY = 32
sizeZ = 5
xy_ratio = 1.0
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ,xy_ratio)

# boundary 
nBoundaryPoints = 128
N_FreeSpaceModes = 32
boundary_radius = 0.011
boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2*boundary_radius)

#=====track bunch through SC calculator ============
sc_length = 0.5
calc2p5d.trackBunch(b,sc_length,boundary)

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(lattice_length, b, "orbit_mpi_bunch_from_pyorbut_output.dat")

print("STOP.")
