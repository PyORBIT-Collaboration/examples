#------------------------------------------------
#ORBIT_MPI and pyORBIT benchmark without spacecharge
#------------------------------------------------

import sys
import math

import orbit_mpi

from bunch import Bunch

from orbit.teapot import TEAPOT_Lattice
from spacecharge import SpaceChargeCalc2p5D
from spacecharge import Boundary2D

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications

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

#------------------------------
#Make a Teapot Lattice
#------------------------------

print "Generate Lattice."
lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD("./LATTICES/TEST_SC_LATTICE.LAT","FODO")
#lattice.readMAD("./LATTICES/Q_0p125.LAT","FODO")

print "lattice length=",lattice.getLength()

# dump initial bunch for ORBIT_MPI input
bunch_pyorbit_to_orbit(lattice.getLength(), b, "orbit_mpi_bunch_input.dat")

#make 2.5D space charge calculator
sizeX = 32
sizeY = 32
sizeZ = 5
xy_ratio = 10.0
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ,xy_ratio)

# boundary 
nBoundaryPoints = 128
N_FreeSpaceModes = 32
boundary_radius = 0.11
boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2*boundary_radius)

#=====track bunch through SC Node============
sc_path_length_min = 0.05
scNode_arr = scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d,boundary)

ACC_TURNS = 1
print("Tracking.")
for i in range(ACC_TURNS):
	lattice.trackBunch(b)
	print "Turn ",i	
	
# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(lattice.getLength(), b, "orbit_mpi_bunch_from_pyorbit_output.dat")

print("STOP.")
