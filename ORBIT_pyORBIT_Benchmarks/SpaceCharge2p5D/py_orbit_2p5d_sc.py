#------------------------------------------------
#ORBIT_MPI and pyORBIT benchmark without spacecharge
#------------------------------------------------

import sys
import math
import random

import orbit_mpi

from bunch import Bunch

from orbit.teapot import TEAPOT_Lattice
from spacecharge import SpaceChargeCalc2p5D
from spacecharge import Boundary2D

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit, bunch_orbit_to_pyorbit
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications

print "Start."
#------------------------------
#Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
b.mass(0.93827231)

INTENSITY = 1.0e+14
nMaxMacroParticles = 100000
macro_size = INTENSITY / nMaxMacroParticles

b.macroSize(macro_size)

energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

#------------------------------
#Make a Teapot Lattice
#------------------------------

print "Generate Lattice."
lattice = TEAPOT_Lattice("sc_lattice")
lattice.readMAD("./LATTICES/Q_0p125.LAT","FODO")
#lattice.readMAD("./LATTICES/TEST_SC_LATTICE.LAT","FODO")

print "lattice length=",lattice.getLength()


"""
#if user want to generagte the bunch by himself
r_bunch = 0.05
l_bunch = lattice.getLength()
x_offset = 0.04
y_offset = 0.00
count = 0
while(count < nMaxMacroParticles):
	x = r_bunch*2.0*(0.5-random.random())
	y = r_bunch*2.0*(0.5-random.random())
	z = l_bunch*2.0*(0.5-random.random())
	if(x**2+y**2 < r_bunch**2):
		b.addParticle(x+x_offset,0.,y+y_offset,0.,z,0.)
		count += 1
"""

# read ORBIT bunch
#bunch_orbit_to_pyorbit(lattice.getLength(),energy, "./DISTRIBUTIONS/Bm_KV_offset_40mm_Uniform",b)
bunch_orbit_to_pyorbit(lattice.getLength(),energy, "./DISTRIBUTIONS/Bm_KV_Uniform",b)
#bunch_orbit_to_pyorbit(lattice.getLength(),energy, "./DISTRIBUTIONS/Bm_50mm_Uniform",b)
#bunch_orbit_to_pyorbit(lattice.getLength(),energy, "./DISTRIBUTIONS/Bm_50mm_offset_40mm_Uniform",b)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "pyorbit_initial.dat")
#make 2.5D space charge calculator
sizeX = 256
sizeY = sizeX
sizeZ = 1
xy_ratio = 1.0
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ,xy_ratio)

# boundary 
nBoundaryPoints = 32
N_FreeSpaceModes = 10
boundary_radius = 0.11
boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2.0*boundary_radius)

#=====track bunch through SC Node============
sc_path_length_min = 0.05
scNode_arr = []
scNode_arr = scLatticeModifications.setSC2p5DAccNodes(lattice,sc_path_length_min,calc2p5d,boundary)
#scNode_arr = scLatticeModifications.setSC2p5DAccNodes(lattice,sc_path_length_min,calc2p5d)
print "number of SC nodes =",len(scNode_arr)

TURNS = 100
print("Tracking.")
for i in range(TURNS):
	lattice.trackBunch(b)
	print "Turn ",i	
	
# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(lattice.getLength(), b, "pyorbit_final.dat")

print("STOP.")
