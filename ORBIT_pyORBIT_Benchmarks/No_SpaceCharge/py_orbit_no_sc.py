#------------------------------------------------
#ORBIT_MPI and pyORBIT benchmark without spacecharge
#------------------------------------------------

import sys
import math

import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch

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

#------------------------------
#Make a Teapot Lattice
#------------------------------

print "Generate Lattice."
lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD("./LATTICES/Q_0p125.LAT","FODO")

# set the number of sections in quads to the same as for ORBIT_MPI
for acc_elem in lattice.getNodes():
	if(acc_elem.getType() == "quad teapot"):
		acc_elem.setnParts(5)

print "lattice length=",lattice.getLength()

# dump initial bunch for ORBIT_MPI input
bunch_pyorbit_to_orbit(lattice.getLength(), b, "orbit_mpi_bunch_input.dat")

#=====track bunch ============
ACC_TURNS = 1
print("Tracking.")
for i in range(ACC_TURNS):
	lattice.trackBunch(b)
	print "Turn ",i

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(lattice.getLength(), b, "orbit_mpi_bunch_from_pyorbut_output.dat")

print("STOP.")
