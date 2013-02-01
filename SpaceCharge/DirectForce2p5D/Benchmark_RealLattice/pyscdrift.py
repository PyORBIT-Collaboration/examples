
import math
import sys

from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.space_charge.directforce2p5d import directforceAccNodes, directforceLatticeModifications
from spacecharge import SpaceChargeForceCalc2p5D

print "Start."

#=====Make a Teapot style lattice======

lattice = teapot.TEAPOT_Ring()
print "Read MAD."
lattice.readMAD("../MAD_Lattice/SNSring_pyOrbitBenchmark.LAT","RING")
print "Lattice=",lattice.getName()," length [m] =",lattice.getLength()," nodes=",len(lattice.getNodes())


#------------------------------
#Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "Benchmark_SpaceCharge"

total_macroSize=1.0e+14
b.mass(0.93827231)
energy = 1.0 #Gev
b.readBunch("KV.dat", 10000)
b.getSyncParticle().kinEnergy(energy)
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize/nParticlesGlobal)

#-----------------------------------
# Add Direct Force Space Charge node
#-----------------------------------

nMacrosMin = 1
sc_path_length_min = 0.00000001
sizeX = 9   #number of grid points in horizontal direction
sizeY = 9   #number of grid points in vertical direction
sizeZ = 1   #number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeForceCalc2p5D(sizeX,sizeY,sizeZ)
#calc2p5d.trackBunch(b, lattice_length)
directforceLatticeModifications.setDirectForce2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)

#
#------------------------------------

paramsDict = {}
paramsDict["bunch"]= b

lattice.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch1.dat")
lattice.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch2.dat")
lattice.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch3.dat")
lattice.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch4.dat")
lattice.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch5.dat")
for i in range(95):
	lattice.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch100.dat")


print "Stop."



