##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit

print "Start."

#=====Main bunch parameters============

intensity = 1.0e14
turns = 1000.0
macrosperturn = 260
macrosize = intensity/turns/macrosperturn

b = Bunch()
b.mass(0.93827231)
b.macroSize(macrosize)
energy = 1.0 #Gev
b.readBunch("kv.dat")
b.getSyncParticle().kinEnergy(energy)

paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b
lostbunch.addPartAttr("LostParticleAttributes") 

#=====Make a Teapot style lattice======

teapot_latt = teapot.TEAPOT_Ring()
print "Read MAD."
teapot_latt.readMAD("MAD_Lattice/SNSring_pyOrbitBenchmark.LAT","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())
lattlength = teapot_latt.getLength()

bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "kv_orbit.dat")

print "Stop."



