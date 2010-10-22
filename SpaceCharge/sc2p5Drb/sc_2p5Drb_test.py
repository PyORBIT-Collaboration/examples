#-----------------------------------------------------
#Creates Space Charge Calculator based on Rick Baartman 
#suggestion and tracks the test bunch through the calculator. 
#-----------------------------------------------------


import sys
import math
import random

import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5Drb

print "Start."

sizeX = 128
sizeY = 128
sizeZ = 10
calc2p5d = SpaceChargeCalc2p5Drb(sizeX,sizeY,sizeZ)


macroSize = 1.0e+13
energy = 1.0
b = Bunch()
b.macroSize(macroSize)
b.getSyncParticle().kinEnergy(energy)
#b.addPartAttr("macrosize")

bunch_radius = 0.005
bunch_length = 200.0

nParts = 100000

for ip in range(nParts):
	r = bunch_radius*math.sqrt(random.random())
	phi = 2*math.pi*random.random()
	x = r*math.sin(phi)
	y = r*math.cos(phi)
	z = bunch_length*(1.0 - 2*random.random())
	"""
	z = 0.5*bunch_length*math.sqrt(random.random())
	if(random.random() > 0.5):
		z = z - bunch_length/2
	else:
		z = bunch_length/2 - z
	"""
	b.addParticle(x,0.,y,0.,z,0.)

print "bunchSize = ",b.getSize()
print "macroSize=",b.macroSize()
print "mass=",b.mass()

pipe_radius = 0.010
slice_length = 0.1

b.dumpBunch("pyorbit_bunch_test_in.dat")

calc2p5d.trackBunch(b,slice_length,pipe_radius)

b.dumpBunch("pyorbit_bunch_test_out.dat")

print "Stop."

