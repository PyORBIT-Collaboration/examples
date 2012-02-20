##############################################################
# This script sets up a collimator and collimates a bunch.
##############################################################

import math
import sys
from bunch import Bunch
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject
from collimator import Collimator
print "Start."

length = 0.5
ma = 3
density_fac = 1.0
shape = 1
a = 0.01
b = 0
c = 0
d = 0
angle = 0

collimator =Collimator(length, ma, density_fac, shape, a, b, c, d, angle)

#------------------------------
#Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "Benchmark_Collimator"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)

#=====track bunch through Collimator Node============

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes") 

collimator.collimateBunch(b, lostbunch)
b.dumpBunch("collimatedbunch.dat")
lostbunch.dumpBunch("lostbunch.dat")
print "Stop."


