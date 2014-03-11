##############################################################
# This script test the particle id numbers.
##############################################################

import math
import sys
from bunch import Bunch
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject
from orbit.bunch_utils import ParticleIdNumber
print "Start."

#------------------------------
#Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "test particle id number"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)

#------------------------------
# Small Test Bunch init
#------------------------------
btest = Bunch()
print "Read Bunch."
runName = "test particle id number"

btest.mass(0.93827231)
btest.macroSize(1.0e+1)
energy = 1.0 #Gev
btest.readBunch("parts.dat",100)
btest.getSyncParticle().kinEnergy(energy)

particleidnumber = ParticleIdNumber()

particleidnumber.addParticleIdNumbers(b, 0)
particleidnumber.addParticleIdNumbers(btest)

#add bunches together
btest.addParticlesTo(b)

#=====track bunch through Collimator Node============

b.dumpBunch("idbunchmain.dat")
btest.dumpBunch("idbunchtest.dat")
print "Stop."


