##############################################################
# This script will test the functionality of the EnergyAperture
# class of the "aperture" module
##############################################################

import math
import sys

from bunch import Bunch
from aperture import EnergyAperture

print "Start."

#------------------------------
#Main Bunch init
#------------------------------

b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)


dE_spread = 0.010 

nParts = 10
for ind in range(nParts):
	x = 0.
	xp = 0.
	y = 0.
	yp = 0.
	z = 0.
	dE = (2*dE_spread/nParts)*(ind - nParts/2)
	b.addParticle(x,xp,y,yp,z,dE)

lostbunch = Bunch()

#==== make EnergyAperture class instance

energyAperture = EnergyAperture()

#==== check set get parameters methods
energyAperture.setPosition(111.0)
print "energyAperture            pos =",energyAperture.getPosition()
energyAperture.setMinMaxEnergy(-0.005,+0.005)
print "energyAperture min max energy =",energyAperture.getMinMaxEnergy()

#=====track bunch through the EnergyAperture ============
print "Tracking..."

#---- this  will collect lost particles in the lost bunch
energyAperture.checkBunch(b,lostbunch)

#---- if you do not care about the lost particles you can do this:
#energyAperture.checkBunch(b)

print "==============init bunch          =========="
b.dumpBunch()
print "==============lost particles bunch=========="
lostbunch.dumpBunch()

print "Stop."
