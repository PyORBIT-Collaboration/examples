##############################################################
# This script will test the functionality of the PhaseAperture
# class of the "aperture" module
##############################################################

import math
import sys

from bunch import Bunch
from aperture import PhaseAperture

print "Start."

#------------------------------
#Main Bunch init
#------------------------------

rf_frequency = 402.5e+6
c = 2.99792458e+8

b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)
beta = b.getSyncParticle().beta()

lambda_rf = c/rf_frequency
phase_to_z_coeff = lambda_rf*beta/360.

nParts = 10
for ind in range(nParts):
	x = 0.
	xp = 0.
	y = 0.
	yp = 0.
	z = (360./nParts)*(ind - nParts/2)*phase_to_z_coeff
	dE = 0.
	b.addParticle(x,xp,y,yp,z,dE)

lostbunch = Bunch()

#==== make PhaseAperture class instance

phaseAperture = PhaseAperture(2*rf_frequency)

#==== check set get parameters methods
print "phaseAperture  1   frequency =",phaseAperture.getRfFrequency()
phaseAperture.setRfFrequency(rf_frequency)
print "phaseAperture  2   frequency =",phaseAperture.getRfFrequency()
phaseAperture.setPosition(111.0)
print "phaseAperture            pos =",phaseAperture.getPosition()
phaseAperture.setMinMaxPhase(-100.,+100.)
print "phaseAperture min max phases =",phaseAperture.getMinMaxPhase()

#=====track bunch through the PhaseAperture ============
print "Tracking..."

#---- this  will collect lost particles in the lost bunch
phaseAperture.checkBunch(b,lostbunch)

#---- if you do not care about the lost particles you can do this:
#phaseAperture.checkBunch(b)

print "==============init bunch          =========="
b.dumpBunch()
print "==============lost particles bunch=========="
lostbunch.dumpBunch()

print "Stop."
