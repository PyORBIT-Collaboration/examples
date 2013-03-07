#-------------------------------------------------------------------------
# This script reads the SuperFish file and creates the SuperFishFieldSource
#  and the Runge-Kutta 3D tracker. The tracker tracks the bunch of particles
#  through the RF cavity.
#--------------------------------------------------------------------------

import sys
import math

from bunch import Bunch
from spacecharge import Grid2D
from orbit.sns_linac.rf_field_readers import SuperFish_3D_RF_FieldReader, RF_AxisFieldAnalysis
from trackerrk4 import RungeKuttaTracker

# from linac import the 
from linac import SuperFishFieldSource

fReader = SuperFish_3D_RF_FieldReader()
fReader.readFile("data/scl_medium_beta_rf_cav_field.dat")
(grid2D_Ez,grid2D_Er,grid2D_H) = fReader.makeGrid2DFileds_EzErH()

fieldSource = SuperFishFieldSource()
fieldSource.setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H)

#----------------------------------------------
# RF field parameters 
#----------------------------------------------
rf_freq = 805.0e+6        # in Hz
zSimmetric = 0            # it is not symmetric
zOrientation = +1         # the cavity is oriented as in the input file
amplitude = 1000000.           # the initial amplitude. It means nothing.
phase = 30.*math.pi/180.  # the initial phase
fieldCenterPos = (grid2D_Ez.getMinX()+grid2D_Ez.getMaxX())/2.0
directionZ = +1
symmetry = 0

fieldSource = SuperFishFieldSource()
fieldSource.setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H)
fieldSource.setFrequency(rf_freq)
fieldSource.setAmplitude(amplitude)
fieldSource.setPhase(phase)
fieldSource.setFieldCenterPos(fieldCenterPos)
fieldSource.setDirectionZ(directionZ)
fieldSource.setSymmetry(symmetry)

print "frequnecy=",fieldSource.getFrequency()
print "amplitude=",fieldSource.getAmplitude()
print "phase=",fieldSource.getPhase()*180./math.pi
print "fieldCenterPos = ",fieldSource.getFieldCenterPos()
print "directionZ=",fieldSource.getDirectionZ()
print "symmetry = ",fieldSource.getSymmetry()
print "min max Z =",(grid2D_Ez.getMinX(),grid2D_Ez.getMaxX())
print "length =",(grid2D_Ez.getMaxX()-grid2D_Ez.getMinX())
#-------Bunch definition ------------------
b = Bunch()
print "Part. m=",b.mass()
print "Part. q=",b.charge()

TK = 0.1856           # in [GeV]
E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())
c_light = 2.99792458e+8


print "TK[GeV] = ",TK
print "P[GeV/c] = ",P

b.addParticle(0.,0.,0.,0.,0.,0.)
b.compress()
syncPart = b.getSyncParticle()
syncPart.kinEnergy(TK)
print "initial syncPart (px,py,pz) =",(syncPart.px(),syncPart.py(),syncPart.pz())

length = grid2D_Ez.getMaxX()-grid2D_Ez.getMinX()

tracker = RungeKuttaTracker(length)

tracker.entrancePlane(0,0,-1.,grid2D_Ez.getMinX())
tracker.exitPlane(0,0,1.,-grid2D_Ez.getMaxX())

print "Entrance plane (a,b,c,d)=",tracker.entrancePlane()
print "Exit     plane (a,b,c,d)=",tracker.exitPlane()
print "Length[m]=",tracker.length()

#b.dumpBunch()

tracker.trackBunch(b,fieldSource)

b.dumpBunch()



print "=========================================="
print "Done."

