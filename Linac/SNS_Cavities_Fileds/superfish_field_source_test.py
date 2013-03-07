#-------------------------------------------------------------------------
# This script reads the SuperFish file and creates the SuperFishFieldSource
# which is a wrapper around C++ class. This class is a field source for 
# 3D Runge-Kutta tracker.
#--------------------------------------------------------------------------

import sys
import math

from spacecharge import Grid2D

from orbit.sns_linac.rf_field_readers import SuperFish_3D_RF_FieldReader, RF_AxisFieldAnalysis

# from linac import the 
from linac import SuperFishFieldSource

fReader = SuperFish_3D_RF_FieldReader()
fReader.readFile("data/scl_medium_beta_rf_cav_field.dat")
(grid2D_Ez,grid2D_Er,grid2D_H) = fReader.makeGrid2DFileds_EzErH()


rf_freq = 805.0e+6
zSimmetric = 0
zOrientation = +1
amplitude = 10.
phase = 30.*math.pi/180.
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
print "max min Z =",(grid2D_Ez.getMinX(),grid2D_Ez.getMaxX())

(Ex,Ey,Ez,Bx,By,Bz) = fieldSource.getElectricMagneticField(0.01,0.01,0.27,0.)
print "(Ex,Ey,Ez,Bx,By,Bz)=",(Ex,Ey,Ez,Bx,By,Bz)

"""
#memory leak check
count = 1
while(1 < 2):
	(grid2D_Ez,grid2D_Er,grid2D_H) = fReader.makeGrid2DFileds_EzErH()
	fieldSource.setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H)
	(grid2D_Ez,grid2D_Er,grid2D_H) = fieldSource.getGrid2D_Fields()
	count += 1
	if(count % 10 == 0):
		print "i=",count
"""

print "=========================================="
print "Done."

