#-------------------------------------------------------------------------
# This script reads the SuperFish file and creates the SuperFishFieldSource
# and the Runge-Kutta 3D tracker. The tracker tracks the bunch of particles
# through the RF cavity. The tracker will be divided onto the Nparts parts.
# The results should be the same no matter how many parts we use.
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
#fReader.readFile("data/scl_high_beta_rf_cav_field.dat")
(grid2D_Ez,grid2D_Er,grid2D_H) = fReader.makeGrid2DFileds_EzErH()

fieldSource = SuperFishFieldSource()
fieldSource.setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H)

#----------------------------------------------
# RF field parameters 
#----------------------------------------------
rf_freq = 805.0e+6                # in Hz
zSimmetric = 0                    # it is not symmetric
zOrientation = -1                 # the cavity is oriented as in the input file
amplitude = 20.0e+6               # the initial amplitude. It is just a number.
phase = (270.-90.)*math.pi/180.  # the initial phase
time_init = 0.                    # initial time

fieldSource = SuperFishFieldSource()
fieldSource.setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H)
fieldSource.setFrequency(rf_freq)
fieldSource.setAmplitude(amplitude)
fieldSource.setPhase(phase)
fieldSource.setDirectionZ(zOrientation)
fieldSource.setSymmetry(zSimmetric)
fieldSource.setTimeInit(time_init)

print "frequnecy=",fieldSource.getFrequency()
print "amplitude=",fieldSource.getAmplitude()
print "phase=",fieldSource.getPhase()*180./math.pi
print "fieldCenterPos = ",fieldSource.getFieldCenterPos()
print "directionZ=",fieldSource.getDirectionZ()
print "symmetry = ",fieldSource.getSymmetry()
print "min max Z =",(grid2D_Ez.getMinX(),grid2D_Ez.getMaxX())
print "min max R =",(grid2D_Ez.getMinY(),grid2D_Ez.getMaxY())
print "length =",(grid2D_Ez.getMaxX()-grid2D_Ez.getMinX())
print "length of the filed =",fieldSource.getLength()
print "average filed [MV/m] =",fieldSource.getAvgField()/1.0e+6
print "AvgField*Length [kV] =",fieldSource.getAvgField()*fieldSource.getLength()/1000.
print "initial time [sec] =",fieldSource.getTimeInit()

#-------Bunch definition ------------------
b = Bunch()
print "Part. m=",b.mass()
print "Part. q=",b.charge()
syncPart = b.getSyncParticle()

TK = 0.400           # in [GeV]
E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())
c_light = 2.99792458e+8
lmbd = c_light/fieldSource.getFrequency()

syncPart.kinEnergy(TK)

print "TK[GeV] = ",TK
print "P[GeV/c] = ",P

print "lambda [mm] =",lmbd*1000.
print "beta*lambda/360 deg [mm/deg] =",lmbd*syncPart.beta()*1000./360.

b.addParticle(0.002,0.001,0.002,0.001,0.005,0.)
b.addParticle(0.0,0.0,0.0,0.,-0.005,0.)

b.compress()

print "initial syncPart (px,py,pz) =",(syncPart.px(),syncPart.py(),syncPart.pz())

length = grid2D_Ez.getMaxX()-grid2D_Ez.getMinX()

tracker = RungeKuttaTracker(length)
#-------------------------------------------------------------------------------
# for the symmetric fields (if zSimmetric == +1) the grid2D has only z = 0,z_max
#-------------------------------------------------------------------------------
if(fieldSource.getSymmetry() == 1):
	tracker.entrancePlane(0,0,-1.,-grid2D_Ez.getMaxX())
else:
	tracker.entrancePlane(0,0,-1.,grid2D_Ez.getMinX())
tracker.exitPlane(0,0,1.,-grid2D_Ez.getMaxX())
tracker.spatialEps(0.0000001)
tracker.stepsNumber(60)

print "Entrance plane (a,b,c,d)=",tracker.entrancePlane()
print "Exit     plane (a,b,c,d)=",tracker.exitPlane()
print "Length[m]=",tracker.length()
	
b1 = Bunch()
b.copyBunchTo(b1)

nParts = 10
s_start = grid2D_Ez.getMinX()
if(fieldSource.getSymmetry() == 1):
	s_start = - grid2D_Ez.getMaxX()
s_end = grid2D_Ez.getMaxX()
s_step = (s_end - s_start)/nParts

time_ini = 0.
syncPart = b1.getSyncParticle()

for i in range(nParts):
	s0 = s_start + i*s_step
	s1 = s0+s_step
	tracker.entrancePlane(0,0,-1.,s0)
	tracker.exitPlane(0,0,1.,-s1)
	#print "============================================================ i=",i
	#print "Entrance plane (a,b,c,d)=",tracker.entrancePlane()
	#print "Exit     plane (a,b,c,d)=",tracker.exitPlane()	
	fieldSource.setTimeInit(syncPart.time())
	tracker.trackBunch(b1,fieldSource)

b1.dumpBunch()
	
print "=========================================="
print "Done."
