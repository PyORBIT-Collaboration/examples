#-------------------------------------------------------------------------
# This script reads the SuperFish file and creates the SuperFishFieldSource
# and the Runge-Kutta 3D tracker. The tracker tracks the bunch of particles
# through the RF cavity. This example is for 1.5 and 1.8 cm MEBT bunchers.
# The focusing effect is calculated from the acceleration and compared with
# the Runge-Kutta bunch tracking.
# The maximal acceleration coincides with the zero focusing and separates 
# transition from de-focusing (xp/x > 0) to focusing(xp/x < 0).
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
#fReader.readFile("data/mebt_1.8cm_field.dat")
fReader.readFile("data/mebt_1.5cm_field.dat")
(grid2D_Ez,grid2D_Er,grid2D_H) = fReader.makeGrid2DFileds_EzErH()

fieldSource = SuperFishFieldSource()
fieldSource.setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H)

#----------------------------------------------
# RF field parameters 
#----------------------------------------------
rf_freq = 402.50e+6        # in Hz
zSimmetric = 1             # it is symmetric
zOrientation = 1           # the cavity is oriented as in the input file
amplitude = 2.0e+6         # the initial amplitude. It is just a number.
phase = 0.*math.pi/180.    # the initial phase
time_init = 0.             # initial time

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
print "length of Grid2D =",(grid2D_Ez.getMaxX()-grid2D_Ez.getMinX())
print "initial time [sec] =",fieldSource.getTimeInit()

#-------Bunch definition ------------------
b = Bunch()
print "Part. m=",b.mass()
print "Part. q=",b.charge()

TK = 0.0025           # in [GeV]
E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())
c_light = 2.99792458e+8
lmbd = c_light/fieldSource.getFrequency()

syncPart = b.getSyncParticle()
syncPart.kinEnergy(TK)

print "TK[GeV] = ",TK
print "P[GeV/c] = ",P

b.addParticle(0.001,0.0,0.000,0.,0.,0.)

b.compress()

print "initial syncPart (px,py,pz) =",(syncPart.px(),syncPart.py(),syncPart.pz())

length = grid2D_Ez.getMaxX()-grid2D_Ez.getMinX()

tracker = RungeKuttaTracker(length)
#--------------------------------------------------------------------------------
# for the symmetric fields (if zSimmetric == +1) the grid2D has only z = 0,z_max
#--------------------------------------------------------------------------------
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

eKin_min = 1.0e+40
eKin_max =-1.0e+40

k_r_min = 1.0e+40
k_r_max =-1.0e+40

phase_start = 0.
phase_step = 1.0
n_steps = int(360./phase_step) + 1
print "#  phase[deg]        Ek[MeV]   k_x    n_steps" 
for i in range(n_steps):
	phase_dgr = phase_start + i*phase_step
	phase =phase_dgr*math.pi/180.
	fieldSource.setPhase(phase)

	b1 = Bunch()
	b.copyBunchTo(b1)
	tracker.trackBunch(b1,fieldSource)
	n_rk4_steps = tracker.stepsNumber()
	k_r = b1.xp(0)/b1.x(0)
	eKin_out = b1.getSyncParticle().kinEnergy()
	if(eKin_min > eKin_out):
		eKin_min = eKin_out
	if(eKin_max < eKin_out):
		eKin_max = eKin_out
	if(k_r_min > k_r):
		k_r_min = k_r
	if(k_r_max < k_r):
		k_r_max = k_r
	print " %3d "%i," %5.1f "%phase_dgr," %12.6f "%(eKin_out*1000)," %12.6f "%k_r," %3d "%n_rk4_steps

	
E0TL = 	(eKin_max - eKin_min)/2.0
beta = syncPart.beta()
gamma = syncPart.gamma() 
mass = syncPart.mass()
print " E0TL [keV] = %12.6f "%(E0TL*1000*1000)
k_r_theory = math.pi*E0TL/(mass*gamma**3*beta**3*lmbd)
print "focusing coef. theory xp/r    = %12.5g "%k_r_theory 
print "focusing coef. SuperFish xp/r = %12.5g "%((k_r_max - k_r_min)/2.0)

	
print "=========================================="
print "Done."
