#! /usr/bin/env python

"""
This script is a test for RfGapThreePointTTF gap model. 
This model uses T,T',S,S' transit time factors (TTF) 
calculated for the second order of polynomial defined by three points 
of the electric field on the axis Em,E0,Ep for -dz,0.,+dz positions.
This class does not use the Polynomial class. Instead it uses
just E(z) = E0*(1+a*z+b*z^2).
"""

import sys
import math
import random

from bunch import Bunch

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF, RfGapThreePointTTF 
from orbit.sns_linac import Drift

# The classes for the field on axis of the cavity
from orbit.sns_linac.rf_field_readers import SuperFish_3D_RF_FieldReader, RF_AxisFieldAnalysis

def makePhaseNear(phase, phase0):                                                                                                                                                          
	""" It will add or substruct any amount of 360. from phase to get close to phase0 """                                                                                              
	n = int(phase0/360.)                                                                                                                                                               
	phase = phase%360.                                                                                                                                                                 
	min_x = 1.0e+38                                                                                                                                                                    
	n_min = 0                                                                                                                                                                          
	for i0 in range(5):                                                                                                                                                                
		i = i0 - 3                                                                                                                                                                 
		d = math.fabs(phase + 360.*(i+n) - phase0)                                                                                                                                 
		if(d < min_x):                                                                                                                                                             
			n_min = i                                                                                                                                                          
			min_x = d                                                                                                                                                          
	return (phase + 360.*(n_min+n)) 

rf_gap = RfGapThreePointTTF()

rf_frequency = 400.0*1.0e+6   # in Hz

#read the RF cavity field
fReader = SuperFish_3D_RF_FieldReader()
fReader.readFile("../SNS_Cavities_Fields/data/scl_medium_beta_rf_cav_field.dat")

#This particular cavity is not symmetric
zSimmetric = 0
spline = fReader.getAxisEz(zSimmetric)

rf_analysis = RF_AxisFieldAnalysis(spline)
root_pos_arr = rf_analysis.rootAnalysis()
center_pos_arr = rf_analysis.gapCentersAnalysis()
print "roots=",root_pos_arr
print "centers=",center_pos_arr

spline = rf_analysis.getNormilizedSpline()

n_gap_steps = 10
z_min = spline.x(0)
z_max = spline.x(spline.getSize() - 1)
step_size = (z_max - z_min)/(n_gap_steps - 1)

dz = step_size/2.0
drift = Drift()
drift.setLength(dz)

def trackBunch(b,cavity_amp,phase_dgr):
	phase = phase_dgr*math.pi/180.	
	time_init = b.getSyncParticle().time()
	for i in range(n_gap_steps-1):
		zm = z_min + i*step_size
		z0 = zm + dz
		zp = z0 + dz
		Em = cavity_amp*spline.getY(zm)
		E0 = cavity_amp*spline.getY(z0)
		Ep = cavity_amp*spline.getY(zp)
		drift.setLength(dz)
		drift.trackBunch(b)
		time_gap = b.getSyncParticle().time()
		delta_phase = 2*math.pi*(time_gap - time_init)*rf_frequency
		rf_gap.trackBunch(b,dz,Em,E0,Ep,rf_frequency,phase+delta_phase)
		#print "====debug Python level========="
		#b.dumpBunch()
		drift.setLength(dz)
		drift.trackBunch(b)
		#b.dumpBunch()
	return b

#---------------------------------------
#---- let's make bunch ---------
#---------------------------------------
b = Bunch()
print "Part. m=",b.mass()
print "Part. q=",b.charge()
TK = 0.1856                    # in [GeV]
syncPart = b.getSyncParticle()
syncPart.kinEnergy(TK)
beta = syncPart.beta()
gamma = syncPart.gamma() 
mass = syncPart.mass()

c_light = 2.99792458e+8
lmbd = c_light/rf_frequency

#---cavity field
E0 =  20.0e+6       # average field in V/m

b.addParticle(0.001,0.0,0.000,0.,0.,0.)

eKin_min = 1.0e+40
eKin_max =-1.0e+40

k_r_min = 1.0e+40
k_r_max =-1.0e+40

phase_start = 0.
phase_step = 1.0
n_steps = int(360./phase_step) + 1
phase_maxE = 0.
print "#  phase[deg]        Ek[MeV]   k_x " 
for i in range(n_steps):
	phase = phase_start + i*phase_step
	b1 = Bunch()
	b.copyBunchTo(b1)
	b1 = trackBunch(b1,E0,phase)
	k_r = b1.xp(0)/b1.x(0)
	eKin_out = b1.getSyncParticle().kinEnergy()
	if(eKin_min > eKin_out):
		eKin_min = eKin_out
	if(eKin_max < eKin_out):
		eKin_max = eKin_out
		phase_maxE = phase
	if(k_r_min > k_r):
		k_r_min = k_r
	if(k_r_max < k_r):
		k_r_max = k_r
	print " %3d "%i," %5.1f "%phase," %12.6f "%(eKin_out*1000)," %12.6f "%k_r

E0TL = 	(eKin_max - eKin_min)/2.0
beta = syncPart.beta()
gamma = syncPart.gamma() 
mass = syncPart.mass()
print "maximal Ekin at RF phase=",makePhaseNear(phase_maxE,0.)
print " E0TL [MeV] = %12.6f "%(E0TL*1000)
k_r_theory = math.pi*E0TL/(mass*gamma**3*beta**3*lmbd)
print "focusing coef. theory     xp/r = %12.5g "%k_r_theory 
print "focusing coef. RF TTF Gap xp/r = %12.5g "%((k_r_max - k_r_min)/2.0)

print "=========================================="
print "Stop."

