#! /usr/bin/env python

"""
This script is a test for RfGapTTF gap model. 
This model uses T,T',S,S' transit time factors (TTF).
It will read the TTF polynomials from the external file.
This file was created by "SNS_Cavities_Fields/rf_ttf_generator.py".
At this moment this script is not parallel.
"""

import sys
import math
import random

from bunch import Bunch

from orbit_utils import Polynomial

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF
from orbit.sns_linac import Drift

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

#------read the parameters from the external file 
fl_in = open("../SNS_Cavities_Fields/data/scl_medium_beta_rf_cav_field_t_tp_s_sp.dat","r")
lns = fl_in.readlines()
fl_in.close()

rf_freq = float(lns[0].split()[1])
(beta_min,beta_max) = (float(lns[1].split()[1]),float(lns[1].split()[2]))
n_gaps = int(lns[2].split()[1])

def split_arr_func(ln):
	res_arr = ln.split()
	val_arr = []
	for st in res_arr[1:]:
		val_arr.append(float(st)) 
	return val_arr

gap_border_points = split_arr_func(lns[3])
gap_positions = split_arr_func(lns[4])
gap_lengths = split_arr_func(lns[5])
gap_E0_amplitudes = split_arr_func(lns[6])
gap_E0L_amplitudes = split_arr_func(lns[7])

def split_polynom_coeffs(ln):
	res_arr = ln.split()
	coeff_arr = []
	count = 0
	for i in range(4,len(res_arr),3):
		val = float(res_arr[i])
		coeff_arr.append(val)
		#print "debug i=",count," val=",val
		count += 1
	poly = Polynomial(len(coeff_arr)-1)
	for i in range(len(coeff_arr)):
		poly.coefficient(i,coeff_arr[i])
	return poly

#---------------------------------------
# We set T, T', S, and S'. The T'=dT/d(kappa) and S'=dS/d(kappa).
# where kappa = 2*PI*frequency/(c*beta)
# The T' and S' are set up as separate polynomials, because
# the accuracy of calculating a derivative from the polynomial
# fitting is very low.
#---------------------------------------		
		
rf_gap_ttf_arr = []
for i_gap in range(n_gaps):
	rf_gap_ttf = RfGapTTF()
	polyT = split_polynom_coeffs(lns[8+i_gap])
	polyTp = split_polynom_coeffs(lns[9+i_gap])
	polyS = split_polynom_coeffs(lns[10+i_gap])
	polySp = split_polynom_coeffs(lns[11+i_gap])
	rf_gap_ttf.setT_TTF(polyT)
	rf_gap_ttf.setTp_TTF(polyTp)
	rf_gap_ttf.setS_TTF(polyS)
	rf_gap_ttf.setSp_TTF(polySp)
	gap_length = gap_lengths[i_gap]
	relative_amplitude = gap_E0_amplitudes[i_gap]
	rf_gap_ttf.setParameters(polyT,polyTp,polyS,polySp,beta_min,beta_max,rf_freq,gap_length,relative_amplitude)
	rf_gap_ttf_arr.append(rf_gap_ttf)

#--------directions of the cavity can be +1 or -1
directionZ = +1
if(directionZ < 0):
	for i_gap in range(len(gap_border_points)):
		gap_border_points[i_gap] = - gap_border_points[i_gap]
	for i_gap in range(len(gap_positions)):
		gap_positions[i_gap] = - gap_positions[i_gap]
	gap_border_points.reverse()
	gap_positions.reverse()
	gap_E0_amplitudes.reverse()
	gap_E0L_amplitudes.reverse()
	gap_lengths.reverse()
	rf_gap_ttf_arr.reverse()
		
print "debug ==========================================="
for i_gap in range(n_gaps):	
	rf_gap_ttf = rf_gap_ttf_arr[i_gap]
	print "=============== RF Gap index=",i_gap
	print "beta min/max= ",rf_gap_ttf.getBetaMinMax()
	print "rf_frequency= ",rf_gap_ttf.getFrequency()
	print "gap_length= ",rf_gap_ttf.getLength()
	print "relative_amplitude= ",rf_gap_ttf.getRelativeAmplitude()	
	
#----------------------------------------------------------
# RF Cavity tracking through the set of RF gaps and drifts
#----------------------------------------------------------
rf_cavity_mode = 1

drift = Drift()

def RF_Cavity_Track(b,E0,phase_dgr):
	phase = phase_dgr*math.pi/180.
	time_init = 0.
	for i_gap in range(n_gaps):
		rf_gap_ttf = rf_gap_ttf_arr[i_gap]
		drfit_1_length = gap_positions[i_gap] - gap_border_points[i_gap]
		drift_2_length = gap_border_points[i_gap+1] - gap_positions[i_gap]
		drift.setLength(drfit_1_length)
		drift.trackBunch(b)
		if(i_gap == 0): time_init = b.getSyncParticle().time()
		time_gap = b.getSyncParticle().time()
		delta_phase = 2*math.pi*(time_gap - time_init)*rf_gap_ttf.getFrequency()
		delta_phase += rf_cavity_mode*math.pi*(i_gap%2)
		rf_gap_ttf.trackBunch(b,E0,phase+delta_phase)
		#print "debug rf_pahse =",  makePhaseNear((phase+delta_phase)*180.0/math.pi,0.)," e_kin=",b.getSyncParticle().kinEnergy()
		drift.setLength(drift_2_length)
		drift.trackBunch(b)		
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
lmbd = c_light/rf_freq

#---cavity field
E0 =  20.0e+6       # average field in V/m

#---let's calculate the approximate maximal energy gain
energy_gain = 0.
for i_gap in range(n_gaps):
	rf_gap_ttf = rf_gap_ttf_arr[i_gap]
	amp = rf_gap_ttf.getRelativeAmplitude()
	length = rf_gap_ttf.getLength()
	polyT = rf_gap_ttf.getT_TTF()
	kappa = 2*math.pi*rf_gap_ttf.getFrequency()/(c_light*beta)
	ttf_t = polyT.value(kappa)
	energy_gain += E0*ttf_t*amp*length

print "Approximate maximal energy gain [MeV] = ",energy_gain/1.0e+6

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
	b1 = RF_Cavity_Track(b1,E0,phase)
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

