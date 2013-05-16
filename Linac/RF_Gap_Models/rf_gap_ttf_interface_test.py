#! /usr/bin/env python

"""
This script is a test for RfGapTTF gap model settings. 
This model uses T,T',S,S' transit time factors (TTF)
to calculate the 6D coordinates transformation in the RF gap. 
This test includes the memory leak test.
"""

import sys
import math
import random

from bunch import Bunch

from orbit_utils import Polynomial

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

rf_gap_ttf = RfGapTTF()
T_ttf = rf_gap_ttf.getT_TTF()
S_ttf = rf_gap_ttf.getS_TTF()
#---------------------------------------
# We set T, T', S, and S'. The T'=dT/d(cappa) and S'=dS/d(cappa).
# where cappa = 2*PI*frequency/(c*beta)
# The T' and S' are set up as separate polynomials, because
# the accuracy of calculating a derivative from the polynomial
# fitting is very low.
#---------------------------------------

polyT = Polynomial(4)
polyT.coefficient(2,2.0)
rf_gap_ttf.setT_TTF(polyT)

polyS = Polynomial(5)
polyS.coefficient(3,3.0)
rf_gap_ttf.setS_TTF(polyS)

polyTp = Polynomial(4)
polyTp.coefficient(3,2.0)
rf_gap_ttf.setT_TTF(polyTp)

polySp = Polynomial(5)
polySp.coefficient(1,3.0)
rf_gap_ttf.setS_TTF(polySp)

beta_min = 0.5
beta_max = 0.9
rf_frequency = 805.0e+6
gap_length = 0.22
relative_amplitude = 0.89

print "===========second set======================="
rf_gap_ttf.setParameters(polyT,polyTp,polyS,polySp,beta_min,beta_max,rf_frequency,gap_length,relative_amplitude)

print "===========second get======================="
T_ttf = rf_gap_ttf.getT_TTF()
S_ttf = rf_gap_ttf.getS_TTF()
Tp_ttf = rf_gap_ttf.getTp_TTF()
Sp_ttf = rf_gap_ttf.getSp_TTF()

print "========================================"
order = T_ttf.order()
for i in range(order+1):
	print "T_ttf i=",i," coef=",T_ttf.coefficient(i)

print "========================================"
order = S_ttf.order()
for i in range(order+1):
	print "S_ttf i=",i," coef=",S_ttf.coefficient(i)

print "========================================"
order = Tp_ttf.order()
for i in range(order+1):
	print "Tp_ttf i=",i," coef=",Tp_ttf.coefficient(i)

print "========================================"
order = Sp_ttf.order()
for i in range(order+1):
	print "Sp_ttf i=",i," coef=",Sp_ttf.coefficient(i)
print "========================================"

print "beta min/max= ",rf_gap_ttf.getBetaMinMax()
print "rf_frequency= ",rf_gap_ttf.getFrequency()
print "gap_length= ",rf_gap_ttf.getLength()
print "relative_amplitude= ",rf_gap_ttf.getRelativeAmplitude()

print "===========memory leak check======================"
count = 0

rf_gap_ttf = RfGapTTF()

while(1 < 2):
	count += 1

	rf_gap_ttf = RfGapTTF()
	
	polyT = Polynomial(4)
	polyT.coefficient(2,2.0)
	T_ttf = rf_gap_ttf.setT_TTF(polyT)
	
	polyS = Polynomial(5)
	polyS.coefficient(3,3.0)
	S_ttf = rf_gap_ttf.setS_TTF(polyS)
	
	polyT = Polynomial(4)
	polyT.coefficient(2,2.0)
	T_ttf = rf_gap_ttf.setTp_TTF(polyT)
	
	polyS = Polynomial(5)
	polyS.coefficient(3,3.0)
	S_ttf = rf_gap_ttf.setSp_TTF(polyS)	
	
	
	T_ttf = rf_gap_ttf.getT_TTF()
	S_ttf = rf_gap_ttf.getS_TTF()
	Tp_ttf = rf_gap_ttf.getTp_TTF()
	Sp_ttf = rf_gap_ttf.getSp_TTF()

	if(count % 100000 == 0): print "count=",count


