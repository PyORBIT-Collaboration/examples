#! /usr/bin/env python

"""
This script is a test for RfGapTTF gap model. This model uses T,T',S,S' transit time factors (TTF)
to calculate the 6D coordinates transformation in the RF gap. 
"""

import sys
import math
import random

from bunch import Bunch

from orbit_utils import Polynomial

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

rf_gap_ttf = RfGapTTF()

#---------------------------------------
# We have to set only T and S. The T'=dT/d(cappa) and S'=dS/d(cappa)
# will be calculated internaly.
# cappa = 2*PI*frequency/(c*beta)
#---------------------------------------

polyT = Polynomial(4)
polyT.coefficient(2,2.0)
#rf_gap_ttf.setT_TTF(polyT)

polyS = Polynomial(5)
polyS.coefficient(3,3.0)
#rf_gap_ttf.setS_TTF(polyS)

beta_min = 0.5
beta_max = 0.9
rf_frequency = 805.0e+6
gap_length = 0.22
relative_amplitude = 0.89

rf_gap_ttf.setParameters(polyT,polyS,beta_min,beta_max,rf_frequency,gap_length,relative_amplitude)

T_ttf = rf_gap_ttf.getT_TTF()
S_ttf = rf_gap_ttf.getS_TTF()

print "========================================"
order = T_ttf.order()
for i in range(order+1):
	print "T_ttf i=",i," coef=",T_ttf.coefficient(i)

print "========================================"
order = S_ttf.order()
for i in range(order+1):
	print "S_ttf i=",i," coef=",S_ttf.coefficient(i)
print "========================================"

print "beta min/max= ",rf_gap_ttf.getBetaMinMax()
print "rf_frequency= ",rf_gap_ttf.getFrequency()
print "gap_length= ",rf_gap_ttf.getLength()
print "relative_amplitude= ",rf_gap_ttf.getRelativeAmplitude()

"""
count = 0
while(1 < 2):
	count += 1
	T_ttf = rf_gap_ttf.getT_TTF()
	S_ttf = rf_gap_ttf.getS_TTF()
	if(count % 1000000 == 0): print "count=",count
"""
	
