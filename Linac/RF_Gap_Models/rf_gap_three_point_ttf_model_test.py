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

rf_gap = RfGapThreePointTTF()

rf_frequency = 400.0*1.0e+6   # in Hz
dz = 0.01                     # in m
Em = 12.0e+6                  # in V/m
E0 = 13.2e+6                  # in V/m
Ep = 14.0e+6                  # in V/m

#kappa = 2*PI*frequency/(c*beta)
c_light = 2.99792458e+8
beta = 0.5
kappa = 2*math.pi*rf_frequency/(c_light*beta)
print "kappa =",kappa

a_param = (Ep-Em)/(2*E0*dz)
b_param = (Ep+Em-2*E0)/(2*E0*dz*dz)

def FieldFunction(z):
	return E0*(1. + a_param*z + b_param*z*z)

print "Em=",Em," E_Func(-dz)= ",FieldFunction(-dz)
print "E0=",E0," E_Func(  0)= ",FieldFunction(0.)
print "Ep=",Ep," E_Func(+dz)= ",FieldFunction(+dz)

E0L = E0*(2*dz+(2./3.)*b_param*dz*dz*dz)


print "=========================================="
print "Stop."

