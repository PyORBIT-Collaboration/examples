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

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator

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

def T_TTF(kappa, n = 32):
	f = Function()
	z_step = 2.0*dz/(n-1)
	for i in range(n):
		z = -dz + i*z_step
		f.add(z,FieldFunction(z)*math.cos(kappa*z))
	spline = SplineCH()
	spline.compile(f)	
	integrator = GaussLegendreIntegrator(n)	
	integrator.setLimits(-dz,+dz)
	return integrator.integral(spline)/E0L

def Tp_TTF(kappa, n = 32):
	f = Function()
	z_step = 2.0*dz/(n-1)
	for i in range(n):
		z = -dz + i*z_step
		f.add(z,-z*FieldFunction(z)*math.sin(kappa*z))
	spline = SplineCH()
	spline.compile(f)	
	integrator = GaussLegendreIntegrator(n)	
	integrator.setLimits(-dz,+dz)
	return integrator.integral(spline)/E0L

def S_TTF(kappa, n = 32):
	f = Function()
	z_step = 2.0*dz/(n-1)
	for i in range(n):
		z = -dz + i*z_step
		f.add(z,FieldFunction(z)*math.sin(kappa*z))
	spline = SplineCH()
	spline.compile(f)	
	integrator = GaussLegendreIntegrator(n)	
	integrator.setLimits(-dz,+dz)
	return integrator.integral(spline)/E0L

def Sp_TTF(kappa, n = 32):
	f = Function()
	z_step = 2.0*dz/(n-1)
	for i in range(n):
		z = -dz + i*z_step
		f.add(z,z*FieldFunction(z)*math.cos(kappa*z))
	spline = SplineCH()
	spline.compile(f)	
	integrator = GaussLegendreIntegrator(n)	
	integrator.setLimits(-dz,+dz)
	return integrator.integral(spline)/E0L

print "==============================================="
res_T = T_TTF(kappa,256)
res_T_test = rf_gap.getT_TTF(dz,a_param,b_param,kappa)
print "T = %16.9g "%res_T," cpp_model= %16.9g"%res_T_test

res_Tp = Tp_TTF(kappa,256)
res_Tp_test = rf_gap.getTp_TTF(dz,a_param,b_param,kappa)
print "Tp= %16.9g "%res_Tp," cpp_model= %16.9g"%res_Tp_test

res_S = S_TTF(kappa,256)
res_S_test = rf_gap.getS_TTF(dz,a_param,b_param,kappa)
print "S = %16.9g "%res_S," cpp_model= %16.9g"%res_S_test

res_Sp = Sp_TTF(kappa,256)
res_Sp_test = rf_gap.getSp_TTF(dz,a_param,b_param,kappa)
print "Sp= %16.9g "%res_Sp," cpp_model= %16.9g"%res_Sp_test

print "=========================================="
print "Stop."

