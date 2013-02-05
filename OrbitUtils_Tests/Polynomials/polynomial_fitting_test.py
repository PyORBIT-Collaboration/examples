#----------------------------------------------
# This is an example of a polymonial fitting 
#----------------------------------------------

import sys
import math

from orbit_utils import Function
from orbit_utils import SplineCH

from orbit.utils.fitting import PolynomialFit

f = Function()
for i in range(15):
	x = 0.1*i
	y = -1+1.0*x**2+2.0*x**3 + 0.01*x**4  
	f.add(x,y)
	
print "polynomial for fitting = y = -1+1.0*x**2+2.0*x**3 + 0.01*x**4"
print "================Function===================================="
poly = PolynomialFit(4)
poly.fitFunction(f)
[coef_arr,err_arr] = poly.getCoefficientsAndErr()
for i in range(len(coef_arr)):
	print "i=",i," coef = %8.5f +- %8.5f"%(coef_arr[i],err_arr[i])
	
print "================SplineCH===================================="
spline = SplineCH()
spline.compile(f)
poly.fitSpline(f)	
[coef_arr,err_arr] = poly.getCoefficientsAndErr()
for i in range(len(coef_arr)):
	print "i=",i," coef = %8.5f +- %8.5f"%(coef_arr[i],err_arr[i])
