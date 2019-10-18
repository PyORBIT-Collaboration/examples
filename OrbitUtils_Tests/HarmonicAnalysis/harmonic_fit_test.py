#-------------------------------------------------------------------------------
#.
#-------------------------------------------------------------------------------

import sys
import math

from orbit_utils import Function
from orbit_utils import HarmonicData

import numpy as np
import scipy
from scipy.optimize import minimize

f = Function()

for ind in range(0,360,30):
	x = 1.0*ind
	y = 2.0*math.cos((math.pi/180.)*(x+ 25.))+ 4.0*math.cos((math.pi/180.)*(4*x+ 35.)) + 0.5
	y_err = 0.01*abs(y)
	f.add(x,y,y_err)

order = 4

harmonic_data = HarmonicData(order,f)

x_arr = [0.8,2.1,25.2,0.,0.,0.,0.,4.3,35.4]
for x_ind in range(len(x_arr)):
	harmonic_data.parameter(x_ind,x_arr[x_ind])

harmonic_data.parameter(0,0.5)
harmonic_data.parameter(1,2.1)
harmonic_data.parameter(2,25.2)
harmonic_data.parameter(7,4.3)
harmonic_data.parameter(8,35.4)

class HarmonicFitFunction:
	def __init__(self,harmonic_data):
		self.harmonic_data = harmonic_data
		
	def getDiff2(self,x_arr):
		for x_ind in range(len(x_arr)):
			self.harmonic_data.parameter(x_ind,x_arr[x_ind])
		return harmonic_data.sumDiff2()
		
print "==================================="
harm_fit_func = HarmonicFitFunction(harmonic_data)
x0 = np.array(x_arr)
res = minimize(harm_fit_func.getDiff2, x0, method='nelder-mead', options={'maxiter': 1500, 'xtol': 1e-8, 'disp': True})
print res.x
print "==================================="
"""
x0 = np.array(x_arr)
res = scipy.optimize.fmin(harm_fit_func.getDiff2, x0, xtol= 1e-8, maxiter=1500)
print res
print "==================================="
"""

for ind in range(harmonic_data.dataSize()):
	x = harmonic_data.valueX(ind)
	y = harmonic_data.valueY(ind)
	y_err = harmonic_data.valueErr(ind)
	y_fit = harmonic_data.fitValueY(x)
	print "ind= %3d "%ind, " (x,y+-y_err,y_fit) = ( %8.1f , %+8.5f +- %8.5f ,%8.5f)"%(x,y,y_err,y_fit)
	
	