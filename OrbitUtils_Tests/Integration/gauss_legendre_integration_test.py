#-------------------------------------------------------------------------
# This is an example of usage of GaussLegendreIntegrator class.
# This class implements the Gauss Legendre method of integration.
# The integration object should be a Function or SplineCH instance. 
#--------------------------------------------------------------------------

import sys
import math


from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator

integrator = GaussLegendreIntegrator(4)
print "Number of integral points =",integrator.getnPoints()
print "The integral of sin(x) from 0. to pi/2 is 1."
#------------------------------------------
# integral of sin(x) from 0. to pi/2 is 1.
#------------------------------------------

x_min = 0.
x_max = math.pi/2 
integrator.setLimits(x_min,x_max)


f = Function()
n = 10
for i in range(n):
	x = x_min + i*(x_max - x_min)/(n-1)
	f.add(x,math.sin(x))
	

res = integrator.integral(f)
print "For Function integral =",res," error=",math.fabs(res-1.0)

spline = SplineCH()
spline.compile(f)
res = integrator.integral(spline)
print "For SplineCH integral =",res," error=",math.fabs(res-1.0)
print "=============== Integration points and weights ================"
point_weight_arr = integrator.getPointsAndWeights()
for (x,w) in point_weight_arr:
	print "x = %8.5f   w = %12.5g "%(x,w) 
print "Done."



