#------------------------------------------------------
# This is an example of Function and derivative 
# functionality of this class
#-------------------------------------------------------

import sys
import math
import random

from orbit_utils import Function

def FF(x):
	return math.sin(x)

def FFP(x):
	return math.cos(x)
	
f = Function()

eps = 1.0e-5
n = 100
step = 2*math.pi/n
for i in range(n):
	x = step*i+ random.uniform(-eps*step,eps*step);
	y = FF(x)
	f.add(x,y)

f.dump()

# let's check  Const step tolerance setting procedure
f.setStepEps(4*eps)
print "Const step tolerance =",f.getStepEps()
res = f.setConstStep(1)
print "Function has a const step on x-variable res=",res

print "==========================================="
print " x   abs(y-y_fit)  abs(dy/dx - dy/dx_fit) "

n = 20
step = 0.8*(f.getMaxX() - f.getMinX())/(n-1)
y_dev_max = 0.
yp_dev_max = 0.
for j in range(n):
	x = f.getMinX() + 0.0001 + j*step
	y = f.getY(x)
	yp = f.getYP(x)
	y_th = FF(x)
	yp_th = FFP(x)
	st = "%6.5f "%x + " %10.5f  %10.5f "%(abs(y-y_th),abs(yp-yp_th))
	print st



