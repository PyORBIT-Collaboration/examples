#------------------------------------------------------
# This is an example of Function.
# It is a containers for (x,y) (or (x,y,err) ) points
#-------------------------------------------------------

import sys
import math

from orbit_utils import Function

f = Function()

def FF(x):
	return math.sin(x)

n = 10
step = 2*math.pi/n
for i in range(n):
	x = step*i+0.1*((1.0*i)/n)**2;
	y = FF(x)
	err = y*0.001
	f.add(x,y,err)

f.dump()

print "======================================="
for i in range(n):
	x = f.x(i)
	y = f.y(i)
	err = f.getYErr(x)
	print "i="," x,y,err= %10.7f %10.7f %10.7f  "%(x,y,err*1000)


