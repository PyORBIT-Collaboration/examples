#-----------------------------------------------------
#Grid2D gradient test
#-----------------------------------------------------

import sys
import math

from spacecharge import Grid1D


print "Start."

sizeZ = 100
zMin =  1.0
zMax = +4.0

grid1D = Grid1D(sizeZ,zMin,zMax)

def Func(z):
	return 1.0/(z*z)

def FuncGrad(z):
	return -2.0/(z*z*z)


for iz in range(sizeZ):
	z = grid1D.getGridZ(iz)
	val = Func(z)
	grid1D.setValue(val,iz)

diff_max = 0.
max_point = z
for iz in range(sizeZ):
	z = grid1D.getGridZ(iz)
	grad0 = FuncGrad(z)
	grad = grid1D.calcGradient(z)
	diff = math.sqrt((grad - grad0)**2)
	if(diff > diff_max):
		diff_max = diff
		max_point = z
		
		
print "max diff =",diff_max
print "max point =",max_point
z = max_point
print "Func at max point =",Func(z)
print "Grad at max point =",FuncGrad(z)

print "Stop."

