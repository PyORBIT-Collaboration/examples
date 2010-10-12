#-----------------------------------------------------
#Grid2D gradient test
#-----------------------------------------------------

import sys
import math

from spacecharge import Grid2D


print "Start."

sizeX = 100
sizeY = 100
xMin =  2.0
xMax = +4.0
yMin =  1.0
yMax = +4.0

grid2D = Grid2D(sizeX,sizeY,xMin,xMax,yMin,yMax)

def Func(x,y):
	return 1.0/(x*x+y*y)

def FuncGrad(x,y):
	return (-2*x/(x*x+y*y)**2,-2*y/(x*x+y*y)**2)


for ix in range(sizeX):
	x = grid2D.getGridX(ix)
	for iy in range(sizeY):
		y = grid2D.getGridY(iy)
		val = Func(x,y)
		grid2D.setValue(val,ix,iy)

diff_max = 0.
max_point = (0.,0.)
for ix in range(sizeX):
	x = grid2D.getGridX(ix) + 0.5*(grid2D.getMaxX() - grid2D.getMinX())/grid2D.getSizeX()
	for iy in range(sizeY):
		y = grid2D.getGridY(iy) + 0.5*(grid2D.getMaxY() - grid2D.getMinY())/grid2D.getSizeY()
		(gradX0,gradY0) = FuncGrad(x,y)
		(gradX,gradY) = grid2D.calcGradient(x,y)
		diff = math.sqrt((gradX - gradX0)**2 + (gradY - gradY0)**2)
		if(diff > diff_max):
			diff_max = diff
			max_point = (x,y)
		
		
print "max diff =",diff_max
print "max point =",max_point
(x,y) = max_point
print "Func at max point =",Func(x,y)
print "Grad at max point =",FuncGrad(x,y)

print "Stop."

