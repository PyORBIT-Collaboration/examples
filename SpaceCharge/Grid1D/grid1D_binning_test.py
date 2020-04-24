#-----------------------------------------------------
# Grid2D values binning test
# User can change number of grid points (sizeZ) and the type of binning 
# (grid1D.binValue or grid1D.binValueSmoothed)
# The test should be performed for sizeZ = 1,2,3,  20
# You expect to see uniform distribution
#-----------------------------------------------------

import sys
import math

from spacecharge import Grid1D


print "Start."

sizeZ = 3
zMin =  1.0
zMax = +4.0

extent = 1.0

grid1D = Grid1D(sizeZ,zMin,zMax)


def FuncTest(z):
	y = 2.0 + 0.5*(z-1.) + z**2
	return y

#-------------------------------------------
# Bin the uniformly distributed values along z
#-------------------------------------------
nPlotPoints = 10000
x_polt_min = zMin - extent
x_polt_max = zMax + extent
step = (x_polt_max - x_polt_min)/(nPlotPoints - 1)

value = 1.0

for ind_p in range(nPlotPoints):
	x = x_polt_min + ind_p*step
	#grid1D.binValue(value,x)
	grid1D.binValueSmoothed(value,x)

#---------------------------------
# Plot the density
#---------------------------------

x_arr = []
y_arr = []
for ind in range(sizeZ):
	x = grid1D.getGridZ(ind)
	y = grid1D.getValueOnGrid(ind)
	print "debug z=",x," rho=",y
	x_arr.append(x)
	y_arr.append(y)
		
		
#---------------------------
# Plot part
#---------------------------
import matplotlib.pyplot as plt	
	
plt.plot(x_arr,y_arr)
plt.ylabel('rho(x)')
plt.xlabel('x')

plt.show()