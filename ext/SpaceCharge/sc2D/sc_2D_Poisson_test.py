#-----------------------------------------------------
#Creates Grid2D for charge density and potential
#Boundary2D instance for solving Poisson problem
#for a charged string (2D case) 
#Copmares with the exact result.
#-----------------------------------------------------
import sys
import math

import orbit_mpi

from spacecharge import Grid2D
from spacecharge import Boundary2D

print "Start."

nBinX = 200
nBinY = 200
xSize = 10.
ySize = 10.
nBoundaryPoints = 100
BoundaryShape = "Circle"
N_FreeSpaceModes = 20

boundary = Boundary2D(nBinX,nBinY,xSize,ySize,nBoundaryPoints,BoundaryShape,N_FreeSpaceModes)

gridRho = Grid2D(boundary)
gridPhi = Grid2D(boundary)

chrage_pos_x = 2.5
chrage_pos_y = 0.0
charge = 1.0
gridRho.binValue(charge,chrage_pos_x,chrage_pos_y)

gridRho.findPotential(gridPhi)

r_test = 4.0
n_angle_steps = 10
angle_step = 360./(n_angle_steps - 1)
print "  i    x       y         phi         phi_theory    ratio phi/theory  "
for i in xrange(n_angle_steps):
	angle = math.pi*i*angle_step/180.
	x = r_test*math.cos(angle)
	y = r_test*math.sin(angle)
	phi = gridPhi.getValue(x,y)
	dist = (chrage_pos_x - x)*(chrage_pos_x - x) + (chrage_pos_y - y)*(chrage_pos_y - y)
	dist = math.sqrt(dist)
	phi_th = -math.log(dist)
	print "",i," %7.4f  %7.4f  %12.5g  %12.5g  %12.5g  "%(x,y,phi,phi_th,phi/phi_th) 

print "Stop."

