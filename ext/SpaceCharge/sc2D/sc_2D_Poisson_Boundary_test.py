#-----------------------------------------------------
#Creates Grid2D for charge density and potential
#Boundary2D instance for solving Poisson problem
#for a charged string (2D case) with a boundary induced potential
#Compares with the exact result.
# phi(r) = ln(abs(r-a)/abs(r-a')) where r,a,a' are vectors
#           abs(a') = R^2/abs(a)  a' is parallel to a
#           a - specifies the position of the charge
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
chrage_pos_a = math.sqrt(chrage_pos_x*chrage_pos_x + chrage_pos_y*chrage_pos_y)
charge = 1.0
gridRho.binValue(charge,chrage_pos_x,chrage_pos_y)

R = boundary.getBoundarySizeX()/2.0
a_prime = R*R/chrage_pos_a
a_prime_x = a_prime*chrage_pos_x/chrage_pos_a
a_prime_y = a_prime*chrage_pos_y/chrage_pos_a

gridRho.findPotential(gridPhi)
gridPhi.addBoundaryPotential()

#-----potential delta-------------------------------
#our potential on the wall is zero, but the exact 
#solution is not zero. There is a constant difference
#between these potentials
#---------------------------------------------------
dist = (chrage_pos_x - R)*(chrage_pos_x - R) + chrage_pos_y*chrage_pos_y
dist = math.sqrt(dist)
dist_prime = (a_prime_x - R)*(a_prime_x - R) + a_prime_y*a_prime_y
dist_prime = math.sqrt(dist_prime)
phi_delta = -math.log(dist/dist_prime)	


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
	dist_prime = (a_prime_x - x)*(a_prime_x - x) + (a_prime_y - y)*(a_prime_y - y)
	dist_prime = math.sqrt(dist_prime)
	phi_th = -math.log(dist/dist_prime)	- phi_delta
	ratio = 0.
	if(phi_th != 0.): ratio = phi/phi_th
	print "",i," %7.4f  %7.4f  %12.5g  %12.5g  %12.5g  "%(x,y,phi,phi_th,ratio) 

print "Stop."

