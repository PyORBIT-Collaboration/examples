#-----------------------------------------------------
# Creates Space Charge Calculator
# test Boundary
#-----------------------------------------------------
import sys
import math
import random

import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5D
from spacecharge import Boundary2D

print "Start."

# Make a SC solver
sizeX = 64
sizeY = 64
sizeZ = 20
long_avg_n = 3
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)
calc2p5d.setLongAveragingPointsN(long_avg_n)

# Make bunch
b = Bunch()
bunch_radius = 0.005
bunch_length = 500e-9*3e8*0.87502565
nParts = 100000

#Use random radius
for ip in range(nParts):
	r = bunch_radius*math.sqrt(random.random())
	phi = 2*math.pi*random.random()
	x = r*math.sin(phi)
	y = r*math.cos(phi)
	z = bunch_length*0.5*(1.0 - 2*random.random())
	b.addParticle(x,0.,y,0.,z,0.)

macroSize = 1.56e+13
energy = 0.08
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(macroSize/nParticlesGlobal)
b.getSyncParticle().kinEnergy(energy)	

# Make a boundary
nBoundaryPoints = 100
N_FreeSpaceModes = 20
R_Boundary = 0.0049
boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes)

for i in xrange(nBoundaryPoints):
	x = R_Boundary*math.cos((2.0*math.pi/(nBoundaryPoints-1))*i)
	y = R_Boundary*math.sin((2.0*math.pi/(nBoundaryPoints-1))*i)
	boundary.setBoundaryPoint(i,x,y)
boundary.initialize()
b_maxx = boundary.getMaxX()
b_minx = boundary.getMinX()
b_maxy = boundary.getMaxY()
b_miny = boundary.getMinY()
print "MaxX=",b_maxx," MinX=",b_minx," MaxY=",b_maxy," MinY=",b_miny

# Shape boundary test
#Choose from Circle/Ellipse/Rectangle
shapeboundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",0.007,0.007)
print "shape name=",boundary.getShapeName()

# Set SC node parameters
pipe_radius = 0.010
slice_length = 0.1

# Track and analysis
#analysis
rhoGrid = calc2p5d.getRhoGrid()
phiGrid = calc2p5d.getPhiGrid()
#x = bunch_radius/2.0
x = R_Boundary/2.0
y = 0.
r = math.sqrt(x*x+y*y)

rho_theory = macroSize*4.0/(math.pi*rhoGrid.getSizeX()*rhoGrid.getSizeY())
phi_theory = macroSize*r**2/(bunch_radius**2)
grad_theory = macroSize*2*r/(bunch_radius**2)

#without boundary
print "without boundary."
calc2p5d.trackBunch(b,slice_length,pipe_radius)

rho = rhoGrid.getValue(x,y)
phi = 2*(phiGrid.getValue(x,y) - phiGrid.getValue(0.,0.))
(ex,ey) = phiGrid.calcGradient(x,y)
grad = 2*math.sqrt(ex*ex+ey*ey)

print "r=",r," rho  = %12.5g "%rho,"  rho_theory = %12.5g "%rho_theory
print "r=",r," phi  = %12.5g "%phi,"  phi_theory = %12.5g "%phi_theory
print "r=",r," grad = %12.5g "%grad," grad_theory = %12.5g "%grad_theory

#with boundary
print "with boundary."
calc2p5d.trackBunch(b,slice_length,pipe_radius,boundary)
#calc2p5d.trackBunch(b,slice_length,pipe_radius,shapeboundary)

rho = rhoGrid.getValue(x,y)
phi = 2*(phiGrid.getValue(x,y) - phiGrid.getValue(0.,0.))
(ex,ey) = phiGrid.calcGradient(x,y)
grad = 2*math.sqrt(ex*ex+ey*ey)

print "r=",r," rho  = %12.5g "%rho,"  rho_theory = %12.5g "%rho_theory
print "r=",r," phi  = %12.5g "%phi,"  phi_theory = %12.5g "%phi_theory
print "r=",r," grad = %12.5g "%grad," grad_theory = %12.5g "%grad_theory

xyp_arr = []
for ip in range(b.getSize()):
	xp = b.xp(ip)
	yp = b.yp(ip)
	xyp_arr.append((xp,yp))
	
# Plot
import Gnuplot
gXYp = Gnuplot.Gnuplot()
gXYp.title('momentum')
gXYp.xlabel('xp')
gXYp.ylabel('yp')
gXYp('set data style line')
gXYp.plot(xyp_arr)

raw_input('Please press return to continue...\n')
