#-----------------------------------------------------
#Creates Grid1D Grid2D for charge density and Poisson Solver
#for a charged string (2D case) 
#Copmares with the exact result.
#-----------------------------------------------------
#The one very useful property of the FFT Poisson solver - scalability 
#you can change the absolute value of steps X and Y and you will get 
#the right results as soon you keep stepX/stepY constant.
#

import sys
import math

import orbit_mpi

from spacecharge import Grid1D
from spacecharge import Grid2D
from bunch import Bunch
from spacecharge import PoissonSolverFFT2D
from spacecharge import SpaceChargeCalc2p5D
from spacecharge import Boundary2D

print "Start."

sizeX = 10
sizeY = 10
xMin = -5.0
xMax = +5.0
yMin = -5.0
yMax = +5.0

sizeZ = 10
zMin = -5.0
zMax = +5.0

nParts = 5
charge = 1.0

nBoundaryPoints = 50
N_FreeSpaceModes = 20
R_Boundary = 1.0

slicelen = 1.0

b = Bunch()
b.readBunch("Bm_Parts_ini_0")
print "bunchSize = ",b.getSize()
syncPart = b.getSyncParticle()
energy = 1.0
syncPart.kinEnergy(energy)

calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ,1)
print "without boundary."
calc2p5d.trackBunch(b,slicelen)
b.dumpBunch("test_result_0.out")

boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes)

for i in xrange(nBoundaryPoints):
	x = R_Boundary*math.cos((2.0*math.pi/(nBoundaryPoints-1))*i)
	y = R_Boundary*math.sin((2.0*math.pi/(nBoundaryPoints-1))*i)
	boundary.setBoundaryPoint(i,x,y)
boundary.initialize()

boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",xMax-xMin)
print "shape name=",boundary.getShapeName()

calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ,1)
print "with boundary."
calc2p5d.trackBunch(b,slicelen*2,boundary)
b.dumpBunch("test_result_1.out")

print "Stop."

