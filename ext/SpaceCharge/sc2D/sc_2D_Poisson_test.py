#-----------------------------------------------------
#Creates Grid2D for charge density and potential
#Boundary2D instance for solving Poisson problem
#for a charged string (2D case) 
#-----------------------------------------------------
import sys

from spacecharge import Grid2D
from spacecharge import Boundary2D

print "Start."

nBinX = 100
nBinY = 200
xSize = 10.
ySize = 10.
nBoundaryPoints = 100
BoundaryShape = "Circle"
N_FreeSpaceModes = 20

boundary = Boundary2D(nBinX,nBinY,xSize,ySize,nBoundaryPoints,BoundaryShape,N_FreeSpaceModes)

gridRho = Grid2D(boundary)
gridPhi = Grid2D(boundary)

gridRho.binValue(1.0,2000.,0.)


print "Stop."

sys.exit(1)
