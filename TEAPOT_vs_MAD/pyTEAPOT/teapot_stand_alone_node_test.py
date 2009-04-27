##############################################################
# This script shows how to generate transport matrix for 
# stand alone TEAPOT elements. There is no need to make a  
# TEAPOT lattice
##############################################################
import math
import sys

from orbit.teapot import BaseTEAPOT, DriftTEAPOT, QuadTEAPOT
from orbit.teapot_base import MatrixGenerator
from orbit_utils import Matrix
from bunch import Bunch

#---PRINT Function for Matrix
def printM(m):
	print "----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("(%1d,%1d)=% 6.5e "%(i,j,m.get(i,j))),
		print ""	
		

b = Bunch()
Ekin = 1.0 # in GeV
b.getSyncParticle().kinEnergy(Ekin)

#define TEAPOT drift
node0 = DriftTEAPOT("drift")
node0.setLength(1.0)

#define TEAPOT quad
node1 = QuadTEAPOT("quad")
node1.setLength(1.0)
node1.addParam("kq",0.5)
	
matrixGenerator = MatrixGenerator()

#========matrix for drift =====
m = Matrix(6,6)
matrixGenerator.initBunch(b)
node0.trackBunch(b)
matrixGenerator.calculateMatrix(b,m)	
print "drift matrix L=",node0.getLength()
printM(m)


#========matrix for quad =====
matrixGenerator.initBunch(b)
node1.trackBunch(b)
matrixGenerator.calculateMatrix(b,m)	
print "quad matrix L=",node1.getLength()," kq=",node1.getParam("kq")
printM(m)

