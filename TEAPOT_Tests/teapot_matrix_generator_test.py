from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot_base import MatrixGenerator
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit_utils import Matrix

print "Start."

b = Bunch()

syncPart = b.getSyncParticle()
energy = 1.0                          #energy in GeV
syncPart.kinEnergy(energy)

lattice = TEAPOT_Lattice("test_lattice")

elem1 = teapot.DriftTEAPOT("drift1")
elem2 = teapot.QuadTEAPOT("quad")
elem3 = teapot.DriftTEAPOT("drift1")

#lattice.addNode(elem1)
lattice.addNode(elem2)
#lattice.addNode(elem3)

#-----------------------------
# Set TEAPOT nodes parameters
#-----------------------------
elem1.setLength(1.0)
elem2.setLength(1.0)
elem3.setLength(1.0)

elem2.setnParts(5)
elem2.addParam("kq",-0.5)

lattice.initialize()
print "lattice length=",lattice.getLength()

transp_matrix = Matrix(6,6)

matrixGenerator = MatrixGenerator()

matrixGenerator.initBunch(b)
lattice.trackBunch(b)
matrixGenerator.calculateMatrix(b,transp_matrix)


def printM(m):
	print "----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("(%1d,%1d)= % 6.5e "%(i,j,m.get(i,j))),
		print ""	
		
printM(transp_matrix)

print "Stop."

