##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and calculates the set of 
# transport matrices for each element. At the end, it will print
# the one turn matrix for the ring.
# The methods used:
#  1. An action in a action container and a track method of 
#     TEAPOT nodes
#  2. A transport matrix generation for the whole TEAPOT lattice
#
#  3. A trackBunch method for each TEAPOT node
# 
#  The 1. method is the safest, because it does not track the bunch
#  through possible Space Charge nodes, RF Nodes with nonzero voltage
#  etc. If you know that you do not have them, they shoild give us the
#  same results.
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import BaseTEAPOT
from orbit.teapot import RingRFTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
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

#---Transformation of the matrix to the MAD variables
#---The momentum and beta are different 
#---MAD uses -ct instead of z, and dE/p instead of dE
def transormToMAD(m,p,beta):
	m.set(0,5,m.get(0,5)*p)
	m.set(1,5,m.get(1,5)*p)
	m.set(4,5,-m.get(4,5)*p/beta)
	m.set(4,0,-m.get(4,0)/beta)
	m.set(4,1,-m.get(4,1)/beta)
	
print "Start."

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("../MAD/LATTICE","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()

#=====track action approach============
class MatrixTracker:
	def __init__(self, Ekin = 1.0):
		self.b = Bunch()
		self.b.getSyncParticle().kinEnergy(Ekin)
		self.matrixGenerator = MatrixGenerator()
		self.matrixArr = []
		self.index = 0
		
	def _twissAction(self,paramsDict):
		node = paramsDict["node"]
		bunch = paramsDict["bunch"]
		s_start = paramsDict["position"]
		s_stop = s_start + node.getLength(node.getActivePartIndex())
		paramsDict["position"] = s_stop
		if(isinstance(node, BaseTEAPOT) == True and isinstance(node,RingRFTEAPOT) == False): 
			#print "i=",self.index," ind=",node.getActivePartIndex()," name=",node.getName(),
			#print " type=",node.getType()," L=",node.getLength(node.getActivePartIndex())
			self.matrixGenerator.initBunch(bunch)
			node.track(paramsDict)
			transp_matrix = Matrix(6,6)
			self.matrixGenerator.calculateMatrix(bunch,transp_matrix)
			self.matrixArr.append((s_start,s_stop,transp_matrix))
			m = transp_matrix
			#print "detX=", (m.get(0,0)*m.get(1,1) - m.get(0,1)*m.get(1,0))
			#printM(transp_matrix)
			self.index = self.index + 1
			#if(self.index > 1000): sys.exit(1)
				
	def getTransportMatrixArray(self,teapot_lattice):
		accContainer = AccActionsContainer()
		accContainer.addAction(self._twissAction,AccActionsContainer.BODY)
		paramsDict = {}
		paramsDict["bunch"] = self.b
		paramsDict["position"] = 0.
		self.matrixArr = []
		self.index = 0
		teapot_lattice.trackActions(accContainer,paramsDict)
		return self.matrixArr
		
	def getOneTurnMatrix(self,teapot_lattice):
		self.matrixArr = self.getTransportMatrixArray(teapot_lattice)
		final_matrix = Matrix(6,6)
		final_matrix.unit()
		for (s0,s1,matrix) in self.matrixArr:
			final_matrix = matrix.mult(final_matrix)
		return final_matrix

matrixTracker = MatrixTracker(1.0)
momentum = matrixTracker.b.getSyncParticle().momentum()
beta = matrixTracker.b.getSyncParticle().beta()

matrixArr = matrixTracker.getTransportMatrixArray(teapot_latt)
print "n marix =",len(matrixArr)
final_matrix = matrixTracker.getOneTurnMatrix(teapot_latt)
transormToMAD(final_matrix,momentum,beta)
printM(final_matrix)

#=====trackBunch method for the whole lattice approach============
print "============= Matrix for the whole lattice =========="
transp_matrix = Matrix(6,6)
matrixGenerator = MatrixGenerator()
matrixGenerator.initBunch(matrixTracker.b)
teapot_latt.trackBunch(matrixTracker.b)
matrixGenerator.calculateMatrix(matrixTracker.b,transp_matrix)
transormToMAD(transp_matrix,momentum,beta)
printM(transp_matrix)

#=====trackBunch method for each TEAPOT element approach============
print "Total number of TEAPOT Base Nodes=",len(teapot_latt.getNodes())
printNames = []
res_matrix = Matrix(6,6)
res_matrix.unit()
n_max = len(teapot_latt.getNodes())
#n_max = 74
length = 0.
for i in range(n_max):
	node = teapot_latt.getNodes()[i]
	#print "i=",i," name =",node.getName()," L=",length	
	m = Matrix(6,6)
	matrixGenerator.initBunch(matrixTracker.b)
	node.trackBunch(matrixTracker.b)
	matrixGenerator.calculateMatrix(matrixTracker.b,m)
	#m_prt0 = res_matrix.copy()
	#transormToMAD(m_prt0,momentum,beta)
	#printM(m_prt0)	
	res_matrix = m.mult(res_matrix)
	#transormToMAD(m,momentum,beta)
	#printM(m)		
	#m_prt = res_matrix.copy()	
	#transormToMAD(m_prt,momentum,beta)
	#printM(m_prt)	
	length = length + node.getLength()
	
print "====res for TEAPOT n elements=",n_max," L=",length
transormToMAD(res_matrix,momentum,beta)
printM(res_matrix)


m = transp_matrix
det_x = m.get(0,0)*m.get(1,1) - m.get(0,1)*m.get(1,0)
print "determinant det(Mx)=",det_x 
det_y = m.get(2,2)*m.get(3,3) - m.get(2,3)*m.get(3,2)
print "determinant det(My)=",det_y 

cos_phi_x = (m.get(0,0)+m.get(1,1))/2.0
cos_phi_y = (m.get(2,2)+m.get(3,3))/2.0

if(math.fabs(cos_phi_x) >= 1.0 or  math.fabs(cos_phi_x) >= 1.0):
	print "Unstable MAD lattice!"
	print "Stop."
	sys.exit(1)

nux = math.acos(cos_phi_x)/(2*math.pi)
nuy = math.acos(cos_phi_y)/(2*math.pi)


print "fractional tune nux=",nux
print "fractional tune nuy=",nuy

print "Stop."


