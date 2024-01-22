#-----------------------------------------------------
# This example tracks 4 particles through a quadrupole.
# The result is compared with an analytical formula.
#-----------------------------------------------------
import sys
import math
import random


from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot

from bunch import Bunch

def getTransportMatrix(kq, charge, length):
	"""
	Return the quad 4x4 transport matrix for (x,x',y,y').
	kq = (dB/dr)*(1/Brho)
	where Brho = 3.33564*sqrt(Tkin*(Tkin + 2*mass)) 
	Tkin and mass in GeV
	dB/dr in [T/m]
	charge - charge in electron charge. It can be >0 or <0
	length - length of quad in meters
	"""
	kqc = kq*charge
	sqrt_kq  = math.sqrt(abs(kqc))
	kqlength = sqrt_kq * length
	cs = math.cos(kqlength)
	ss = math.sin(kqlength)
	ch = math.cosh(kqlength)
	sh = math.sinh(kqlength)
	#---- transport matrix
	m = []
	m.append([1.,0.,0.,0.])
	m.append([0.,1.,0.,0.])
	m.append([0.,0.,1.,0.])
	m.append([0.,0.,0.,1.])
	if(kqc == 0.): return m
	if(kqc > 0.):
		m[0][0] = cs
		m[0][1] = ss / sqrt_kq
		m[1][0] = -ss * sqrt_kq
		m[1][1] = cs
		m[2][2] = ch
		m[2][3] = sh / sqrt_kq
		m[3][2] = sh * sqrt_kq
		m[3][3] = ch
	else:
		m[0][0] = ch
		m[0][1] = sh / sqrt_kq
		m[1][0] = sh * sqrt_kq
		m[1][1] = ch
		m[2][2] = cs
		m[2][3] = ss / sqrt_kq
		m[3][2] = -ss * sqrt_kq
		m[3][3] = cs
	return m

def trackBunchThroughTransportMatrix(bunch,matrix):
	bunch_out = Bunch()
	bunch.copyEmptyBunchTo(bunch_out)
	for ind in range(bunch.getSize()):
		(x,xp,y,yp) = (bunch.x(ind),bunch.xp(ind),bunch.y(ind),bunch.yp(ind))
		vector_in = [x,xp,y,yp]
		vector_out = [0.,0.,0.,0]
		for i in range(4):
			vector_out[i] = 0.
			for j in range(4):
				vector_out[i] += matrix[i][j]*vector_in[j]
		[x_out,xp_out,y_out,yp_out] = vector_out
		bunch_out.addParticle(x_out,xp_out,y_out,yp_out,bunch.z(ind),bunch.dE(ind))
	return bunch_out

#=====================================================
#   START OF THE SCRIPT
#=====================================================

print "Start."

#=====Ion and bunch parameters=====================

ion_charge = 31
ion_mass = 197
ion_energy = 0.0032 #Kinetic energy per nucleon, GeV/nucl  

mass = ion_mass*0.9314940954
charge = ion_charge
kin_energy = ion_mass*ion_energy 

b_init = Bunch()
b_init.mass(mass)
b_init.charge(charge)
b_init.getSyncParticle().kinEnergy(kin_energy)

#--------------------------------------------------
# This bunch will have only four particles
# 1. with non-zero initial postion in x-direction
# 2. with non-zero initial angle in x-direction
# 3. with non-zero initial postion in y-direction
# 4. with non-zero initial angle in y-direction
#--------------------------------------------------

#---- transverse initial deviation in meters, angle in radians
transv_deviation = 0.0001
transv_angle = 0.0002

#----- coordinates (x,x',y,y',z,dE)
b_init.addParticle(transv_deviation,0.,0.,0.,0.,0.)
b_init.addParticle(0.,transv_angle,0.,0.,0.,0.)
b_init.addParticle(0.,0.,transv_deviation,0.,0.,0.)
b_init.addParticle(0.,0.,0.,transv_angle,0.,0.)

b_init.dumpBunch()

#-----------------------------------------
# Let's make a lattice with 1 Quad
#-----------------------------------------
print "========================================="
quad_Length = 1.0    #m
quad_k1 = 1.0        #normalized field gradient k1=G/b_ro, [1/m/m]
n_Splits = 1;

teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")

for i in range(n_Splits):
	elem = teapot.QuadTEAPOT("quad"+str(i))
	elem.setLength(quad_Length/n_Splits)
	elem.addParam("kq",quad_k1)
        teapot_lattice.addNode(elem)

teapot_lattice.initialize()

#---- dictionary with node positions position_dict[node] = [start_pos,end_pos]
position_dict = teapot_lattice.getNodePositionsDict()

print "========Lattice==========================="
nodes = teapot_lattice.getNodes()
for index,node in enumerate(nodes):
	st  = "i=" + str(index)
	st += " node=" + node.getName()
	st += " type=" + node.getType()
	st += " length=" + str(node.getLength())
	st += " position (start,stop) =" + str(position_dict[node])
	print st

print "========================================="

lattice = teapot_lattice
lattice_length = teapot_lattice.getLength()
print "Line total length[m]= %6.2f "%(lattice_length)
print "========================================="


#-----------------------------------------
# Let's track the bunch through the lattice
#-----------------------------------------

#---- we need this parameter because by default the TEAPOT elements will
#---- assume the charge equal +1 as for proton. There will be no assumptions 
#---- about mass and energy

b = Bunch()
b_init.copyBunchTo(b)

paramsDict = {}
paramsDict["useCharge"]=1

lattice.trackBunch(b,paramsDict)

#--------------------------------------------------------------
# Let's compare bunches after TEAPOT and Matrix tracking
#--------------------------------------------------------------

print "======Final Bunch after TEAPOT Tracking =========================="
b.dumpBunch()

matrx = getTransportMatrix(quad_k1, charge, quad_Length)
bunch_matrix = trackBunchThroughTransportMatrix(b_init,matrx)
print "======Final Bunch after MATRIX Tracking =========================="
bunch_matrix.dumpBunch()

print "=================STOP========================"