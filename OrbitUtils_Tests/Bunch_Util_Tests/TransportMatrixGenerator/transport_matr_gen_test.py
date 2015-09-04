##############################################################
# This script tests the transport matrix generator from
# two bunches that are the states of the bunch at two
# different places in the lattice. The matrix is a 7x7 matrix
# that transforms the particles coordinates.
##############################################################

import math
import sys
import time
import random

from bunch import Bunch
from orbit.bunch_utils import ParticleIdNumber

import orbit_utils
from orbit_utils import bunch_utils_functions

from orbit_utils import Matrix

#------------------------------------------------------
def printM(m, st=""):
	print st,"----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("m(" + str(i) + "," + str(j)+")="+"%12.5g"%m.get(i,j) + " "),
		print ""	
#------------------------------------------------------

print "Start."

#------------------------------
# Initial Bunch
#------------------------------
nParticles = 10
b_in = Bunch()
for i in range(nParticles):
	(x,xp,y,yp,z,dE) = (random.random(),random.random(),random.random(),random.random(),random.random(),random.random())
	b_in.addParticle(x,xp,y,yp,z,dE)
	
#----set up ids
ParticleIdNumber.addParticleIdNumbers(b_in)

b_out = Bunch()
b_in.copyEmptyBunchTo(b_out)

#-----set up initial transport matix, later we have to extract it from bunches
mtrxA_init = Matrix(7,7)
mtrxA_init.zero()
for ix in range(6):
	mtrxA_init.set(ix,6,1.0*ix)
mtrxA_init.set(6,6,1.0)
for ix in range(6):
	mtrxA_init.set(ix,ix,1.*(ix+1))
printM(mtrxA_init, "Init M ")
#-----------------------------------------------------------------------------
#------ put particles from b_in bunch to b_out 
noise_level = 0.000001
for i in range(nParticles):
	coord_arr = (b_in.x(i),b_in.xp(i),b_in.y(i),b_in.yp(i),b_in.z(i),b_in.dE(i))
	coord_arr_res = [0.,0.,0.,0.,0.,0.]
	for ix in range(6):
		coord_arr_res[ix] = mtrxA_init.get(ix,6)
		for iy in range(6):
			coord_arr_res[ix] += mtrxA_init.get(ix,iy)*coord_arr[iy]
		coord_arr_res[ix] += noise_level*random.random()
	[x,xp,y,yp,z,dE] = coord_arr_res
	b_out.addParticle(x,xp,y,yp,z,dE)
	part_id = b_in.partAttrValue("ParticleIdNumber", i, 0)
	b_out.partAttrValue("ParticleIdNumber", i, 0, part_id)

b_out.deleteParticleFast(int(nParticles*0.5))
b_out.deleteParticleFast(int(nParticles*0.3))

#-----get the matrix
mtrxA = Matrix(7,7)
n_part_analysis = bunch_utils_functions.trasportMtrx(b_in,b_out,mtrxA)
print "Total N=",n_part_analysis
printM(mtrxA, "Transp. M ")

print "Stop."


