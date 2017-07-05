##############################################################
# This script tests the transport matrix generator from
# two bunches that are the states of the bunch at two
# different places in the lattice. The matrix is a 7x7 matrix
# that transforms the particles coordinates.
# Both bunches have Id number particles Attributes to define
# the relations between particles from two bunches.
#
# This example also includes the Twiss filtering function.
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
nParticles = 10000
b_in = Bunch()
b_in.addPartAttr("macrosize")
b_in_has_ms = b_in.hasPartAttr("macrosize")
for i in range(nParticles):
	(x,xp,y,yp,z,dE) = (random.random(),random.random(),random.random(),random.random(),random.random(),random.random())
	b_in.addParticle(x,xp,y,yp,z,dE)
	if(b_in_has_ms):
		b_in.partAttrValue("macrosize", i, 0, random.random())
b_in.addParticle(0.,0.,0.,0.,0.,0.)	
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

mtrxA_init.set(0,1,1.5*.001)
mtrxA_init.set(2,3,1.6*.001)
mtrxA_init.set(4,5,1.7*.001)

mtrxA_init.set(1,0,1.4*.001)
mtrxA_init.set(3,2,1.3*.001)
mtrxA_init.set(5,4,1.2*.001)

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
	if(b_in_has_ms):
		m_size = b_in.partAttrValue("macrosize", i, 0)
		b_out.partAttrValue("macrosize", i, 0, m_size)

#---- let's remove 20% of the particles from the out-bunch (they are lost!)
n_deleted = int(nParticles*0.2)
for ind in range(n_deleted):
	ind = int(nParticles*random.random())
	if(ind > nParticles-1): ind = nParticles-1
	b_out.deleteParticleFast(ind)

b_out.compress()
b_in.compress()

print "debug before filtering n_part in ==in == bunch=",b_in.getSizeGlobal()
print "debug before filtering n_part in ==out== bunch=",b_out.getSizeGlobal()

#----- Twiss filtering. 
#----- bunch_utils_functions.bunchTwissFiltering(b_in,b_bad,coeff_x,coeff_y,coeff_z)
#----- b_bad - collection of removed macro-particles
#----- coeff_x the cutt-off coefficients for the value (if coeff_x < 0 no filtering)
#------(gamma*x^2+2*alpha*x*xp+beta*xp^2)/(2*emittance)
#------ for the x-direction
bunch_bad = Bunch()
bunch_utils_functions.bunchTwissFiltering(b_in,bunch_bad,1.0,-1.0,-1.0)
bunch_utils_functions.bunchTwissFiltering(b_out,bunch_bad,5.0,-1.0,-1.0)

print "debug after  filtering n_part in ==in == bunch=",b_in.getSizeGlobal()
print "debug after  filtering n_part in ==out== bunch=",b_out.getSizeGlobal()


count = 0
while(1 < 2):

	#-----get the matrix
	mtrxA = Matrix(7,7)
	n_part_analysis = bunch_utils_functions.transportMtrx(b_in,b_out,mtrxA,1,1,1)

	#==== we can use transport matrix generator function without Twiss weights for macroparticles 
	#n_part_analysis = bunch_utils_functions.transportMtrx(b_in,b_out,mtrxA)

	printM(mtrxA, "Transp. M ")
	print "Total N=",n_part_analysis," count=",count
	count += 1
	sys.exit()
print "Stop."


