##############################################################
# This script tests the transport matrix generator from
# one bunch that has the Initial Coordinates particles Attrubute.
# The coordinates in the attributes are considered as initial ones.
# The matrix is a 7x7 matrix that transforms the initial particles 
# coordinates to the final ones that are in the bunch.
#
# This example also includes the Twiss filtering function.
##############################################################

import math
import sys
import time
import random

from bunch import Bunch

import orbit_utils
from orbit_utils import bunch_utils_functions
from bunch_utils_functions import copyCoordsToInitCoordsAttr
from bunch_utils_functions import swapInitCoordsAttrAndCoords
from bunch_utils_functions import transportMtrxFromInitCoords

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

#---- Copy the coordinates to the initial coordinates array of Particles Attribute
copyCoordsToInitCoordsAttr(b_in)

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
#------ track coordinates of the particle with the transport matrix 
noise_level = 0.000001
nParticles = b_in.getSize()
for i in range(nParticles):
	coord_arr = (b_in.x(i),b_in.xp(i),b_in.y(i),b_in.yp(i),b_in.z(i),b_in.dE(i))
	coord_arr_res = [0.,0.,0.,0.,0.,0.]
	for ix in range(6):
		coord_arr_res[ix] = mtrxA_init.get(ix,6)
		for iy in range(6):
			coord_arr_res[ix] += mtrxA_init.get(ix,iy)*coord_arr[iy]
		coord_arr_res[ix] += noise_level*random.random()
	[x,xp,y,yp,z,dE] = coord_arr_res
	b_in.x(i,x)
	b_in.xp(i,xp)
	b_in.y(i,y)
	b_in.yp(i,yp)
	b_in.z(i,z)
	b_in.dE(i,dE)

#---- let's remove 20% of the particles from the out-bunch (they are lost!)
n_deleted = int(nParticles*0.2)
for ind in range(n_deleted):
	ind = int(nParticles*random.random())
	if(ind > nParticles-1): ind = nParticles-1
	b_in.deleteParticleFast(ind)

b_in.compress()

print "debug before filtering n_part in ==in == bunch=",b_in.getSizeGlobal()

#----- Twiss filtering. 
#----- bunch_utils_functions.bunchTwissFiltering(b_in,b_bad,coeff_x,coeff_y,coeff_z)
#----- b_bad - collection of removed macro-particles
#----- coeff_x the cutt-off coefficients for the value (if coeff_x < 0 no filtering)
#------(gamma*x^2+2*alpha*x*xp+beta*xp^2)/(2*emittance)
#------ for the x-direction
bunch_bad = Bunch()
bunch_utils_functions.bunchTwissFiltering(b_in,bunch_bad,5.0,-1.0,-1.0)

print "debug after  filtering n_part in ==in == bunch=",b_in.getSizeGlobal()

count = 0
while(1 < 2):

	#-----get the matrix
	mtrxA = Matrix(7,7)
	n_part_analysis = transportMtrxFromInitCoords(b_in,mtrxA,1,1,1)

	#==== we can use transport matrix generator function without Twiss weights for macroparticles 
	#n_part_analysis = transportMtrxFromInitCoords(b_in,mtrxA)

	printM(mtrxA, "Transp. M ")
	print "Total N=",n_part_analysis," count=",count
	count += 1

print "Stop."


