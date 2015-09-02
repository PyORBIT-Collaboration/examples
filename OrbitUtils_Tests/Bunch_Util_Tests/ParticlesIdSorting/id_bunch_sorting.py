##############################################################
# This script test the bunch particles sorting according to 
# the particle id numbers.
##############################################################

import math
import sys
import time

from bunch import Bunch
from orbit.bunch_utils import ParticleIdNumber

import orbit_utils
from orbit_utils import bunch_utils_functions

print "Start."

#------------------------------
#Main Bunch init
#------------------------------
nParticles = 100000
b = Bunch()
(x,xp,y,yp,z,dE) = (0.,0.,0.,0.,0.,0.)
for i in range(nParticles):
	b.addParticle(x,xp,y,yp,z,dE)
	
#----set up ids
ParticleIdNumber.addParticleIdNumbers(b)

#----- sorting speed measurements
while(1 < 2):

	n_parts = b.getSize()
	for i in range(n_parts):
		b.partAttrValue("ParticleIdNumber", i, 0, n_parts - i)
	
	time_start = time.clock()
	#=====sort particles in the bunch============
	bunch_utils_functions.bunchSortId(b)
	
	tm = (time.clock() - time_start)
	print "debug sorting speed =",tm
	
#---- dump bunch into the file
b.dumpBunch("sorted_bunch.dat")
print "Stop."


