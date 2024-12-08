##############################################################
# This script tests set the macrosizes to bunch or particles
##############################################################

import math
import sys
import time

from bunch import Bunch

print "Start."

#------------------------------
#Main Bunch init
#------------------------------

nParticles = 5
bunch = Bunch()
(x,xp,y,yp,z,dE) = (0.,0.,0.,0.,0.,0.)
for ind in range(nParticles):
	bunch.addParticle(x,xp,y,yp,z,dE)
	
#---- The macro_size for Bunch instance
#---- after that each added or exiting macro-particle in the bunch
#---- will have a charge of 1000 electrons. So, If you add a macro-particle
#---- to bunch, you will increase the total charge of the bunch.
#---- 1.602176634e-19 - charge of 1 electron (abs. value)
macro_size = 1000*1.602176634e-19
bunch.macroSize(macro_size)
print "================init bunch==============="
bunch.dumpBunch()
print "========================================="

#---- add particles macrosize attribut
#---- after that macro-particles can have different macro-sizes.
#---- by default the value is zero.
#---- After you setup macrosize attribute to particles any space-charge will be 
#---- defined by these macrosizes, the clobal macro-size for the bunch object 
#---- will be ingnored.
bunch.addPartAttr("macrosize")
print "================bunch with macrosize=0 particles attribute ==============="
bunch.dumpBunch()
print "=========================================================================="
atribute_array_index = 0
for ind in range(nParticles):
	#---- set macrosize for particle with index=ind
	macro = macro_size*ind
	bunch.partAttrValue("macrosize", ind, atribute_array_index, macro)
  #---- get macrosize for particle with index=ind
	new_macro_size = bunch.partAttrValue("macrosize", ind, atribute_array_index)
	print "ind =",ind," size=",new_macro_size
	
print "================bunch with macrosize particles attribute ==============="
bunch.dumpBunch()
print "=========================================================================="

print "Stop."
