import math
import sys
from bunch import Bunch
from foil import Foil
from injection import InjectParts
print "Start."

xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050

foilparams = (xmin, xmax, ymin, ymax)

#------------------------------
#Bunch init
#------------------------------
b = Bunch()
runName = "Test_Injection"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)

lostfoilbunch = Bunch()
lostfoilbunch.addPartAttr("LostParticleAttributes") 

#------------------------------
#Initial Distribution Functions
#------------------------------

xFunc = JohoTransverse()
yFunc = JohoTransverse()
lFunc = JohoLongitudinal()

#------------------------------
# Inject some particles
#------------------------------

nparts = 10

addParticles(nparts, b, lostfoilbunch, foilparams, xFunc, yFunc, lFunc)

b.dumpBunch()
print "Stop."



