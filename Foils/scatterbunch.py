
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
# Below is 1000 times the width of normal foil but will do only one turn.
thick = 400

foil = Foil(xmin, xmax, ymin, ymax, thick)

#------------------------------
#Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "Benchmark_Collimator"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)

#=====track bunch through Foil============

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes") 

foil.traverseFoilFullScatter(b, lostbunch)
b.dumpBunch("scatteredbunch.dat")
lostbunch.dumpBunch("lostbunch.dat")

print "Stop."



