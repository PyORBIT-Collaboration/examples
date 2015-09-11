#------------------------------------------------
#pyORBIT error and correction example
#------------------------------------------------

import sys
import math
from pylab import *	


import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.errors import AddErrorNode
from orbit.errors import AddErrorSet

from orbit.teapot import TEAPOT_MATRIX_Lattice


from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

# ATTENTION !!! The python packet numpy and scipy are required
from orbit.orbit_correction import orbit, correction   



print "Start."
#---------------------------------------------Bunch init---------------------------------------------
b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0)

energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)
#---------------------------------------------Bunch init---------------------------------------------

print "Generate Lattice."
#---------------------------------------------Make a Teapot Lattice----------------------------------
lattice = TEAPOT_Lattice("lattice")
lattice.readMAD("sis18.lat","SIS18")
#---------------------------------------------Make a Teapot Lattice----------------------------------

print "INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES"
#---------------------------------------------ORBIT ERRORS-------------------------------------------
# WE INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES; dx, dy = HOR AND VER DISPLACEMENT OF QUADRUPOLES
setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "TransDisp"
paramsDict["sample"]      = "Uniform"
paramsDict["maximum"]        = 0.5
paramsDict["minimum"]       = 0.0
paramsDict["dx"]       = 0.005
paramsDict["dy"]       = 0.007

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=50)
# ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
#---------------------------------------------ORBIT ERRORS-------------------------------------------


print "CALCULATE DISTORTED ORBIT AND PLOT"
#---------------------------------------------CALCULATE DISTORTED ORBIT------------------------------
OrbitX, OrbitY = orbit(lattice,b).get_orbit()
#---------------------------------------------CALCULATE DISTORTED ORBIT------------------------------

#---------------------------------------------PLOT DISTORTED ORBIT-----------------------------------
x = []
y = []
s = []
for i in xrange(len(OrbitX)):
	s.append(OrbitX[i][0])
	x.append(OrbitX[i][1])
	y.append(OrbitY[i][1])

plot(s,x,'r-', label="x")
plot(s,y,'b-', label="y")
#---------------------------------------------PLOT DISTORTED ORBIT-----------------------------------

print "CORRECTED ORBIT"
#---------------------------------------------CORRECTED ORBIT----------------------------------------
corr = correction(lattice,b)
corr.orbit_corr()
#---------------------------------------------CORRECTED ORBIT----------------------------------------

print "CALCULATE CORRECTED ORBIT AND PLOT"
#---------------------------------------------CALCULATE CORRECTED ORBIT------------------------------
OrbitX_corr, OrbitY_corr = orbit(lattice,b).get_orbit()
#---------------------------------------------CALCULATE CORRECTED ORBIT------------------------------

#---------------------------------------------PLOT CORRECTE ORBIT------------------------------------
x_corr = []
y_corr = []
for i in xrange(len(OrbitX_corr)):
	x_corr.append(OrbitX_corr[i][1])
	y_corr.append(OrbitY_corr[i][1])
	

plot(s,x_corr,'r--', label="x corr")
plot(s,y_corr,'b--', label="y corr")
#---------------------------------------------PLOT CORRECTE ORBIT------------------------------------

xlabel("s [m]")
ylabel("[m]")
legend(loc="lower right")
savefig("plot.pdf")
show()
print("STOP.")



