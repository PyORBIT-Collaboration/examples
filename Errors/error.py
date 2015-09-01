#------------------------------------------------
#pyORBIT error module benchmark
#------------------------------------------------

import sys
import math

import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.errors import AddErrorNode
from orbit.errors import AddErrorSet

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

print "Start."
#------------------------------
#Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
b.readBunch("pyorbit_bunch_input.dat")
b.mass(0.93827231)
b.macroSize(1.0)

energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

#------------------------------
#Make a Teapot Lattice
#------------------------------

print "Generate Lattice."
lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD("./LATTICES/Test.LAT","TEST")

setDict = {}
paramsDict = {}

"""
#####################################################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "TransDisp"
paramsDict["sample"]   = "Fixed"
paramsDict["dx"]       = 0.1
paramsDict["dy"]       = 0.1

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "LongDisp"
paramsDict["sample"]   = "Fixed"
paramsDict["ds"]       = 0.3

#####################################################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "XYRot"
paramsDict["sample"]   = "Fixed"
paramsDict["angle"]    = 0.2

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "XSRot"
paramsDict["sample"]   = "Fixed"
paramsDict["angle"]    = 0.1

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "YSRot"
paramsDict["sample"]   = "Fixed"
paramsDict["angle"]    = 0.1

###################################

positioni = 2.3
positionf = 2.4
paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "KickField"
paramsDict["sample"]   = "Fixed"
paramsDict["fracerr"]  = 0.2

###################################

positioni = 0.6
positionf = 0.7
paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "SolenoidField"
paramsDict["sample"]   = "Fixed"
paramsDict["fracerr"]  = 10.0

###################################

positioni = 3.6
positionf = 3.7
paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "MultipoleField"
paramsDict["sample"]   = "Fixed"
paramsDict["fracerr"]  = 50.0

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"
paramsDict["sample"]   = "Fixed"
paramsDict["fracerr"]  = 0.5

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "BendField"
paramsDict["sample"]   = "Fixed"
paramsDict["fracerr"]  = 0.2

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]  = "BendDisplacementError"
paramsDict["subtype"]  = "XDisp"
paramsDict["sample"]   = "Fixed"
paramsDict["disp"]     = 0.005

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]  = "BendDisplacementError"
paramsDict["subtype"]  = "YDisp"
paramsDict["sample"]   = "Fixed"
paramsDict["disp"]     = 0.005

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]  = "BendDisplacementError"
paramsDict["subtype"]  = "LongDisp"
paramsDict["sample"]   = "Fixed"
paramsDict["disp"]     = 0.005

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "sbend"
paramsDict["subtype"]     = "XY"
paramsDict["sample"]      = "Fixed"
paramsDict["angle"]       = 0.005

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "sbend"
paramsDict["subtype"]     = "XS"
paramsDict["sample"]      = "Fixed"
paramsDict["angle"]       = 0.05

###################################

positioni = 1.1
positionf = 1.9
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "sbend"
paramsDict["subtype"]     = "YS"
paramsDict["sample"]      = "Fixed"
paramsDict["angle"]       = 0.05

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "straight"
paramsDict["subtype"    ] = "XY"
paramsDict["sample"]      = "Fixed"
paramsDict["angle"]       = 0.2

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "straight"
paramsDict["subtype"]     = "XS"
paramsDict["sample"]      = "Fixed"
paramsDict["angle"]       = 0.1

###################################

positioni = 2.8
positionf = 3.2
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "straight"
paramsDict["subtype"]     = "YS"
paramsDict["sample"]      = "Fixed"
paramsDict["angle"]       = 0.1

###################################

ENode = AddErrorNode(lattice, positioni, positionf, paramsDict)

#####################################################################

setDict["elementtype"] = "drift"
setDict["elementtype"] = "kick"
setDict["elementtype"] = "soln"
setDict["elementtype"] = "mult"
setDict["elementtype"] = "quad"
setDict["elementtype"] = "sbend"

setDict["ringline"] = "ring"
setDict["ringline"] = "line"

###################################

paramsDict["sample"]      = "Fixed"
paramsDict["sample"]      = "Uniform"
paramsDict["minimum"]     = -1.0
paramsDict["maximum"]     =  1.0
paramsDict["sample"]      = "Gaussian"
paramsDict["mean"]        = 0.0
paramsDict["sigma"]       = 1.0

###################################

ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict)

#####################################################################

"""

positioni = 0.0
positionf = 6.0
paramsDict["errtype"]     = "RotationError"
paramsDict["elementtype"] = "straight"
paramsDict["subtype"]     = "YS"
paramsDict["sample"]      = "Gaussian"
paramsDict["mean"]        = 0.0
paramsDict["sigma"]       = 1.0
paramsDict["angle"]       = 0.1

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict)

z = 0.0
print "Start lattice, z = ", z

# set the number of sections in quads to the same as for ORBIT_MPI
for acc_elem in lattice.getNodes():
	z += acc_elem.getLength()
	print "Node = ", acc_elem.getName()," type = ", acc_elem.getType(),\
	    " L = ", acc_elem.getLength(), " N child nodes = ",\
	    acc_elem.getNumberOfChildren(), " z = ", z
	if(acc_elem.getType() == "quad teapot"):
		acc_elem.setnParts(5)

print "lattice length=",lattice.getLength()

# dump initial bunch for ORBIT_MPI input
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch_input.dat")

#=====track bunch ============
ACC_TURNS = 1
print("Tracking.")
for i in range(ACC_TURNS):
	lattice.trackBunch(b)
	print "Turn ",i

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch_output.dat")

print("STOP.")
