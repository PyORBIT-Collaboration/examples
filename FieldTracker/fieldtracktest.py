
import math
import sys
from bunch import Bunch
from fieldtracker import FieldTracker
from orbit.parsers.field_parser import Field_Parser3D


print "Start."


#------------------------------
#Main Bunch init
#------------------------------

b = Bunch()
print "Read Bunch."
runName = "Test_FieldTracker"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)

#=====track bunch through Foil============

myparser = Field_Parser3D()
myparser.parse("testfile", 1, 1, 1)

mytracker = FieldTracker(0)
mytracker.trackBunch(b)
b.dumpBunch("final.dat")

print "Stop."



