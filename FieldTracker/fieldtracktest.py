
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
print "FILE NEEDS TO NOT CONTAIN A HEADER OR BLANK LINES"
data = myparser.parse("BIGDATA~", -20.0,20.0,-12.0, 12.0 ,0.0, 250.0 , 0.5 , 0.5 , 0.5)

# for i in range(10):
#   print "At Coordinates " , i , ", " , i , ", " , i 
#   print BXGrid.getValueOnGrid(i,i,i), BYGrid.getValueOnGrid(i,i,i), BZGrid.getValueOnGrid(i,i,i)

## Variables to test BGrid3D
BXGrid = data[0]
BYGrid = data[1]
BZGrid = data[2]
fieldgrid3DMag = data[3]
XGrid = data[4]
YGrid = data[5]
ZGrid = data[6]

print "X, Y, Z Grid Size: ", len(XGrid), len(YGrid), len(ZGrid)

nXGrid = len(XGrid)
nYGrid = len(YGrid)
nZGrid = len(ZGrid)

xField3D = 0.2
yField3D = 0.2
zField3D = 0.2

BxField3D = 0.0
ByField3D = 0.0
BzField3D = 0.0

zsymmetry = 1



mytracker = FieldTracker(0)
mytracker.trackBunch(b)

FieldArr = mytracker.BGrid3D(xField3D,yField3D,zField3D,
XGrid,YGrid,ZGrid,
nXGrid, nYGrid, nZGrid,
BxField3D, ByField3D, BzField3D,
BXGrid,  BYGrid, BZGrid,
zsymmetry)





b.dumpBunch("final.dat")

print "Stop."



