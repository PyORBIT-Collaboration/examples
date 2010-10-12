#-----------------------------------------------------
#Creates Space Charge Calculator based on Rick Baartman 
#suggestion and tracks the test bunch through the calculator. 
#-----------------------------------------------------


import sys
import math

import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5Drb

print "Start."

sizeX = 128
sizeY = 128
sizeZ = 2
calc2p5d = SpaceChargeCalc2p5Drb(sizeX,sizeY,sizeZ)

charge = 1.0

b = Bunch()
b.addPartAttr("macrosize")

macroSize = 1.0e-13

energy = 1.0
b.getSyncParticle().kinEnergy(energy)

b.addParticle( 0.,0.,0.,0.,0.,0.)
b.addParticle( 0.001,0.,0.,0.,0.,0.)
b.addParticle(-0.001,0.,0.,0.,0.,0.)
b.addParticle( 0.,0., 0.001,0.,0.,0.)
b.addParticle( 0.,0.,-0.001,0.,0.,0.)
b.addParticle( 0.,0.,0.,0., 0.001,0.)
b.addParticle( 0.,0.,0.,0.,-0.001,0.)

b.partAttrValue("macrosize",0,0,macroSize)
b.partAttrValue("macrosize",1,0,0.)
b.partAttrValue("macrosize",2,0,0.)
b.partAttrValue("macrosize",3,0,0.)
b.partAttrValue("macrosize",4,0,0.)
b.partAttrValue("macrosize",5,0,0.)
b.partAttrValue("macrosize",6,0,0.)

nParts = b.getSize()
macroSize = macroSize/nParts
b.macroSize(macroSize)

print "bunchSize = ",b.getSize()
print "macroSize=",b.macroSize()
print "mass=",b.mass()

pipe_radius = 0.005
slice_length = 0.1

b.dumpBunch("pyorbit_bunch_test.in")

calc2p5d.trackBunch(b,slice_length,pipe_radius)

b.dumpBunch("pyorbit_bunch_test.out")

print "Stop."

