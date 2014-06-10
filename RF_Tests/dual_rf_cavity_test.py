import sys

sys.path.append("/u/sappel/py-orbit/py")

from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from rfcavities import Harmonic_Cav
from rfcavities import Dual_Harmonic_Cav


#///////////////////////////////////////////////////////////
from orbit.rf_cavities import RFNode, RFLatticeModifications
#///////////////////////////////////////////////////////////

import math

print "Start."

b = Bunch()
b.addParticle(1.0e-3,0.0,0.0,0.0,0.0,0.0)
b.addParticle(0.0,1.0e-3,0.0,0.0,0.0,0.0)
b.addParticle(0.0,0.0,1.0e-3,0.0,0.0,0.0)
b.addParticle(0.0,0.0,0.0,1.0e-3,0.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,1.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,0.0,1.0e-3)
b.compress()

syncPart = b.getSyncParticle()
energy = 1.0                          #energy in GeV
#p = syncPart.energyToMomentum(energy)
#syncPart.pz(p)
syncPart.kinEnergy(energy)

lattice = AccLattice("test_lattice")

elem0 = teapot.DriftTEAPOT("drift0")

lattice.addNode(elem0)

#-----------------------------
# Set TEAPOT nodes parameters
#-----------------------------

elem0.setLength(4.0)

lattice.initialize()

#///////////////////////////////////////////////////////////
ZtoPhi = 2.0 * math.pi / lattice.getLength()
RFHNum = 1.0
RatioRFHNum = 2.0
RFVoltage = 0.1
RatioVoltage = 0.5
RFPhase = 1.0
RFPhase2 = 1.0
length = 0.0
name = "dual_harmonic_rfnode"
position = 0.0

rf_node = RFNode.Dual_Harmonic_RFNode(ZtoPhi, RFHNum, RatioRFHNum, RFVoltage, RatioVoltage, RFPhase, RFPhase2, length, name)
RFLatticeModifications.addRFNode(lattice, position, rf_node)



print "Lattice length = ", lattice.getLength()
print "ZtoPhi = ", ZtoPhi
print "RFHNum = ", RFHNum
print "RFVoltage = ", RFVoltage
print "RatioVoltage = ", RatioVoltage
print "RFPhase = ", RFPhase
#///////////////////////////////////////////////////////////

print "==============BEFORE============================"
b.dumpBunch()
print "=========================================="

#=====track action ============
def bodyAction(paramsDict):
    node = paramsDict["node"]
    node.track(paramsDict)


accContainer = AccActionsContainer()
accContainer.addAction(bodyAction,AccActionsContainer.BODY)

paramsDict = {}
paramsDict["bunch"] = b

for i in range(10):
    lattice.trackActions(accContainer,paramsDict)

print "=============AFTER============================="
b.dumpBunch()
print "=========================================="

print "lattice length=",lattice.getLength()
print "beta=",b.getSyncParticle().beta()
print "TEAPOT time[sec]=",b.getSyncParticle().time()
print "SIMPLE time[sec]=",lattice.getLength()/(b.getSyncParticle().beta()*2.99792458e+8)
print "Stop."


