#! /usr/bin/env python

"""
This is a test script to check the functionality of the 
linac acc lattice.
"""

import sys

from orbit.sns_linac import SimplifiedLinacParser
from orbit.sns_linac import LinacLatticeFactory, LinacAccLattice

from bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

parser = SimplifiedLinacParser("sns_linac.xml")
linacTree = parser.getLinacStructureTree()
print "======================================="
print "Total length=",linacTree.getLength()
print "======================================="
sequences = linacTree.getSeqs()
totalLength = 0.
for seq in sequences:
	totalLength +=  seq.getLength()	
	print "seq=",seq.getName()," L=",seq.getLength(),"  total length=",totalLength

lattFactory = 	LinacLatticeFactory(linacTree)
accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh"])

b = Bunch()
syncPart = b.getSyncParticle()
#set H- mass
b.mass(0.9382723 + 2*0.000511)
b.charge(-1.0)
syncPart.kinEnergy(0.0025)

#set up design
paramsDict = {"test_pos":0.,"count":0}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")

def test_action(paramsDict):
	node = paramsDict["node"]
	length = node.getLength()
	pos = paramsDict["test_pos"] + length
	paramsDict["test_pos"] = pos	
	bunch = paramsDict["bunch"]
	eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3	
	if(node.getName().find(":Rg") >= 0):
		paramsDict["count"]	+= 1
		print "test debug i=",paramsDict["count"]," node=",node.getName()," pos=",(pos - length/2)," Ekin[MeV]=",eKin 

actionContainer.addAction(test_action, AccActionsContainer.EXIT)


accLattice.trackDesignBunch(b, paramsDict = paramsDict, actionContainer = actionContainer)

sys.exit(1)

