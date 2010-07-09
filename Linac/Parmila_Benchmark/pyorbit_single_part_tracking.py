#! /usr/bin/env python

"""
This script will track the bunch with one particle through the SNS Linac
and print the coordinates into the file.
"""

import sys
import math
import random

from orbit.sns_linac import SimplifiedLinacParser
from orbit.sns_linac import LinacLatticeFactory, LinacAccLattice

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D


from bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

parser = SimplifiedLinacParser("../SNS_Linac_XML/sns_linac.xml")
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
lattFactory.setMaxDriftLength(0.03)
accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh"])
accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6"])

print "Acc Lattice is ready. "
#set H- mass
#self.bunch.mass(0.9382723 + 2*0.000511)
bunch = Bunch()
bunch.mass(0.939294)
bunch.charge(-1.0)
bunch.getSyncParticle().kinEnergy(0.0025)
bunch.addParticle(0.000,0.0,0.000,0.0,0.001,0.0)

#set up design
accLattice.trackDesignBunch(bunch)

print "Design tracking completed."

#track through the lattice 
paramsDict = {"test_pos":0.,"count":0}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")

print " N  name position[m]  x[mm]  xp[mrad]    y[mm]  yp[mrad]   z[mm]   dE[keV]  eKin[MeV]  "
file_out = open("pyorbit_trajectory_ekin.dat","w")
file_out.write(" N  name position[m]  x[mm]  xp[mrad]    y[mm]  yp[mrad]   z[mm]   dE[keV]  eKin[MeV]  \n")

def action_entrance(paramsDict):
	bunch = paramsDict["bunch"]
	node = paramsDict["node"]	
	length = node.getLength()
	#print "debug ============= entr xp=",bunch.xp(0), "   name=",node.getName()," L=",length
	if(isinstance(paramsDict["parentNode"],AccLattice)):
		pos = paramsDict["test_pos"]
		(x,xp,y,yp,z,dE) = (bunch.x(0)*1000.,bunch.xp(0)*1000.,bunch.y(0)*1000.,bunch.yp(0)*1000.,bunch.z(0)*1000.,bunch.dE(0)*1000.*1000.)
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f  %9.6f  %9.6f  %9.6f  %9.3f  %9.6f  %9.5f   %12.6f  "%(paramsDict["count"],node.getName(),pos,x,xp,y,yp,z,dE,eKin)
		file_out.write(s +"\n")
		print s		
	

def action_exit(paramsDict):
	node = paramsDict["node"]
	length = node.getLength()
	pos = paramsDict["test_pos"] + length
	paramsDict["test_pos"] = pos	
	bunch = paramsDict["bunch"]
	#print "debug ============= exit xp=",bunch.xp(0), "   name=",node.getName()," L=",length
	if(isinstance(paramsDict["parentNode"],AccLattice)):	
		paramsDict["count"]	+= 1
		(x,xp,y,yp,z,dE) = (bunch.x(0)*1000.,bunch.xp(0)*1000.,bunch.y(0)*1000.,bunch.yp(0)*1000.,bunch.z(0)*1000.,bunch.dE(0)*1000.*1000.)
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f  %9.6f  %9.6f  %9.6f  %9.3f  %9.6f  %9.5f   %12.6f  "%(paramsDict["count"],node.getName(),pos,x,xp,y,yp,z,dE,eKin)
		file_out.write(s +"\n")
		print s	
	
actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)
accLattice.trackBunch(bunch, paramsDict = paramsDict, actionContainer = actionContainer)

file_out.close()









