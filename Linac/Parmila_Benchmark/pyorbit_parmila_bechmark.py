#! /usr/bin/env python

"""
This script will track the bunch through the SNS MEBT in pyORBIT and 
will generate the intermediate file for PARMILA
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

from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

random.seed(100)

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
lattFactory.setMaxDriftLength(0.5)
accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh"])
#accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1"])

print "Acc Lattice is ready. "
#-----TWISS Parameters at the entrance of the MEBT ---------------
# transverse emittances are normalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.0025 # in [GeV]
mass = 0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta
frequency = 402.5e+6
v_light = 2.99792458e+8  # in [m/sec]


emittX = 0.21
emittY = 0.21
emittZ = 7.6e-7

emittX = (emittX/(beta*gamma))*1.0e-6
emittY = (emittY/(beta*gamma))*1.0e-6
emittZ = emittZ*1.0e-9*v_light*beta

alphaX = -1.962
betaX = 0.0183*1.0e+3*1.0e-2

alphaY = 1.768
betaY = 0.0161*1.0e+3*1.0e-2

# Parmila betaZ in [deg/MeV] and we need [m/GeV]
alphaZ = 0.0196
betaZ = (772.8/360.)*(v_light*beta/frequency)*1.0e+3

print " aplha beta emitt X=",alphaX,betaX,emittX
print " aplha beta emitt Y=",alphaY,betaY,emittY
print " aplha beta emitt Z=",alphaZ,betaZ,emittZ

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

xal_emittZ = emittZ/(gamma**3*beta**2*mass)
xal_betaZ = betaZ*(gamma**3*beta**2*mass)
print "XAL Twiss Longitudinal parameters alpha=",alphaZ," beta=", xal_betaZ," emittZ =",xal_emittZ
print "==============================================="

print "Start Bunch Generation."
bunch_gen = SNS_Linac_BunchGenerator(twissX,twissY,twissZ)

bunch_in = bunch_gen.getBunch(nParticles = 5000, distributorClass = WaterBagDist3D)

bunch_gen.dumpParmilaFile(bunch_in, phase_init = -45.0, fileName = 	"parmila_bunch.txt")
print "Bunch Generation completed."
#set up design
accLattice.trackDesignBunch(bunch_in)

print "Design tracking completed."

#track through the lattice 
paramsDict = {"test_pos":0.,"count":0}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")

twiss_analysis = TwissAnalysis(3)



print "   N           node           position         sizeX       sizeY    sizeZ    sizeXP   sizeYP   size_dE   eKin "
file_out = open("pyorbit_sizes_ekin.dat","w")
file_out.write(" N           node   position  sizeX  sizeY  sizeZ  sizeXP  sizeYP sizedE  eKin Nparts \n")

def action_entrance(paramsDict):
	if(isinstance(paramsDict["parentNode"],AccLattice)):
		node = paramsDict["node"]
		pos = paramsDict["test_pos"]
		bunch = paramsDict["bunch"]
		twiss_analysis.init()
		for i in range(bunch.getSize()):
			(x,xp,y,yp,z,dE) = (bunch.x(i),bunch.xp(i),bunch.y(i),bunch.yp(i),bunch.z(i),bunch.dE(i))
			twiss_analysis.account((x,xp,y,yp,z,dE))
		x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
		y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
		z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
		xp_rms = math.sqrt(twiss_analysis.getTwiss(0)[2]*twiss_analysis.getTwiss(0)[3])*1000.
		yp_rms = math.sqrt(twiss_analysis.getTwiss(1)[2]*twiss_analysis.getTwiss(1)[3])*1000.
		dE_rms = math.sqrt(twiss_analysis.getTwiss(2)[2]*twiss_analysis.getTwiss(2)[3])*1000. 
		#emittX = twiss_analysis.getTwiss(0)[3]*1000.0*1000.0	*bunch.getSyncParticle().gamma()*bunch.getSyncParticle().beta()
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f  %5.3f  %5.3f   %5.3f  %5.3f  %5.3f  %7.5f  %10.6f   %8d "%(paramsDict["count"],node.getName(),pos,x_rms,y_rms,z_rms,xp_rms,yp_rms,dE_rms,eKin,bunch.getSize())
		file_out.write(s +"\n")
		print s	
	

def action_exit(paramsDict):
	node = paramsDict["node"]
	length = node.getLength()
	pos = paramsDict["test_pos"] + length
	paramsDict["test_pos"] = pos	
	if(isinstance(paramsDict["parentNode"],AccLattice)):	
		bunch = paramsDict["bunch"]
		paramsDict["count"]	+= 1
		twiss_analysis.init()
		for i in range(bunch.getSize()):
			(x,xp,y,yp,z,dE) = (bunch.x(i),bunch.xp(i),bunch.y(i),bunch.yp(i),bunch.z(i),bunch.dE(i))
			twiss_analysis.account((x,xp,y,yp,z,dE))
		x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
		y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
		z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
		xp_rms = math.sqrt(twiss_analysis.getTwiss(0)[2]*twiss_analysis.getTwiss(0)[3])*1000.
		yp_rms = math.sqrt(twiss_analysis.getTwiss(1)[2]*twiss_analysis.getTwiss(1)[3])*1000.
		dE_rms = math.sqrt(twiss_analysis.getTwiss(2)[2]*twiss_analysis.getTwiss(2)[3])*1000. 
		#emittX = twiss_analysis.getTwiss(0)[3]*1000.0*1000.0	*bunch.getSyncParticle().gamma()*bunch.getSyncParticle().beta()
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f  %5.3f  %5.3f   %5.3f  %5.3f  %5.3f  %7.5f  %10.6f   %8d "%(paramsDict["count"],node.getName(),pos,x_rms,y_rms,z_rms,xp_rms,yp_rms,dE_rms,eKin,bunch.getSize())
		file_out.write(s +"\n")
		print s	
	
actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)
accLattice.trackBunch(bunch_in, paramsDict = paramsDict, actionContainer = actionContainer)

file_out.close()









