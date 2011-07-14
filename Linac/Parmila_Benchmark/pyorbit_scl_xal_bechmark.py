#! /usr/bin/env python

"""
This script will track the bunch through the SNS SCL in pyORBIT
"""

import sys
import math
import random
import time

from orbit.sns_linac import SimplifiedLinacParser
from orbit.sns_linac import LinacLatticeFactory, LinacAccLattice

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D


from bunch import Bunch
from bunch import BunchTwissAnalysis

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
lattFactory.setMaxDriftLength(0.01)
#accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh"])
accLattice = lattFactory.getLinacAccLattice(["SCLMed",])

#-----------------------------------------------------
# Set up Space Charge Acc Nodes
#-----------------------------------------------------
from orbit.space_charge.sc2p5d import setSC2p5DrbAccNodes
from orbit.space_charge.sc3d import setSC3DAccNodes, setUniformEllipsesSCAccNodes
from spacecharge import SpaceChargeCalc2p5Drb, SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D
sc_path_length_min = 0.01

"""
#the 2p5rb Space Charge nodes
sizeX = 32
sizeY = 32
sizeZ = 32
long_avg_n = 5
calc2p5d = SpaceChargeCalc2p5Drb(sizeX,sizeY,sizeZ)
calc2p5d.setLongAveragingPointsN(long_avg_n)

pipe_radius = 0.015
space_charge_nodes = setSC2p5DrbAccNodes(accLattice,sc_path_length_min,calc2p5d,pipe_radius)
"""

"""
# set of uniformly charged ellipses Space Charge
nEllipses = 1
calcUnifEllips = SpaceChargeCalcUnifEllipse(nEllipses)
space_charge_nodes = setUniformEllipsesSCAccNodes(accLattice,sc_path_length_min,calcUnifEllips)
"""


# set FFT 3D Space Charge
sizeX = 64
sizeY = 64
sizeZ = 64
calc3d = SpaceChargeCalc3D(sizeX,sizeY,sizeZ)
space_charge_nodes =  setSC3DAccNodes(accLattice,sc_path_length_min,calc3d)


"""
max_sc_length = 0.
min_sc_length = accLattice.getLength()
for sc_node in space_charge_nodes:
	scL = sc_node.getLengthOfSC()
	if(scL > max_sc_length): max_sc_length = scL
	if(scL < min_sc_length): min_sc_length = scL
print "maximal SC length =",max_sc_length,"  min=",min_sc_length
"""


print "Acc Lattice is ready. "
#-----TWISS Parameters at the entrance of the SCL ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.185627 # in [GeV]
mass = 0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta
frequency = 402.5e+6
v_light = 2.99792458e+8  # in [m/sec]

"""
#Data from XAL
X alpha, beta emitt=  -1.573126   8.131378   3.34159E-7
Y alpha, beta emitt=  0.684698   8.634511   3.82927E-7
Z alpha, beta emitt=  0.213766   6.071008   2.99985E-7
z_size= 1.349522632222224  delta_E= 111054.43305322397
#the sign of alphaZ should be changed until we fix the pyORBIT Z meaning
"""


emittX = 0.334159
emittY = 0.382927
emittZ = 0.299985

emittX = (emittX)*1.0e-6
emittY = (emittY)*1.0e-6
emittZ = (emittZ)*1.0e-6*gamma**3*beta**2*mass

alphaX = -1.573126
betaX = 8.131378

alphaY = 0.684698
betaY = 8.634511

# betaZ=z_size/emittZ
alphaZ = -0.213766
betaZ = (1.349522e-3)**2/emittZ

print " aplha beta emitt X=",alphaX,betaX,emittX
print " aplha beta emitt Y=",alphaY,betaY,emittY
print " aplha beta emitt Z=",alphaZ,betaZ,emittZ

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

xal_betaZ = betaZ*(gamma**3*beta**2*mass)
xal_emittZ = emittZ/(gamma**3*beta**2*mass)
print "XAL Twiss Longitudinal parameters alpha=",alphaZ," beta=", xal_betaZ," emittZ =",xal_emittZ
print "==============================================="

print "Start Bunch Generation."
bunch_gen = SNS_Linac_BunchGenerator(twissX,twissY,twissZ)

#set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

#set the beam peak current in mA
bunch_gen.setBeamCurrent(38.0)

#bunch_in = bunch_gen.getBunch(nParticles = 20000, distributorClass = WaterBagDist3D)
bunch_in = bunch_gen.getBunch(nParticles = 20000, distributorClass = GaussDist3D)
#bunch_in = bunch_gen.getBunch(nParticles = 20000, distributorClass = KVDist3D)

print "Bunch Generation completed."
#set up design
accLattice.trackDesignBunch(bunch_in)

print "Design tracking completed."

#track through the lattice 
paramsDict = {"test_pos":0.,"count":0}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")

twiss_analysis = BunchTwissAnalysis()

print "   N           node           position         sizeX       sizeY    sizeZ  sizeZdeg  sizeXP   sizeYP   size_dE   eKin Nparts"
file_out = open("pyorbit_scl_sizes_ekin.dat","w")
file_out.write(" N           node   position  sizeX  sizeY  sizeZ  sizeZdeg  sizeXP  sizeYP sizedE  eKin Nparts \n")

def action_entrance(paramsDict):
	if(isinstance(paramsDict["parentNode"],AccLattice)):
		node = paramsDict["node"]
		pos = paramsDict["test_pos"]
		bunch = paramsDict["bunch"]
		twiss_analysis.analyzeBunch(bunch)
		x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
		y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
		z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
		z_rms_deg = bunch_gen.getZtoPhaseCoeff(bunch)*z_rms/1000.0
		xp_rms = math.sqrt(twiss_analysis.getTwiss(0)[2]*twiss_analysis.getTwiss(0)[3])*1000.
		yp_rms = math.sqrt(twiss_analysis.getTwiss(1)[2]*twiss_analysis.getTwiss(1)[3])*1000.
		dE_rms = math.sqrt(twiss_analysis.getTwiss(2)[2]*twiss_analysis.getTwiss(2)[3])*1000. 
		#emittX = twiss_analysis.getTwiss(0)[3]*1000.0*1000.0	*bunch.getSyncParticle().gamma()*bunch.getSyncParticle().beta()
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f    %5.3f  %5.3f   %5.3f    %5.3f    %5.3f  %5.3f  %7.5f  %10.6f   %8d "%(paramsDict["count"],node.getName(),pos,x_rms,y_rms,z_rms,z_rms_deg,xp_rms,yp_rms,dE_rms,eKin,bunch.getSize())
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
		twiss_analysis.analyzeBunch(bunch)
		x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
		y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
		z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
		z_rms_deg = bunch_gen.getZtoPhaseCoeff(bunch)*z_rms/1000.0		
		xp_rms = math.sqrt(twiss_analysis.getTwiss(0)[2]*twiss_analysis.getTwiss(0)[3])*1000.
		yp_rms = math.sqrt(twiss_analysis.getTwiss(1)[2]*twiss_analysis.getTwiss(1)[3])*1000.
		dE_rms = math.sqrt(twiss_analysis.getTwiss(2)[2]*twiss_analysis.getTwiss(2)[3])*1000. 
		#emittX = twiss_analysis.getTwiss(0)[3]*1000.0*1000.0	*bunch.getSyncParticle().gamma()*bunch.getSyncParticle().beta()
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f    %5.3f  %5.3f   %5.3f    %5.3f    %5.3f  %5.3f  %7.5f  %10.6f   %8d "%(paramsDict["count"],node.getName(),pos,x_rms,y_rms,z_rms,z_rms_deg,xp_rms,yp_rms,dE_rms,eKin,bunch.getSize())
		file_out.write(s +"\n")
		print s	
	
	
#actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

time_start = time.clock()

accLattice.trackBunch(bunch_in, paramsDict = paramsDict, actionContainer = actionContainer)

time_exec = time.clock() - time_start
print "time[sec]=",time_exec

file_out.close()


