#! /usr/bin/env python

"""
This script will track the bunch through the SNS Linac with an upgraded
for the second target station (STS) SCL linac 
"""

import sys
import math
import random
import time

from orbit.sns_linac import SimplifiedLinacParser
from orbit.sns_linac import LinacLatticeFactory, LinacAccLattice
from linac import BaseRfGap, MatrixRfGap


from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D


from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

random.seed(100)

parser = SimplifiedLinacParser("../SNS_Linac_XML/sns_linac_sts.xml")
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
#accLattice = lattFactory.getLinacAccLattice(["SCLMed","SCLHigh"])
accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4"])

#-----------------------------------------------------
# Set up Space Charge Acc Nodes
#-----------------------------------------------------
from orbit.space_charge.sc3d import setSC3DAccNodes, setUniformEllipsesSCAccNodes
from spacecharge import SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D
sc_path_length_min = 0.02

print "Set up Space Charge nodes. "
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



max_sc_length = 0.
min_sc_length = accLattice.getLength()
for sc_node in space_charge_nodes:
	scL = sc_node.getLengthOfSC()
	if(scL > max_sc_length): max_sc_length = scL
	if(scL < min_sc_length): min_sc_length = scL
print "maximal SC length =",max_sc_length,"  min=",min_sc_length


print "Acc Lattice is ready. "

#-------set up external fields to quads
quads = accLattice.getQuads(accLattice.getSequence("SCLMed"))
quads = accLattice.getQuads()
#for quad in quads:
#	print "quad=",quad.getName()," G[T/m]=",quad.getParam("dB/dr")
	
fl_in = open("./data/sns_linac_quad_fields_matched.dat","r")
lns = fl_in.readlines()
fl_in.close()
quad_name_field_dict = {}
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) == 2):
		quad_name_field_dict[res_arr[0]] = float(res_arr[1])
		
for quad in quads:
	if(	quad_name_field_dict.has_key(quad.getName())):
		field = quad_name_field_dict[quad.getName()]
		quad.setParam("dB/dr",field)
		print "debug quad=",quad.getName()," new field=",field
		
#------- read the SCL cavities phases 
fl_in = open("./data/scl_cavs_phases_sts.dat","r")
lns = fl_in.readlines()
fl_in.close()
cav_phases_dict = {}
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) == 2):
		cav_phases_dict[res_arr[0]] = float(res_arr[1])	
		
#------- read the SCL RF gaps E0TL parameters  (E0TL in MeV) we want GeV 
fl_in = open("./data/scl_rf_gaps_e0tl_sts.dat","r")
lns = fl_in.readlines()
fl_in.close()
rf_gaps_e0tl_dict = {}
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) == 2):
		rf_gaps_e0tl_dict[res_arr[0]] = 0.001*float(res_arr[1])	

#-------correct cavities phases and amplitudes if necessary
#cppGapModel = MatrixRfGap
cppGapModel = BaseRfGap

rf_cavs = accLattice.getRF_Cavities()
for rf_cav in rf_cavs:
	#print "rf_cav=",rf_cav.getName()," amp=",rf_cav.getAmp()," phase=",(rf_cav.getPhase()-math.pi)*180.0/math.pi
	if(cav_phases_dict.has_key(rf_cav.getName())):
		rf_cav.setPhase(cav_phases_dict[rf_cav.getName()]*math.pi/180.)
	rf_gaps = rf_cav.getRF_GapNodes()
	for rf_gap in rf_gaps:
		#print "      rf_gap=",rf_gap.getName()," E0TL=",rf_gap.getParam("E0TL")," phase=",rf_gap.getParam("gap_phase")*180.0/math.pi
		rf_gap.setCppGapModel(cppGapModel())
		if(rf_gaps_e0tl_dict.has_key(rf_gap.getName())):
			rf_gap.setParam("E0TL",rf_gaps_e0tl_dict[rf_gap.getName()])

def getRF_Cav(rf_cavs,rf_name):
	for rf_cav in rf_cavs:
		if(rf_cav.getName() == rf_name): return rf_cav 
	return None 

cav = getRF_Cav(rf_cavs,"MEBT_RF:Bnch03")
if(cav != None):
	amp = cav.getRF_GapNodes()[0].getParam("E0TL")*0.851389
	cav.getRF_GapNodes()[0].setParam("E0TL",amp)
	
cav = getRF_Cav(rf_cavs,"MEBT_RF:Bnch04")
if(cav != None):
	amp = cav.getRF_GapNodes()[0].getParam("E0TL")*1.01944
	cav.getRF_GapNodes()[0].setParam("E0TL",amp)
	
cav = getRF_Cav(rf_cavs,"SCL_RF:Cav01a")
if(cav != None): cav.setPhase(cav.getPhase() -9.73136*math.pi/180.)

cav = getRF_Cav(rf_cavs,"SCL_RF:Cav01c")
if(cav != None): cav.setPhase(cav.getPhase() +23.1382*math.pi/180.)

cav = getRF_Cav(rf_cavs,"SCL_RF:Cav02a")
if(cav != None): cav.setPhase(cav.getPhase() -1.11485*math.pi/180.)

#-----TWISS Parameters at the entrance of MEBT ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.0025 # in [GeV]
mass = 0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta
frequency = 402.5e+6
v_light = 2.99792458e+8  # in [m/sec]

#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
(alphaX,betaX,emittX) = (-1.9620, 0.1831, 0.21)
(alphaY,betaY,emittY) = ( 1.7681, 0.1620, 0.21)
(alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)

alphaZ = -alphaZ

#---make emittances un-normalized XAL units [m*rad]
emittX = 1.0e-6*emittX/(gamma*beta)
emittY = 1.0e-6*emittY/(gamma*beta)
emittZ = 1.0e-6*emittZ/(gamma**3*beta)
print " ========= XAL Twiss ==========="
print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
print " aplha beta emitt[mm*mrad] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

#---- long. size in mm
sizeZ = math.sqrt(emittZ*betaZ)*1.0e+3

#---- transform to pyORBIT emittance[GeV*m]
emittZ = emittZ*gamma**3*beta**2*mass
betaZ = betaZ/(gamma**3*beta**2*mass)

print " ========= PyORBIT Twiss ==========="
print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
print " aplha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

print "Start Bunch Generation."
bunch_gen = SNS_Linac_BunchGenerator(twissX,twissY,twissZ)

#set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

#set the beam peak current in mA
bunch_gen.setBeamCurrent(50.0)

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

eKin = bunch_in.getSyncParticle().kinEnergy()*1.0e+3
twiss_analysis.analyzeBunch(bunch_in)
(alphaX,betaX,emittX) = (twiss_analysis.getTwiss(0)[0],twiss_analysis.getTwiss(0)[1],twiss_analysis.getTwiss(0)[3])
(alphaY,betaY,emittY) = (twiss_analysis.getTwiss(1)[0],twiss_analysis.getTwiss(1)[1],twiss_analysis.getTwiss(1)[3])
(alphaZ,betaZ,emittZ) = (twiss_analysis.getTwiss(2)[0],twiss_analysis.getTwiss(2)[1],twiss_analysis.getTwiss(2)[3])

print " ========= CCL4 exit PyORBIT Twiss =========== eKin=",eKin
print " aplha beta emitt[mm*mrad] X= ( %6.4f , %6.4f , %6.4f ) "%(alphaX,betaX,emittX*1.0e+6)
print " aplha beta emitt[mm*mrad] Y= ( %6.4f , %6.4f , %6.4f ) "%(alphaY,betaY,emittY*1.0e+6)
print " aplha beta emitt[mm*MeV]  Z= ( %6.4f , %6.4f , %6.4f ) "%(alphaZ,betaZ,emittZ*1.0e+6)


