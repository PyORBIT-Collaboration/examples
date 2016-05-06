#! /usr/bin/env python

"""
This script will track the bunch through the SNS Linac with an upgraded
for the second target station (STS) SCL linac.

SCL Second Target Station linac does not have cavities in the SCL modules
24 and 26. 

The model includes the following features:
1. At the beginning you have to generate the bunch at SCL:23d:Rg06 exit by 
   un-commenting the part with 'track bunch to the last RF gap in SCL:Cav23d'.
   The script will generate the bunch and will stop. Then comment this part back.
   The rest of the script will read the bunch from the file and will track it 
   to end of HEBT2.
   
2. User will control which cavities are included in acceleration by changing the
   parts of the script with comments 'definitions of the activated cavities with Amp = 1'
   By default the cavities after SCL:23d have Amp=0, and they are not changing the beam.

3. At the end the script will plot the longitudinal phase space by using 
   the Gnuplot package.
   
4. The script includes apertures for all quads

5. There are no overlapping fields quads in the MEBT.

"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from linac import BaseRfGap, MatrixRfGap

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

# import the utilities
from orbit.utils import phaseNearTargetPhase, phaseNearTargetPhaseDeg


def setSynchPhase(bunch_in,accLattice,cav_name,synchPhaseDeg):
	"""
	This function will find the first RF gap phase to get the average
	phase for all gaps equal to the specified synchronous phase.
	"""
	#print "debug start def setSynchPhase(...)"
	b = Bunch()
	rf_cav = accLattice.getRF_Cavity(cav_name)
	rf_gaps = rf_cav.getRF_GapNodes()
	ind_start = accLattice.getNodeIndex(rf_gaps[0])
	ind_stop = accLattice.getNodeIndex(rf_gaps[len(rf_gaps)-1])
	bunch_in.copyEmptyBunchTo(b)
	b = accLattice.trackDesignBunch(b,None,None,-1,ind_start-1)
	e_kin_in = b.getSyncParticle().kinEnergy()
	e_kin_max = 0
	phase_max = 0.
	for ind in range(-180,180):
		phase_deg = ind*1.0
		rf_cav.setPhase(phase_deg*math.pi/180.)
		bunch_in.copyEmptyBunchTo(b)
		b.getSyncParticle().kinEnergy(e_kin_in)
		b = accLattice.trackDesignBunch(b,None,None,ind_start,ind_stop)
		e_kin_out = b.getSyncParticle().kinEnergy()
		if(e_kin_max < e_kin_out):
			e_kin_max = e_kin_out
			phase_max = phase_deg
	cav_phase = phase_max + synchPhaseDeg
	rf_cav.setPhase(cav_phase*math.pi/180.)
	#print "debug cav=",cav_name," phase=",cav_phase," delta_e max=",(e_kin_max-e_kin_in)/1.0e-3
	return cav_phase


#-------------------------------------------------------------------
#          Start of the script
#-------------------------------------------------------------------

random.seed(100)

names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4"]
names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh"]
names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]

#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.02)

#---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"
xml_file_name = "../sns_linac_xml/sns_sts_linac_with_aprt.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print "Linac lattice is ready. L=",accLattice.getLength()

#-----------------------------------------------------
# Set up Space Charge Acc Nodes
#-----------------------------------------------------
from orbit.space_charge.sc3d import setSC3DAccNodes, setUniformEllipsesSCAccNodes
from spacecharge import SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D
sc_path_length_min = 0.02

print "Set up Space Charge nodes. "

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


print "===== Aperture Nodes ======="
aprtNodes = Add_quad_apertures_to_lattice(accLattice)

#for node in aprtNodes:
#	print "aprt=",node.getName()," pos =",node.getPosition()
	
#------ definitions of the activated cavities with Amp = 1 
sclHigh_accSeq = accLattice.getSequence("SCLHigh")
cav_names = []
cav_names = ["SCL:Cav32a","SCL:Cav32b","SCL:Cav32c","SCL:Cav32d",] + cav_names
#cav_names = ["SCL:Cav31a","SCL:Cav31b","SCL:Cav31c","SCL:Cav31d",] + cav_names
#cav_names = ["SCL:Cav30a","SCL:Cav30b","SCL:Cav30c","SCL:Cav30d",] + cav_names
#cav_names = ["SCL:Cav29a","SCL:Cav29b","SCL:Cav29c","SCL:Cav29d",] + cav_names
#cav_names = ["SCL:Cav28a","SCL:Cav28b","SCL:Cav28c","SCL:Cav28d",] + cav_names
#cav_names = ["SCL:Cav27a","SCL:Cav27b","SCL:Cav27c","SCL:Cav27d",] + cav_names
#cav_names = ["SCL:Cav25a","SCL:Cav25b","SCL:Cav25c","SCL:Cav25d",] + cav_names

rf_cavs = []
for cav_name in cav_names:
	rf_cav = sclHigh_accSeq.getRF_Cavity(cav_name)
	rf_cav.setAmp(1.0)
	rf_cavs.append(rf_cav)
	

#------ definitions synchronous phases for the activated cavities with Amp = 1 	
cav_synch_phases = {}

cav_synch_phases["SCL:Cav25a"] = -20.0
cav_synch_phases["SCL:Cav25b"] = -20.0
cav_synch_phases["SCL:Cav25c"] = -20.0
cav_synch_phases["SCL:Cav25d"] = -20.0

cav_synch_phases["SCL:Cav27a"] = -30.0
cav_synch_phases["SCL:Cav27b"] = +30.0
cav_synch_phases["SCL:Cav27c"] = -30.0
cav_synch_phases["SCL:Cav27d"] = -0.0

cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = +30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -0.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = +30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -0.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = +30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -0.0

#------scl32 added -----start----
cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = +30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] =  +0.0
#------scl32 added -----stop----

"""
#------scl31 added -----start----
cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = +30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = +30.0

cav_synch_phases["SCL:Cav32a"] = -35.0
cav_synch_phases["SCL:Cav32b"] = +35.0
cav_synch_phases["SCL:Cav32c"] = -35.0
cav_synch_phases["SCL:Cav32d"] = -20.0
#------scl31 added -----stop----

#------scl30 added -----start----
cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = +30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = +30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = +30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl30 added -----stop----


#------scl29 added -----start----
cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -10.0
cav_synch_phases["SCL:Cav30d"] = +10.0

cav_synch_phases["SCL:Cav31a"] = -0.0
cav_synch_phases["SCL:Cav31b"] = -0.0
cav_synch_phases["SCL:Cav31c"] = -0.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl29 added -----stop----


#------scl28 added -----start----
cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = -30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -30.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl28 added -----stop----

#------scl27 added -----start----
cav_synch_phases["SCL:Cav27a"] = -30.0
cav_synch_phases["SCL:Cav27b"] = -30.0
cav_synch_phases["SCL:Cav27c"] = -30.0
cav_synch_phases["SCL:Cav27d"] = -30.0

cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = -30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -30.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl27 added -----stop----


#------scl25 added -----start----
cav_synch_phases["SCL:Cav25a"] = -30.0
cav_synch_phases["SCL:Cav25b"] = -30.0
cav_synch_phases["SCL:Cav25c"] = -30.0
cav_synch_phases["SCL:Cav25d"] = -30.0

cav_synch_phases["SCL:Cav27a"] = -30.0
cav_synch_phases["SCL:Cav27b"] = -30.0
cav_synch_phases["SCL:Cav27c"] = -30.0
cav_synch_phases["SCL:Cav27d"] = -30.0

cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = -30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -30.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl25 added -----stop----
"""

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


for cav_name in cav_names:
	cav_synch_phase = cav_synch_phases[cav_name]
	setSynchPhase(bunch_in,accLattice,cav_name,cav_synch_phase)

#set up design
accLattice.trackDesignBunch(bunch_in)

#----- the charge of H- is negative, so the phases are shifted by -180.
cavs = accLattice.getRF_Cavities()
for cav in cavs:
	avg_phase = phaseNearTargetPhaseDeg(cav.getAvgGapPhaseDeg()-180.,0.)
	print "debug cav=",cav.getName()," gaps phaseAvg=",avg_phase," cav_amp=",cav.getAmp()
	
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
		alpha_z = twiss_analysis.getTwiss(2)[0]
		#emittX = twiss_analysis.getTwiss(0)[3]*1000.0*1000.0	*bunch.getSyncParticle().gamma()*bunch.getSyncParticle().beta()
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %5d  %35s  %4.5f    %5.3f  %5.3f   %5.3f    %5.3f    %5.3f  %5.3f  %7.5f  %10.6f   %8d "%(paramsDict["count"],node.getName(),pos,x_rms,y_rms,z_rms,z_rms_deg,xp_rms,yp_rms,dE_rms,eKin,bunch.getSize())
		file_out.write(s +"\n")
		print s," a= %5.3f "%alpha_z
	
	
#actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

time_start = time.clock()


rf_cav = accLattice.getRF_Cavity("SCL:Cav23d")
rf_gaps = rf_cav.getRF_GapNodes()
ind_stop = accLattice.getNodeIndex(rf_gaps[len(rf_gaps)-1])


#------------track bunch to the last RF gap in SCL:Cav23d 
#accLattice.trackBunch(bunch_in, paramsDict = paramsDict, actionContainer = actionContainer,index_start = -1, index_stop = ind_stop)
#bunch_in.dumpBunch("bunch_sts_after_scl_23d_rg06.dat")
#sys.exit(1)
#--------------------------------------------------------

bunch_in.deleteAllParticles()
bunch_in.readBunch("bunch_sts_after_scl_23d_rg06.dat")

#set up design
accLattice.trackDesignBunch(bunch_in,None,None,ind_stop+1)

paramsDict["test_pos"] = accLattice.getNodePositionsDict()[rf_gaps[len(rf_gaps)-1]][0]

accLattice.trackBunch(bunch_in, paramsDict = paramsDict, actionContainer = actionContainer,index_start = (ind_stop+1))

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

bunch_in.dumpBunch("bunch_sts_end_hebt2.dat")

z_deg_arr = []
dE_arr = []
for ind in range(bunch_in.getSize()):
	(z_deg,dE) = (bunch_in.z(ind)*bunch_gen.getZtoPhaseCoeff(bunch_in),bunch_in.dE(ind)*1000.)
	z_deg_arr.append(z_deg)
	dE_arr.append(dE)

import Gnuplot
data = Gnuplot.Data(z_deg_arr,dE_arr,with_='p 19', title='Z-dE')
gp = Gnuplot.Gnuplot(persist = 1)
gp('set grid')
gp('set key left')
gp.title("Long. Phase Space")
gp('set xlabel "phi,deg"')
gp('set ylabel "dE, MeV"')
#gp('set pointsize 1.5')
gp.plot(data)
#raw_input('Please press return to stop:\n')

sys.exit(1)

