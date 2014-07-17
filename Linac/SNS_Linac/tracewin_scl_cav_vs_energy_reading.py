#! /usr/bin/env python

"""
This script will track the bunch through the SNS Linac with an upgraded
for the second target station (STS) SCL linac 
"""

import sys
import math
import random
import time

from orbit.sns_linac import SimplifiedLinacParser,BaseRF_Gap
from orbit.sns_linac import LinacLatticeFactory, LinacAccLattice
from linac import MatrixRfGap

from bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer


def makePhaseNear(phase, phase0):
	""" It will add or substruct any amount of 360. from phase to get close to phase0 """
	n = int(phase0/360.)
	phase = phase%360.
	min_x = 1.0e+38
	n_min = 0
	for i0 in range(5):
		i = i0 - 3
		d = math.fabs(phase + 360.*(i+n) - phase0)
		if(d < min_x):
			n_min = i
			min_x = d
	return (phase + 360.*(n_min+n))



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
accLattice = lattFactory.getLinacAccLattice(["SCLMed","SCLHigh"])

print "Acc Lattice is ready. "
		
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
rf_cavs = accLattice.getRF_Cavities()
for rf_cav in rf_cavs:
	#print "rf_cav=",rf_cav.getName()," amp=",rf_cav.getAmp()," phase=",(rf_cav.getPhase()-math.pi)*180.0/math.pi
	if(cav_phases_dict.has_key(rf_cav.getName())):
		rf_cav.setParam("",cav_phases_dict[rf_cav.getName()]*math.pi/180.)
	rf_gaps = rf_cav.getRF_GapNodes()
	for rf_gap in rf_gaps:
		#print "      rf_gap=",rf_gap.getName()," E0TL=",rf_gap.getParam("E0TL")," phase=",rf_gap.getParam("gap_phase")*180.0/math.pi
		if(rf_gaps_e0tl_dict.has_key(rf_gap.getName())):
			rf_gap.setParam("E0TL",rf_gaps_e0tl_dict[rf_gap.getName()])

def getRF_Cav(rf_cavs,rf_name):
	for rf_cav in rf_cavs:
		if(rf_cav.getName() == rf_name): return rf_cav 
	return None 

rf_cavs_avg_phases_dict = {}
for rf_cav in rf_cavs:
	rf_cavs_avg_phases_dict[rf_cav] = -20.0

cav = getRF_Cav(rf_cavs,"SCL_RF:Cav01a")
if(cav != None): rf_cavs_avg_phases_dict[cav] += -9.73136

cav = getRF_Cav(rf_cavs,"SCL_RF:Cav01c")
if(cav != None): rf_cavs_avg_phases_dict[cav] += 23.1382

cav = getRF_Cav(rf_cavs,"SCL_RF:Cav02a")
if(cav != None): rf_cavs_avg_phases_dict[cav] += -1.11485

trace_win_fl_in = open("./data/trace_win_results.dat","r")
lns = trace_win_fl_in.readlines()[1:]
trace_win_fl_in.close()

trace_win_pos_eKIn_arr = []
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) > 2):
		pos = float(res_arr[0])
		if(pos > 95.605984):
			trace_win_pos_eKIn_arr.append([pos-95.605984,float(res_arr[1])])
			#print "debug pos=",pos," eKIn=",float(res_arr[1])
		
def get_eKin(pos):
	for pos_ind in range(len(trace_win_pos_eKIn_arr)-1):
		if(pos >= trace_win_pos_eKIn_arr[pos_ind][0] and trace_win_pos_eKIn_arr[pos_ind+1][0] >= pos):
			return trace_win_pos_eKIn_arr[pos_ind][1]

eKin_in = trace_win_pos_eKIn_arr[0][1]

node_pos_dict = 	accLattice.getNodePositionsDict()
#print "debug dict=",node_pos_dict

cav_eKin_dict = {}
for rf_cav_ind in range(len(rf_cavs)-1):
	rf_cav = rf_cavs[rf_cav_ind]
	rf_gaps = rf_cav.getRF_GapNodes()
	rf_cav1 = rf_cavs[rf_cav_ind+1]
	rf_gaps1 = rf_cav.getRF_GapNodes()
	(posBefore, posAfter) = node_pos_dict[rf_gaps[len(rf_gaps)-1]]
	(posBefore1, posAfter1) = node_pos_dict[rf_gaps1[0]]
	pos = (posAfter+posBefore1)/2.0
	eKin_out = get_eKin(pos)
	cav_eKin_dict[rf_cav] = [eKin_in,eKin_out]
	eKin_in = eKin_out
cav_eKin_dict[rf_cavs[len(rf_cavs)-1]]	 = [eKin_in,trace_win_pos_eKIn_arr[len(trace_win_pos_eKIn_arr)-1][1]]

for rf_cav in rf_cavs:
	[eKin_in,eKin_out] = cav_eKin_dict[rf_cav]
	#print "debug cav=",rf_cav.getName()," eKin_in=",eKin_in," eKin_out=",eKin_out
	
	
bunch_init = Bunch()
syncPart = bunch_init.getSyncParticle()
#set H- mass
#self.bunch.mass(0.9382723 + 2*0.000511)
bunch_init.mass(0.939294)
bunch_init.charge(-1.0)
syncPart.kinEnergy(trace_win_pos_eKIn_arr[0][1]*0.001)


def getResultsDict():
	bunch = Bunch()
	bunch_init.copyEmptyBunchTo(bunch)

	#set up design
	accLattice.trackDesignBunch(bunch)

	#track through the lattice START SCL with 95.610 
	rf_gaps_eKin_phases_dict = {}
	paramsDict = {"test_pos":95.605984,"count":0,"rf_gap_dict":rf_gaps_eKin_phases_dict}
	actionContainer = AccActionsContainer("Bunch Tracking")
	
	def action_exit(paramsDict):
		node = paramsDict["node"]
		length = node.getLength()
		pos = paramsDict["test_pos"] + length
		paramsDict["test_pos"] = pos	
		if(isinstance(paramsDict["parentNode"],AccLattice)):
			if(isinstance(node,BaseRF_Gap)):
				bunch_inner = paramsDict["bunch"]
				eKin_out = bunch_inner.getSyncParticle().kinEnergy()*1.0e+3
				phase = makePhaseNear(node.getGapPhase()*180./math.pi-180.,0.)
				rf_gaps_eKin_phases_dict = paramsDict["rf_gap_dict"]
				rf_gaps_eKin_phases_dict[node] = [eKin_out,phase]
				#print "debug eKin out=",eKin_out
	actionContainer.addAction(action_exit, AccActionsContainer.EXIT)
	accLattice.trackBunch(bunch, paramsDict = paramsDict, actionContainer = actionContainer)
	return rf_gaps_eKin_phases_dict
	
rf_gap_e0tl_dict = {}
rf_cav_new_phases_dict = {}
for rf_cav in rf_cavs:
	rf_gaps = rf_cav.getRF_GapNodes()
	[eKin_in,eKin_out] = cav_eKin_dict[rf_cav]
	deltaE = 10.
	deltaPhase = 10.	
	while(math.fabs(deltaE) > 0.01):
		rf_gaps_eKin_phases_dict = getResultsDict()
		eKin_out_new = rf_gaps_eKin_phases_dict[rf_gaps[len(rf_gaps)-1]][0]
		deltaE = eKin_out - eKin_out_new
		coeff = deltaE/50.
		for rf_gap in rf_gaps:
			E0TL = rf_gap.getParam("E0TL")
			rf_gap.setParam("E0TL",E0TL*(1.0+coeff))
			#print "debug E0TL=",E0TL," new E0TL=",E0TL*(1.0+coeff)," deltaE=",deltaE," coeff=",coeff
		while(math.fabs(deltaPhase) > 0.1):
			rf_gaps_eKin_phases_dict = getResultsDict()
			phase_gap_avg = 0.
			for rf_gap in rf_gaps:
				phase_gap_avg += rf_gaps_eKin_phases_dict[rf_gap][1]
			phase_gap_avg /= len(rf_gaps)
			deltaPhase = phase_gap_avg - rf_cavs_avg_phases_dict[rf_cav]
			phase = rf_cav.getPhase() - 0.3*deltaPhase*math.pi/180.
			rf_cav.setPhase(phase)
			#print "rf_cav=",rf_cav.getName()," phase_gap_avg=",phase_gap_avg
		rf_gaps_eKin_phases_dict = getResultsDict()
		eKin_out_new = rf_gaps_eKin_phases_dict[rf_gaps[len(rf_gaps)-1]][0]
		deltaE = eKin_out - eKin_out_new
		phase_gap_avg = 0.
		for rf_gap in rf_gaps:
			phase_gap_avg += rf_gaps_eKin_phases_dict[rf_gap][1]
		phase_gap_avg /= len(rf_gaps)
		deltaPhase = phase_gap_avg - rf_cavs_avg_phases_dict[rf_cav]
	#----------------------------------------------
	rf_gaps_eKin_phases_dict = getResultsDict()
	eKin_out_new = rf_gaps_eKin_phases_dict[rf_gaps[len(rf_gaps)-1]][0]
	deltaE = eKin_out - eKin_out_new
	phase_gap_avg = 0.
	for rf_gap in rf_gaps:
		phase_gap_avg += rf_gaps_eKin_phases_dict[rf_gap][1]
	phase_gap_avg /= len(rf_gaps)
	deltaPhase = phase_gap_avg - rf_cavs_avg_phases_dict[rf_cav]	
	E0TL_avg = 0
	for rf_gap in rf_gaps:
		E0TL_avg += rf_gap.getParam("E0TL")
		rf_gap_e0tl_dict[rf_gap] = rf_gap.getParam("E0TL")
	rf_cav_new_phases_dict[rf_cav] = rf_cav.getPhase()*180./math.pi
	E0TL_avg /= len(rf_gaps)
	print "debug cav=",rf_cav.getName()," new phase=",rf_cav.getPhase()*180/math.pi," deltaE=",deltaE," E0TL_avg=",E0TL_avg

	
fl_out = open("scl_sts_pyorbit_cav_phase.dat","w")
for rf_cav in rf_cavs:
	if(rf_cav_new_phases_dict.has_key(rf_cav)):
		fl_out.write(rf_cav.getName()+" %12.5f "%rf_cav_new_phases_dict[rf_cav]+"\n")
fl_out.close()

fl_out = open("scl_sts_pyorbit_rf_gaps_e0tl.dat","w")
for rf_cav in rf_cavs:
	rf_gaps = rf_cav.getRF_GapNodes()
	for rf_gap in rf_gaps:
		if(rf_gap_e0tl_dict.has_key(rf_gap)):
			fl_out.write(rf_gap.getName()+" %12.10f "%rf_gap_e0tl_dict[rf_gap]+"\n")
fl_out.close()	


