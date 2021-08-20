##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################
from ROOT import *

#This is completely broken if strippers are included in closed orbit
import math
import sys
import os

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import KickTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.injection import TeapotInjectionNode
from orbit.injection import addTeapotInjectionNode
from orbit.injection import InjectParts
from orbit.injection import JohoTransverse, JohoLongitudinal, SNSESpreadDist
from orbit.kickernodes import XKicker, YKicker
from orbit.kickernodes import rootTWaveform, flatTopWaveform
from orbit.kickernodes import TeapotXKickerNode, TeapotYKickerNode,addTeapotKickerNode
from orbit.foils import TeapotFoilNode, addTeapotFoilNode
from foil import Foil
from orbit.collimation import TeapotCollimatorNode, addTeapotCollimatorNode
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D, Boundary2D
from spacecharge import LSpaceChargeCalc
from orbit.space_charge.sc1d import addLongitudinalSpaceChargeNode, SC1D_AccNode
from orbit.rf_cavities import RFNode, RFLatticeModifications
from spacecharge import Boundary2D

from KevinPython import notRandom
from KevinPython.printNode import Print_Node
import argparse

from orbit.utils.fitting import Solver
from orbit.utils.fitting import Scorer
from orbit.utils.fitting import SolveStopperFactory
from orbit.utils.fitting import VariableProxy
from orbit.utils.fitting import TrialPoint

from orbit.utils.fitting import SimplexSearchAlgorithm
from orbit_utils import Function

from orbit.teapot import GeneralDipoleNoStripSeperateField
from orbit.teapot import GeneralDipoleStripSeperateField
from orbit.teapot import addDipoleStripperNode

from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from orbit.bunch_generators import TwissContainer
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator
from bunch import BunchTwissAnalysis

from Optimizer_Lattice_Class import OptimizerLattice
from MagneticFieldClass import MagneticField
from ConfigureFileClass import ConfigureFileReader

from KevinPython.changeCharge import Change_Charge_Child

class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self,OL_teapot_latt_full,OL_teapot_latt_partial,OL_inject_start,OL_inject_full,beamLatticeFileName,optimizerSettingsFileName):
		Scorer.__init__(self)
		
		#this might be temp. for now all added dipoles should be able to be added to inject lattice
		self.debug_Crash_If_Stripper_Dipole_Ref_Does_Not_Exist=True
		self.OL_teapot_latt=OL_teapot_latt_full
		self.OL_teapot_latt_full=OL_teapot_latt_full
		self.OL_teapot_latt_partial=OL_teapot_latt_partial
		self.OL_inject_start=OL_inject_start
		self.OL_inject_full=OL_inject_full
		#self.doDipoleKickers=doDipoleKickers
		self.addChicaneFieldToStripper=True
		#forDebuggingShouldAlwaysBeTrue
		self.rescaleChicaneFieldInStripper=True

		self.beamLatticeDictionary=ConfigureFileReader(beamLatticeFileName)
		#self.beamLatticeDictionary.printDictionary()
		self.optimizerSettingsDictionary=ConfigureFileReader(optimizerSettingsFileName)
		#self.optimizerSettingsDictionary.printDictionary()	
		
		self.startingNumber=8
		
		#the fudgeFactor is the spacing used when inserting the stripping dipole immediately before or after a node
		#without it rounding errors can place its end or begginning extending past the edge of reference node. So it is shifted the fudgeFactor from the edge
		self.fudgeFactor=.00000001
		#these are the optimizer settings
		
		#these are target values for closed beam at end of lattice
		self.xTarget=0
		self.pxTarget=0
		self.yTarget=0
		self.pyTarget=0
		if self.optimizerSettingsDictionary.hasKey("target_xOffsetClosed"):
			self.xTarget=float(self.optimizerSettingsDictionary.getValue("target_xOffsetClosed"))
		if self.optimizerSettingsDictionary.hasKey("target_pxOffsetClosed"):
			self.pxTarget=float(self.optimizerSettingsDictionary.getValue("target_pxOffsetClosed"))		
		if self.optimizerSettingsDictionary.hasKey("target_yOffsetClosed"):
			self.yTarget=float(self.optimizerSettingsDictionary.getValue("target_yOffsetClosed"))
		if self.optimizerSettingsDictionary.hasKey("target_pyOffsetClosed"):
			self.pyTarget=float(self.optimizerSettingsDictionary.getValue("target_pyOffsetClosed"))
		
		self.targetFullOffsetDifferenceX=0
		self.targetFullOffsetDifferenceY=0
		self.targetPartOffsetDifferenceX=0
		self.targetPartOffsetDifferenceY=0	
		if self.optimizerSettingsDictionary.hasKey("targetFullOffsetDifferenceX"):
			self.targetFullOffsetDifferenceX=float(self.optimizerSettingsDictionary.getValue("targetFullOffsetDifferenceX"))
		if self.optimizerSettingsDictionary.hasKey("targetFullOffsetDifferenceY"):
			self.targetFullOffsetDifferenceY=float(self.optimizerSettingsDictionary.getValue("targetFullOffsetDifferenceY"))		
		if self.optimizerSettingsDictionary.hasKey("targetPartOffsetDifferenceX"):
			self.targetPartOffsetDifferenceX=float(self.optimizerSettingsDictionary.getValue("targetPartOffsetDifferenceX"))
		if self.optimizerSettingsDictionary.hasKey("targetPartOffsetDifferenceY"):
			self.targetPartOffsetDifferenceY=float(self.optimizerSettingsDictionary.getValue("targetPartOffsetDifferenceY"))		
		#these are intial values for closed beam bunch
		self.xInit=0
		self.pxInit=0
		self.yInit=0
		self.pyInit=0
		self.zInit=0
		self.dEInit=0	

		if self.beamLatticeDictionary.hasKey("xOffsetClosed"):
			self.xInit=float(self.beamLatticeDictionary.getValue("xOffsetClosed"))
		if self.beamLatticeDictionary.hasKey("pxOffsetClosed"):
			self.pxInit=float(self.beamLatticeDictionary.getValue("pxOffsetClosed"))		
		if self.beamLatticeDictionary.hasKey("yOffsetClosed"):
			self.yInit=float(self.beamLatticeDictionary.getValue("yOffsetClosed"))
		if self.beamLatticeDictionary.hasKey("pyOffsetClosed"):
			self.pyInit=float(self.beamLatticeDictionary.getValue("pyOffsetClosed"))
			
		#=====Main bunch parameters============
		intensity = 7.8e13
		#turns = 1000.0
		self.turns=1
		turns = self.turns
		#macrosperturn = 260
		macrosperturn = 1
		macrosize = intensity/turns/macrosperturn
		
		#this is the pencil beam bunch for closed beam
		self.b = Bunch()
		self.mass=0.93827231
		self.b.mass(self.mass)
		self.b.macroSize(macrosize)
		#energy = 1.0 #Gev
		e_kin_ini = float(beamLatticeDictionary.getValue("e_kin_ini"))
		self.energy = e_kin_ini #Gev
		self.b.getSyncParticle().kinEnergy(self.energy)
		#self.b.addParticle(0,0,0,0,0,0)
		self.b.addParticle(self.xInit,self.pxInit,self.yInit,self.pyInit,self.zInit,self.dEInit)
		self.initial_chargeClosed=1
		if self.beamLatticeDictionary.hasKey("initial_chargeClosed"):
			self.initial_chargeClosed=int(self.beamLatticeDictionary.getValue("initial_chargeClosed"))
		self.b.charge(self.initial_chargeClosed)
		
		#this is the injected bunch
		self.b2 = Bunch()
		self.initial_chargeInjection=-1
		if self.beamLatticeDictionary.hasKey("initial_chargeInjection"):
			self.initial_chargeInjection=int(self.beamLatticeDictionary.getValue("initial_chargeInjection"))
		self.nPartsInjection=10000
		if self.beamLatticeDictionary.hasKey("nPartsInjection"):
			self.nPartsInjection=int(self.beamLatticeDictionary.getValue("nPartsInjection"))
			
		sp = self.b.getSyncParticle()
		self.beta= sp.beta()
		self.gamma=sp.gamma()
		c=299792458
		self.rigidity= sp.momentum()/(c/math.pow(10.,9))
		
		self.paramsDict = {}
		self.paramsDict["bunch"]= self.b

		self.paramsDict2 = {}
		self.paramsDict2["bunch"]= self.b2
		
		self.numberOfStripperDipoles=-1
		if self.beamLatticeDictionary.hasKey("numberOfStripperDipoles"):
			self.numberOfStripperDipoles=int(self.beamLatticeDictionary.getValue("numberOfStripperDipoles"))
		else:
			print "numberOfStripperDipoles not in beamLattice config file, exiting"
			sys.exit(0)
		self.magneticFields=[]
		for index in range(self.numberOfStripperDipoles):
			if self.beamLatticeDictionary.hasKey("stripper%d"%(index+1)):
				tempFieldClass=None
				#stripperX=["<magnetic field file>","<Dipole Name>",isStripper<True/False>,<Name of node to place it in reference too>,<position with respect to reference node>]
				#<position with respect to reference node> can be "before"(places it right before reference), "after"(places it right after reference), or float [0-1](places it this fraction into reference)				
				tempFieldArray=self.beamLatticeDictionary.getArray("stripper%d"%(index+1))
				tempFieldClass=MagneticField(tempFieldArray[0],tempFieldArray[1],tempFieldArray[2],tempFieldArray[3],tempFieldArray[4],tempFieldArray[5])
				self.magneticFields.append(tempFieldClass)
			else:
				print "stripper%d not in beamLattice config file, exiting"%(index+1)
				sys.exit(0)			
			
			
		#the initial offset of bunches for injection bunch
		self.xOffsetInjection=0.25671
		self.pxOffsetInjection=-.042
		self.yOffsetInjection=0.046
		self.pyOffsetInjection=0
		if self.beamLatticeDictionary.hasKey("xOffsetInjection"):
			self.xOffsetInjection=float(self.beamLatticeDictionary.getValue("xOffsetInjection"))
		if self.beamLatticeDictionary.hasKey("pxOffsetInjection"):
			self.pxOffsetInjection=float(self.beamLatticeDictionary.getValue("pxOffsetInjection"))		
		if self.beamLatticeDictionary.hasKey("yOffsetClosed"):
			self.yOffsetInjection=float(self.beamLatticeDictionary.getValue("yOffsetInjection"))
		if self.beamLatticeDictionary.hasKey("pyOffsetClosed"):
			self.pyOffsetnjection=float(self.beamLatticeDictionary.getValue("pyOffsetInjection"))
			
		intensity = 7.8e13
		turns=1
		e_kin_ini = self.energy # in [GeV]
		mass =  self.mass # in [GeV]
		gamma = (mass + e_kin_ini)/mass
		beta = math.sqrt(gamma*gamma - 1.0)/gamma
		self.macrosize2 = intensity/turns/self.nPartsInjection
		#print "relat. gamma=",gamma
		#print "relat.  beta=",beta
		
		
		#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
		(alphaZ,betaZ,emittZ) = ( float(self.beamLatticeDictionary.getValue("alphaZInjection")), float(self.beamLatticeDictionary.getValue("alphaZInjection")), float(self.beamLatticeDictionary.getValue("alphaZInjection")))
		
		(alphaX,betaX,emittX) = ( float(self.beamLatticeDictionary.getValue("alphaXInjection")), float(self.beamLatticeDictionary.getValue("betaXInjection")), float(self.beamLatticeDictionary.getValue("emittXInjection")))
		(alphaY,betaY,emittY) = ( float(self.beamLatticeDictionary.getValue("alphaYInjection")), float(self.beamLatticeDictionary.getValue("betaYInjection")), float(self.beamLatticeDictionary.getValue("emittYInjection")))
		
		#---make emittances un-normalized XAL units [m*rad]
		emittX = 1.0e-6*emittX/(gamma*beta)
		emittY = 1.0e-6*emittY/(gamma*beta)
		emittZ = 1.0e-6*emittZ/(gamma**3*beta)
		
		#---- long. size in mm
		sizeZ = math.sqrt(emittZ*betaZ)*1.0e+3
		
		#---- transform to pyORBIT emittance[GeV*m]
		emittZ = emittZ*gamma**3*beta**2*mass
		betaZ = betaZ/(gamma**3*beta**2*mass)
		
		twissX = TwissContainer(alphaX,betaX,emittX)
		twissY = TwissContainer(alphaY,betaY,emittY)
		twissZ = TwissContainer(alphaZ,betaZ,emittZ)
		
		#print "Start Bunch Generation."
		self.bunch_gen = SNS_Linac_BunchGenerator(twissX,twissY,twissZ)			
	def getpxOffsetInjection(self):
		return self.pxOffsetInjection
	def getpyOffsetInjection(self):
		return self.pyOffsetInjection	
	def setpxOffsetInjection(self, pxOffset):
		self.pxOffsetInjection=pxOffset
	def setpyOffsetInjection(self,pyOffset):
		self.pyOffsetInjection=pyOffset		
	def getxOffsetInjection(self):
		return self.xOffsetInjection
	def getyOffsetInjection(self):
		return self.yOffsetInjection	
	def setxOffsetInjection(self, xOffset):
		self.xOffsetInjection=xOffset
	def setyOffsetInjection(self,yOffset):
		self.yOffsetInjection=yOffset						
		
	def getpxOffsetClosed(self):
		return self.pxInit
	def getpyOffsetClosed(self):
		return self.pyInit
	def setpxOffsetClosed(self, pxOffset):
		self.pxInit=pxOffset
	def setpyOffsetClosed(self,pyOffset):
		self.pyInit=pyOffset	
	def getNStripper(self):
		return self.numberOfStripperDipoles		
	def resetBunch2(self):

		
		self.b2 = Bunch()

		#generate initial bunch

		self.b2 = self.bunch_gen.getBunch(nParticles = self.nPartsInjection, distributorClass = WaterBagDist3D)
		#bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
		#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)
	
		xOffset=self.xOffsetInjection
		pxOffset=self.pxOffsetInjection
		yOffset=self.yOffsetInjection
		pyOffset=self.pyOffsetInjection
		#if reading bunch from file the offset should already have been added
		
		for i in range(self.b2.getSize()):
			self.b2.x(i,self.b2.x(i)+xOffset)
			self.b2.px(i,self.b2.px(i)+pxOffset)
			self.b2.y(i,self.b2.y(i)+yOffset)
			self.b2.py(i,self.b2.py(i)+pyOffset)			
		self.b2.mass(self.mass) #mass
		self.b2.macroSize(self.macrosize2)
		energy = self.energy # 1.0 #Gev
		self.b2.getSyncParticle().kinEnergy(energy)
		self.b2.charge(self.initial_chargeInjection)
		self.paramsDict2["bunch"]= self.b2
		firstChicaneFail = Bunch()
		firstChicaneFail.charge(-1)
		secondChicaneFail = Bunch()
		secondChicaneFail.charge(0)
		lostbunch = Bunch()
		self.paramsDict2["lostbunch"]=lostbunch
		self.paramsDict2["firstChicaneFail"]=firstChicaneFail
		self.paramsDict2["secondChicaneFail"]= secondChicaneFail		
	#currently assumes bunch has just one macroparticle
	def resetBunch(self):
		self.b.x(0,self.xInit)
		self.b.px(0,self.pxInit)
		self.b.y(0,self.yInit)
		self.b.py(0,self.pyInit)
		self.b.z(0,self.zInit)
		self.b.dE(0,self.dEInit)
	def setStripperLength(self,stripperToScale,scale):
		self.magneticFields[stripperToScale].setStripperLength(self.magneticFields[stripperToScale].getStripperLength()*scale)
	def setStripperAngle(self,stripperToRotate,offset):
		self.magneticFields[stripperToRotate].setFieldDirection(self.magneticFields[stripperToRotate].getFieldDirection()+offset)	
	def setStripperOffset(self,stripperToOffset,offset):
		self.magneticFields[stripperToOffset].setNodePositionOffset(self.magneticFields[stripperToOffset].getNodePositionOffset()+offset)			
	def makeNewInjectLattice(self,full=False):
		inj_latt_start = teapot.TEAPOT_Ring()
		#print "Read MAD."
		#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
		if not full:
			inj_latt_start.readMAD(self.beamLatticeDictionary.getValue("latticeInjection"),"RING")
		else:
			if self.beamLatticeDictionary.hasKey("latticeInjectionFull"):
				inj_latt_start.readMAD(self.beamLatticeDictionary.getValue("latticeInjectionFull"),"RING")
			else:
				inj_latt_start.readMAD(self.beamLatticeDictionary.getValue("latticeInjection"),"RING")
		self.OL_teapot_latt =OptimizerLattice(inj_latt_start)
		if (self.beamLatticeDictionary.hasKey("doDipoleStrippersInjection") and self.beamLatticeDictionary.getValue("doDipoleStrippersInjection")=="True"):
			self.OL_teapot_latt.setDoDipoleStrippers(True)
		self.magneticFields=[]
		for index in range(self.numberOfStripperDipoles):
			if self.beamLatticeDictionary.hasKey("stripper%d"%(index+1)):
				tempFieldClass=None
				#stripperX=["<magnetic field file>","<Dipole Name>",isStripper<True/False>,<Name of node to place it in reference too>,<position with respect to reference node>]
				#<position with respect to reference node> can be "before"(places it right before reference), "after"(places it right after reference), or float [0-1](places it this fraction into reference)				
				tempFieldArray=self.beamLatticeDictionary.getArray("stripper%d"%(index+1))
				tempFieldClass=MagneticField(tempFieldArray[0],tempFieldArray[1],tempFieldArray[2],tempFieldArray[3],tempFieldArray[4],tempFieldArray[5])
				self.magneticFields.append(tempFieldClass)
			else:
				print "stripper%d not in beamLattice config file, exiting"%(index+1)
				sys.exit(0)				
		#print "boom1"
		#self.changeLattice(self.OL_inject_start)
		#self.initChicanes()
		#print "boom2"
		#self.OL_teapot_latt=self.OL_inject_start				
	def setScaleChicane(self,chicaneToScale,scale):
		nodes=self.OL_teapot_latt.getTeapotLattice().getNodes()
		for i in range(len(self.OL_teapot_latt.getChicaneNodes()[chicaneToScale])):
			if self.OL_teapot_latt.getChicaneNodes()[chicaneToScale][i] >=0:
				nodes[self.OL_teapot_latt.getChicaneNodes()[chicaneToScale][i]].setParam("kx", self.OL_teapot_latt.getChicaneNodeStrength()[chicaneToScale][i]*scale)
			
		self.findStrippers()
		for currentIndex in range(self.numberOfStripperDipoles):
			#make sure we are inside the chicane being scaled
			if self.OL_teapot_latt.getDoDipoleStrippers() and self.magneticFields[currentIndex].getIsInsideChicane() and chicaneToScale==self.magneticFields[currentIndex].getChicaneItIsInside() and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
				#print self.magneticFields[currentIndex].getNodeIndex()
				#print nodes[self.magneticFields[currentIndex].getNodeIndex()].getName()
				xkickerField=nodes[self.magneticFields[currentIndex].getNodeIndex()].getChicaneFieldx()
				ykickerField=nodes[self.magneticFields[currentIndex].getNodeIndex()].getChicaneFieldy()	
				nParts=int(self.magneticFields[currentIndex].getNParts())
				lengthStripper=float(self.magneticFields[currentIndex].getStripperLength())
				stepSize=lengthStripper/nParts
				fieldDirection=float(self.magneticFields[currentIndex].getFieldDirection())
				magneticFieldx= Function()
				magneticFieldy= Function()
				for i in range(nParts):
					x = stepSize*i;
					#y = constantField(x)
					y = float(self.magneticFields[currentIndex].getValueOfField(x))
					magneticFieldx.add(x,y*math.cos(fieldDirection)+xkickerField*scale)
					magneticFieldy.add(x,y*math.sin(fieldDirection)+ykickerField*scale)	
				nodes[self.magneticFields[currentIndex].getNodeIndex()].setFunctionMagneticFieldx(magneticFieldx)
				nodes[self.magneticFields[currentIndex].getNodeIndex()].setFunctionMagneticFieldy(magneticFieldy)
				nodes[self.magneticFields[currentIndex].getNodeIndex()].computeFunctions()							
	#find location of stripper nodes
	def findStrippers(self):
		index=0
		for node in self.OL_teapot_latt.getTeapotLattice().getNodes():
			for currentIndex in range(self.numberOfStripperDipoles):
				nodeToFind=self.magneticFields[currentIndex].getNodeName().strip()
				if node.getName().strip() == nodeToFind:
					#print nodeToFind
					#print index
					self.magneticFields[currentIndex].setNodeIndex(index)
					
			index+=1		
	#find location of chicane nodes
	def findChicanes(self):
		#[0-3] are chicane10-13,
		chicaneNodes=[]
		# [0] is drift DB12 and [1] is drift DB23
		driftNodes=[]
		chicane10Nodes=[]
		chicane11Nodes=[]
		chicane12Nodes=[]
		chicane13Nodes=[]
		db12Nodes=[]
		db23Nodes=[]
		index=0
		for node in self.OL_teapot_latt.getTeapotLattice().getNodes():
			if node.getName().strip() == "DH_A10":
				chicane10Nodes.append(index)
			elif node.getName().strip() == "DH_A11":
				chicane11Nodes.append(index)
			elif node.getName().strip() == "DH_A12":
				chicane12Nodes.append(index)
			elif node.getName().strip() == "DH_A13":
				chicane13Nodes.append(index)
			elif node.getName().strip() == "DB12" or node.getName().strip() == "DB12_Laser":
				db12Nodes.append(index)
			elif node.getName().strip() == "DB23":
				db23Nodes.append(index)				
			index+=1
		if len(chicane10Nodes)<1:
			chicane10Nodes.append(-1)
		chicaneNodes.append(chicane10Nodes)
		if len(chicane11Nodes)<1:
			chicane11Nodes.append(-1)		
		chicaneNodes.append(chicane11Nodes)
		if len(chicane12Nodes)<1:
			chicane12Nodes.append(-1)		
		chicaneNodes.append(chicane12Nodes)
		if len(chicane13Nodes)<1:
			chicane13Nodes.append(-1)		
		chicaneNodes.append(chicane13Nodes)
		driftNodes.append(db12Nodes)
		driftNodes.append(db23Nodes)
		self.OL_teapot_latt.setChicaneNodes(chicaneNodes)
		self.OL_teapot_latt.setDriftNodes(driftNodes)		
	def findChicaneStrength(self):	
		self.findChicanes()
		chicaneNodeStrength=[]
		#chicane10Str=[]
		#chicane11Str=[]
		#chicane12Str=[]
		#chicane13Str=[]
		for chicane in self.OL_teapot_latt.getChicaneNodes():
			chicaneStr=[]
			for i in chicane:
				if i>=0:
					chicaneStr.append(self.OL_teapot_latt.getTeapotLattice().getNodes()[i].getParam("kx"))
			chicaneNodeStrength.append(chicaneStr)
			
		self.OL_teapot_latt.setChicaneNodeStrength(chicaneNodeStrength)
	#initialize the chicanes in the teapot lattice.
	def initScorer(self):
		self.initChicanes()
		self.changeLattice(self.OL_teapot_latt_partial)
		self.initChicanes()
		self.changeLattice(self.OL_inject_start)
		self.initChicanes()
		self.changeLattice(self.OL_inject_full)
		self.initChicanes()		
		if False:
			nodes=self.OL_inject_start.getTeapotLattice().getNodes()
			i = 0
			path_length=0
			for node in nodes:
				path_length=path_length+node.getLength()
				print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%self.OL_inject_start.getTeapotLattice().getNodePositionsDict()[node], " path_length= ",path_length
				i=i+1		
				
			functionX=self.OL_inject_start.getTeapotLattice().getNodes()[3].getFunctionMagneticFieldx()
			currentGraph=TGraph()
			for i in range(functionX.getSize()):
				currentGraph.SetPoint(currentGraph.GetN(),float(functionX.x(i)),float(functionX.y(i)))
			theCanvas=TCanvas("TGraph","TGraph",0,0,500,500)
			currentGraph.Draw("AP")   
			theCanvas.Print("textField.png")
			#theCanvas.SaveAs("temp%d.png"%counter)
			theCanvas.Clear()
						
		
	#initialize the chicanes in the teapot lattice.
	def initChicanes(self):
		self.findChicanes()
		#print self.OL_teapot_latt.getChicaneNodes()
		chicanewave = flatTopWaveform(1.0)
		nodes = self.OL_teapot_latt.getTeapotLattice().getNodes()
		for chicane in range(len(self.OL_teapot_latt.getChicaneNodes())):
			for i in range(len(self.OL_teapot_latt.getChicaneNodes()[chicane])):
				if self.OL_teapot_latt.getChicaneNodes()[chicane][i] >=0:
					nodes[self.OL_teapot_latt.getChicaneNodes()[chicane][i]].setWaveform(chicanewave)
					#print chicane, "i= ",i
					#print self.OL_teapot_latt.getChicaneNodeStrength()[chicane][i]
					nodes[self.OL_teapot_latt.getChicaneNodes()[chicane][i]].setParam("kx",self.OL_teapot_latt.getChicaneNodeStrength()[chicane][i])	
		self.addStripperDipoles()
		self.findChicanes()
		self.findChicaneStrength()
		self.findStrippers()				
	#change the lattice
	def changeLattice(self,OL_lattice):
		self.OL_teapot_latt=OL_lattice
	#find and fill the referenceNodeIndex
	#if index=-1 ie none passed it find them all, otherwise it only does it for specified index
	def findReferenceNode(self,index=-1):
		if index==-1:
			for currentIndex in range(self.numberOfStripperDipoles):
				self.magneticFields[currentIndex].setRefNodeIndex(-1)
				counter=0
				nodeToFind=self.magneticFields[currentIndex].getRefNodeName().strip()
				for node in self.OL_teapot_latt.getTeapotLattice().getNodes():
					if node.getName().strip() == nodeToFind:
						self.magneticFields[currentIndex].setRefNodeIndex(counter)			
					counter+=1	
		else:
			counter=0
			nodeToFind=self.magneticFields[index].getRefNodeName().strip()
			self.magneticFields[index].setRefNodeIndex(-1)	
			for node in self.OL_teapot_latt.getTeapotLattice().getNodes():
				if node.getName().strip() == nodeToFind:
					self.magneticFields[index].setRefNodeIndex(counter)			
				counter+=1				
				
	#addFirstStripper
	def addStripperDipoles(self):
		if self.OL_teapot_latt.getDoDipoleStrippers():
			#self.magneticFields=[]
			for index in range(self.numberOfStripperDipoles):
				self.findChicanes()
				self.findReferenceNode(index)		
				
				#if reference node isnt found skip
				if int(self.magneticFields[index].getRefNodeIndex()) >=0:
					#self.debug_Crash_If_Stripper_Dipole_Ref_Does_Not_Exist
					#calculate where to place stripper dipole
					position=-100.
					if self.magneticFields[index].getNodePosition().lower()=="before":
						position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[int(self.magneticFields[index].getRefNodeIndex())]][0]-float(self.magneticFields[index].getStripperLength())-self.fudgeFactor
					elif self.magneticFields[index].getNodePosition().lower()=="after":
						position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[int(self.magneticFields[index].getRefNodeIndex())]][1]+self.fudgeFactor
					else:
						position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[int(self.magneticFields[index].getRefNodeIndex())]][0]+self.OL_teapot_latt.getTeapotLattice().getNodes()[int(self.magneticFields[index].getRefNodeIndex())].getLength()*float(self.magneticFields[index].getNodePosition())
						
					#add possible offset
					position= position+self.magneticFields[index].getNodePositionOffset()
					#check if we are in kicker or drift
					position_start = position
					position_stop = position + float(self.magneticFields[index].getStripperLength())
					#if self.magneticFields[index].getNodePosition().lower()=="before":
						#position_stop=self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[int(self.magneticFields[index].getRefNodeIndex())]][0]
						
					(node_start_ind,node_stop_ind,z,ind) = (-1,-1, 0., 0)
					for nodeCurrent in self.OL_teapot_latt.getTeapotLattice().getNodes():
						if(position_start >= z and position_start <= z + nodeCurrent.getLength()):
							node_start_ind = ind
						if(position_stop > z and position_stop <= z + nodeCurrent.getLength()):
							node_stop_ind = ind
						ind += 1
						z += nodeCurrent.getLength()	
						
					if node_start_ind!=node_stop_ind:
						#the stripping dipole spans more than 1 node
						print "something is going to be broken1"
						if True:
							nodes=self.OL_inject_start.getTeapotLattice().getNodes()
							i = 0
							path_length=0
							for node in nodes:
								path_length=path_length+node.getLength()
								print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%self.OL_inject_start.getTeapotLattice().getNodePositionsDict()[node], " path_length= ",path_length
								i=i+1						
						print position_start
						print position_stop
						sys.exit(0)
					nodeCurrent=self.OL_teapot_latt.getTeapotLattice().getNodes()[node_start_ind]
					#print nodeCurrent.getName()
					xkickerField=0.
					ykickerField=0.
					#nothing to change because no field in drift
					if(isinstance(nodeCurrent,DriftTEAPOT)):
						self.magneticFields[index].setIsInsideChicane(False)
						#print "stripper dipole is in drift"
					#in kicker so add kicker field to stripping field
					elif (isinstance(nodeCurrent,KickTEAPOT)):
						#print "stripper dipole is in kick node"
						if self.addChicaneFieldToStripper:
							self.magneticFields[index].setIsInsideChicane(True)
							if nodeCurrent.getName().strip() == "DH_A10":
								self.magneticFields[index].setChicaneItIsInside(0)
							elif nodeCurrent.getName().strip() == "DH_A11":
								self.magneticFields[index].setChicaneItIsInside(1)
							elif nodeCurrent.getName().strip() == "DH_A12":
								self.magneticFields[index].setChicaneItIsInside(2)
							elif nodeCurrent.getName().strip() == "DH_A13":
								self.magneticFields[index].setChicaneItIsInside(3)							
							#compute field to create kick
							length=nodeCurrent.getLength()
							kx=nodeCurrent.getParam("kx")
							ykickerField=-kx*self.rigidity/length
							ky=nodeCurrent.getParam("ky")
							xkickerField=ky*self.rigidity/length
					#actually creates the field functions
					nParts=int(self.magneticFields[index].getNParts())
					lengthStripper=float(self.magneticFields[index].getStripperLength())
					stepSize=lengthStripper/nParts
					fieldDirection=float(self.magneticFields[index].getFieldDirection())
					magneticFieldx= Function()
					magneticFieldy= Function()
					for i in range(nParts):
						x = stepSize*i;
						#y = constantField(x)
						y = float(self.magneticFields[index].getValueOfField(x))
						magneticFieldx.add(x,y*math.cos(fieldDirection)+xkickerField)
						magneticFieldy.add(x,y*math.sin(fieldDirection)+ykickerField)
					
					if self.magneticFields[index].getIsStripper()=="True":
						myDipole_DH_A11=GeneralDipoleStripSeperateField(magneticFieldx,magneticFieldy,nParts,lengthStripper,self.gamma,self.beta,self.magneticFields[index].getNodeName(),float(self.magneticFields[index].getFixedStrippingLength()))
					else:
						myDipole_DH_A11=GeneralDipoleNoStripSeperateField(magneticFieldx,magneticFieldy,nParts,lengthStripper,self.gamma,self.beta,self.magneticFields[index].getNodeName())
						if self.magneticFields[index].getFoilTest():
							chargeChangeNode=Change_Charge_Child("foil1",1)
							myDipole_DH_A11.addChildNode(chargeChangeNode,AccNode.EXIT)	
					#print "xkickerField=",xkickerField
					myDipole_DH_A11.setChicaneFieldx(xkickerField)
					myDipole_DH_A11.setChicaneFieldy(ykickerField)
					addDipoleStripperNode(self.OL_teapot_latt.getTeapotLattice(),position,myDipole_DH_A11)			
			
	def getScore(self,trialPoint):
		#self.resetBunch2()
		
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()
		x3 = trialPoint.getVariableProxyArr()[3].getValue()
		x4 = trialPoint.getVariableProxyArr()[4].getValue()
		x5 = trialPoint.getVariableProxyArr()[5].getValue()
		x6 = trialPoint.getVariableProxyArr()[6].getValue()
		x7 = trialPoint.getVariableProxyArr()[7].getValue()
	
		self.setpxOffsetClosed(x6)
		self.setpyOffsetClosed(x7)
		self.resetBunch()		
		self.setpxOffsetInjection(x4)
		self.setpyOffsetInjection(x5)
		self.resetBunch2()
		
		self.changeLattice(self.OL_teapot_latt_full)
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)
		
		self.changeLattice(self.OL_teapot_latt_partial)
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)

		self.changeLattice(self.OL_inject_start)
		if (self.optimizerSettingsDictionary.hasKey("modifyAStripper") and self.optimizerSettingsDictionary.getValue("modifyAStripper")=="True"):
			self.makeNewInjectLattice()
			self.OL_inject_start=self.OL_teapot_latt
			self.changeLattice(self.OL_inject_start)
			#if we are changing length need to remake lattice
			for currentIndex in range(self.numberOfStripperDipoles):
				x8 = trialPoint.getVariableProxyArr()[(self.startingNumber+currentIndex*2)].getValue()
				x9 = trialPoint.getVariableProxyArr()[(self.startingNumber+1+currentIndex*2)].getValue()
				x10 = trialPoint.getVariableProxyArr()[(self.startingNumber+self.numberOfStripperDipoles*2+currentIndex)].getValue()				
				#if (self.optimizerSettingsDictionary.hasKey("floatStripperLength%d"%(currentIndex+1)) and self.optimizerSettingsDictionary.getValue("floatStripperLength%d"%(currentIndex+1))=="True"):
					#self.setStripperLength(currentIndex,x8)
				#need to rotate field
				#if (self.optimizerSettingsDictionary.hasKey("floatStripperAngle%d"%(currentIndex+1)) and self.optimizerSettingsDictionary.getValue("floatStripperAngle%d"%(currentIndex+1))=="True"):
					#self.rotateSecondStripperField(currentIndex,x9)
					
				self.setStripperLength(currentIndex,x8)	
				self.setStripperAngle(currentIndex,x9)
				self.setStripperOffset(currentIndex,x10)
			
			self.initChicanes()
			if False:
				nodes=self.OL_inject_start.getTeapotLattice().getNodes()
				i = 0
				path_length=0
				for node in nodes:
					path_length=path_length+node.getLength()
					print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%self.OL_inject_start.getTeapotLattice().getNodePositionsDict()[node], " path_length= ",path_length
					i=i+1				
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)
		
		self.changeLattice(self.OL_inject_full)
		if (self.optimizerSettingsDictionary.hasKey("modifyAStripper") and self.optimizerSettingsDictionary.getValue("modifyAStripper")=="True"):
			self.makeNewInjectLattice(True)
			self.OL_inject_full=self.OL_teapot_latt
			self.changeLattice(self.OL_inject_full)
			#if we are changing length need to remake lattice
			for currentIndex in range(self.numberOfStripperDipoles):
				x8 = trialPoint.getVariableProxyArr()[(self.startingNumber+currentIndex*2)].getValue()
				x9 = trialPoint.getVariableProxyArr()[(self.startingNumber+1+currentIndex*2)].getValue()
				x10 = trialPoint.getVariableProxyArr()[(self.startingNumber+self.numberOfStripperDipoles*2+currentIndex)].getValue()				
				#if (self.optimizerSettingsDictionary.hasKey("floatStripperLength%d"%(currentIndex+1)) and self.optimizerSettingsDictionary.getValue("floatStripperLength%d"%(currentIndex+1))=="True"):
					#self.setStripperLength(currentIndex,x8)
				#need to rotate field
				#if (self.optimizerSettingsDictionary.hasKey("floatStripperAngle%d"%(currentIndex+1)) and self.optimizerSettingsDictionary.getValue("floatStripperAngle%d"%(currentIndex+1))=="True"):
					#self.rotateSecondStripperField(currentIndex,x9)
					
				self.setStripperLength(currentIndex,x8)	
				self.setStripperAngle(currentIndex,x9)
				self.setStripperOffset(currentIndex,x10)
			
			self.initChicanes()
			if False:
				nodes=self.OL_inject_full.getTeapotLattice().getNodes()
				i = 0
				path_length=0
				for node in nodes:
					path_length=path_length+node.getLength()
					print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%self.OL_inject_full.getTeapotLattice().getNodePositionsDict()[node], " path_length= ",path_length
					i=i+1				
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)			

		
		score=0
		self.OL_teapot_latt_full.getTeapotLattice().trackBunch(self.b, self.paramsDict)
		self.OL_inject_full.getTeapotLattice().trackBunch(self.b2, self.paramsDict2)
		twiss_analysis2 = BunchTwissAnalysis()  
		twiss_analysis2.analyzeBunch(self.b2)	
		(xavg2,xpavg2,yavg2,ypavg2)=(twiss_analysis2.getAverage(0),twiss_analysis2.getAverage(1),twiss_analysis2.getAverage(2),twiss_analysis2.getAverage(3))
		if self.optimizerSettingsDictionary.hasKey("usedClosedScore") and self.optimizerSettingsDictionary.getValue("usedClosedScore")=="True":
			score = score +(self.b.x(0)-self.xTarget)**2 + (self.b.px(0)-self.pxTarget)**2 + (self.b.y(0)-self.yTarget)**2+(self.b.py(0)-self.pyTarget)**2
		if self.optimizerSettingsDictionary.hasKey("useFullOffsetDifferenceX") and self.optimizerSettingsDictionary.getValue("useFullOffsetDifferenceX")=="True":
			score = score +(xavg2-self.b.x(0)-self.targetFullOffsetDifferenceX)**2
		if self.optimizerSettingsDictionary.hasKey("useFullOffsetDifferenceY") and self.optimizerSettingsDictionary.getValue("useFullOffsetDifferenceY")=="True":
			score = score +(yavg2-self.b.y(0)-self.targetFullOffsetDifferenceY)**2
		print "xavg2=",xavg2, " self.b.x(0)=",self.b.x(0)," score=",score, " self.targetFullOffsetDifferenceX=",self.targetFullOffsetDifferenceX
		#print "xpavg2=",xpavg2
		self.resetBunch()
		for i in range(self.turns):
			self.OL_teapot_latt_partial.getTeapotLattice().trackBunch(self.b, self.paramsDict)
			
		self.resetBunch2()
		self.OL_inject_start.getTeapotLattice().trackBunch(self.b2, self.paramsDict2)
		twiss_analysis = BunchTwissAnalysis()  
		twiss_analysis.analyzeBunch(self.b2)
		(xavg,xpavg,yavg,ypavg)=(twiss_analysis.getAverage(0),twiss_analysis.getAverage(1),twiss_analysis.getAverage(2),twiss_analysis.getAverage(3))			
		if self.optimizerSettingsDictionary.hasKey("usePartOffsetDifferenceX") and self.optimizerSettingsDictionary.getValue("usePartOffsetDifferenceX")=="True":
			score = score +(xavg-self.b.x(0)-self.targetPartOffsetDifferenceX)**2
		if self.optimizerSettingsDictionary.hasKey("usePartOffsetDifferenceY") and self.optimizerSettingsDictionary.getValue("usePartOffsetDifferenceY")=="True":
			score = score +(yavg-self.b.y(0)-self.targetPartOffsetDifferenceY)**2	
		
		if self.optimizerSettingsDictionary.hasKey("useParallelScore") and self.optimizerSettingsDictionary.getValue("useParallelScore")=="True":
			score = score+(self.b.px(0)-xpavg)**2 +(self.b.py(0)-ypavg)**2
		if self.optimizerSettingsDictionary.hasKey("useParallelScoreY") and self.optimizerSettingsDictionary.getValue("useParallelScoreY")=="True":
			score = score+(self.b.py(0)-ypavg)**2	
		if self.optimizerSettingsDictionary.hasKey("useParallelScoreX") and self.optimizerSettingsDictionary.getValue("useParallelScoreX")=="True":
			score = score+(self.b.px(0)-xpavg)**2		
		#print "self.b.px(0)=",self.b.px(0), " xpavg=",xpavg," score=",score, " x4=",x4
		#print "self.b.py(0)=",self.b.py(0), " ypavg=",ypavg," score=",score, " x5=",x5
		#print "self.b.py(0)=",self.b.py(0), " ypavg=",ypavg," score=",score, "x10=", trialPoint.getVariableProxyArr()[10].getValue()
		#print "xavg2=",xavg2, " self.b.x(0)=",self.b.x(0)," score=",score, " self.targetFullOffsetDifferenceX=",self.targetFullOffsetDifferenceX
		
		return score	
		
print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
#parser.add_argument("--addChicaneFieldToStripper",type=bool, dest='addChicaneFieldToStripper', default=True, help="Include the chicane fields in the stripper if stripper is inside chicane")
parser.add_argument("--outputDirectory", dest='outputDirectory', default="InjectBeam4_ReverseSecond_ChangeOffset", help="Where to put output")
#parser.add_argument("--inputDirectory", dest='chicaneScaleDirectory', default="InjectBeam4_ReverseSecond", help="Where to get chicane scales from")

parser.add_argument("--optimizerConfigFile", dest='optimizerConfigFile', default="OptimizerConfigFiles/DefaultSettings.txt", help="info on optimizer configing")
parser.add_argument("--beamLatticeFile", dest='beamLatticeFile', default="OptimizerConfigFiles/DefaultBeamLattice.txt", help="infoOnInitalBeamAndLattices")
parser.add_argument("--chicaneScaleDirectory", dest='chicaneScaleDirectory', default="InjectBeam3_ChangeOffset", help="Where to get chicane scales from")
args = parser.parse_args()

startingNumber=8

outputDirectory=args.outputDirectory
inputDirectoryChicaneScales=args.chicaneScaleDirectory

beamLatticeDictionary=ConfigureFileReader(args.beamLatticeFile)
#beamLatticeDictionary.printDictionary()
optimizerSettingsDictionary=ConfigureFileReader(args.optimizerConfigFile)
#optimizerSettingsDictionary.printDictionary()
#doDipoleKickers=args.doDipoleKickers
#outputDirectory=args.outputDirectory
#addChicaneFieldToStripper=args.addChicaneFieldToStripper
#inputDirectoryChicaneScales=args.chicaneScaleDirectory

if not os.path.exists(outputDirectory):
	os.mkdir(outputDirectory)
#=====Make a Teapot style lattice======


#the default chicane kick strength array
#chicaneStrengthArray=[-0.041456,0.052434,0.0298523,-0.0398609]
#chicaneNodes=[29,31,34,36]

doDipoleStrippersInjection=False
doDipoleStrippersClosed=False
#print "beamLatticeDictionary.getValue(\"doDipoleStrippers\")",beamLatticeDictionary.getValue("doDipoleStrippers")
#print "hi=%s"%beamLatticeDictionary.getValue("doDipoleStrippers").strip()
#print beamLatticeDictionary.getValue("doDipoleStrippers")
if (beamLatticeDictionary.hasKey("doDipoleStrippersInjection") and beamLatticeDictionary.getValue("doDipoleStrippersInjection")=="True") or (beamLatticeDictionary.hasKey("doDipoleStrippersClosed") and beamLatticeDictionary.getValue("doDipoleStrippersClosed")=="True"):
	#print "lets go"
	if (beamLatticeDictionary.hasKey("doDipoleStrippersInjection") and beamLatticeDictionary.getValue("doDipoleStrippersInjection")=="True"):
		doDipoleStrippersInjection=True
	if (beamLatticeDictionary.hasKey("doDipoleStrippersClosed") and beamLatticeDictionary.getValue("doDipoleStrippersClosed")=="True"):
		doDipoleStrippersClosed=True

latticeInjectionName="none.txt"
latticeInjectionNameFull="none.txt"
latticeClosedCompareToInjectionName="none.txt"
latticeClosedName="none.txt"
if beamLatticeDictionary.hasKey("latticeInjection"):
	latticeInjectionName=beamLatticeDictionary.getValue("latticeInjection")
if beamLatticeDictionary.hasKey("latticeInjectionFull"):
	latticeInjectionNameFull=beamLatticeDictionary.getValue("latticeInjectionFull")
else:
	latticeInjectionNameFull=latticeInjectionName
if beamLatticeDictionary.hasKey("latticeClosedCompareToInjection"):
	latticeClosedCompareToInjectionName=beamLatticeDictionary.getValue("latticeClosedCompareToInjection")
if beamLatticeDictionary.hasKey("latticeClosed"):
	latticeClosedName=beamLatticeDictionary.getValue("latticeClosed")
	

	
inj_latt_start = teapot.TEAPOT_Ring()
print "Read MAD."
#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
inj_latt_start.readMAD(latticeInjectionName,"RING")
print "Lattice=",inj_latt_start.getName()," length [m] =",inj_latt_start.getLength()," nodes=",len(inj_latt_start.getNodes())

OL_inj_latt_start =OptimizerLattice(inj_latt_start)

inj_latt_full = teapot.TEAPOT_Ring()
print "Read MAD."
#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
inj_latt_full.readMAD(latticeInjectionNameFull,"RING")
print "Lattice=",inj_latt_full.getName()," length [m] =",inj_latt_full.getLength()," nodes=",len(inj_latt_full.getNodes())

OL_inj_latt_full =OptimizerLattice(inj_latt_full)
if doDipoleStrippersInjection:
	OL_inj_latt_start.setDoDipoleStrippers(True)
	OL_inj_latt_full.setDoDipoleStrippers(True)

#OL_inj_latt_start.setFirstDipoleIsStripper(True)
#OL_inj_latt_start.setSecondDipoleIsStripper(True)
#OL_inj_latt_start.setSecondDipoleFixedStripLength(.02)

teapot_latt_full = teapot.TEAPOT_Ring()
teapot_latt_full.readMAD(latticeClosedName,"RING")
#print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

teapot_latt_partial = teapot.TEAPOT_Ring()
teapot_latt_partial.readMAD(latticeClosedCompareToInjectionName,"RING")
#print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())


#Turn off injection kickers if present in clsoed lattices

strength_hkicker10 = 0
strength_hkicker13 = strength_hkicker10
strength_hkicker11 = 0
strength_hkicker12 = strength_hkicker11
strength_vkicker10 = 0
strength_vkicker13 = strength_vkicker10
strength_vkicker11 = 0
strength_vkicker12 = strength_vkicker11	

kickerwave = flatTopWaveform(1.0)

nodes = teapot_latt_full.getNodes()
nodes2 = teapot_latt_partial.getNodes()
hkick10 = None
vkick10 = None
hkick11	= None
vkick11 = None
vkick12 = None
hkick12 = None
vkick13 = None
hkick13	= None		
for index in range(len(nodes)):
	if nodes[index].getName().strip() == "IKICKH_A10":
		hkick10=nodes[index]
	elif nodes[index].getName().strip() == "IKICKV_A10":
		vkick10=nodes[index]
	elif nodes[index].getName().strip() == "IKICKH_A11":
		hkick11=nodes[index]
	elif nodes[index].getName().strip() == "IKICKV_A11":
		vkick11=nodes[index]
	elif nodes[index].getName().strip() == "IKICKV_A12":
		vkick12=nodes[index]
	elif nodes[index].getName().strip() == "IKICKH_A12":
		hkick12=nodes[index]
	elif nodes[index].getName().strip() == "IKICKV_A13":
		vkick13=nodes[index]
	elif nodes[index].getName().strip() == "IKICKH_A13":
		hkick13=nodes[index]	
if hkick10 is not None:
	hkick10.setParam("kx", strength_hkicker10)
	hkick10.setWaveform(kickerwave)
if vkick10 is not None:
	vkick10.setParam("ky", strength_vkicker10)
	vkick10.setWaveform(kickerwave)
if hkick11 is not None:
	hkick11.setParam("kx", strength_hkicker11)
	hkick11.setWaveform(kickerwave)
if vkick11 is not None:
	vkick11.setParam("ky", strength_vkicker11)
	vkick11.setWaveform(kickerwave)
if vkick12 is not None:
	vkick12.setParam("ky", strength_vkicker12)
	vkick12.setWaveform(kickerwave)
if hkick12 is not None:
	hkick12.setParam("kx", strength_hkicker12)
	hkick12.setWaveform(kickerwave)
if vkick13 is not None:
	vkick13.setParam("ky", strength_vkicker13)
	vkick13.setWaveform(kickerwave)
if hkick13 is not None:
	hkick13.setParam("kx", strength_hkicker13)
	hkick13.setWaveform(kickerwave)	

hkick10_2 = None
vkick10_2 = None
hkick11_2 = None
vkick11_2 = None
vkick12_2 = None
hkick12_2 = None
vkick13_2 = None
hkick13_2 = None		
for index in range(len(nodes2)):
	if nodes2[index].getName().strip() == "IKICKH_A10":
		hkick10_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKV_A10":
		vkick10_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKH_A11":
		hkick11_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKV_A11":
		vkick11_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKV_A12":
		vkick12_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKH_A12":
		hkick12_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKV_A13":
		vkick13_2=nodes2[index]
	elif nodes2[index].getName().strip() == "IKICKH_A13":
		hkick13_2=nodes2[index]	
if hkick10_2 is not None:
	hkick10_2.setParam("kx", strength_hkicker10)
	hkick10_2.setWaveform(kickerwave)
if vkick10_2 is not None:
	vkick10_2.setParam("ky", strength_vkicker10)
	vkick10_2.setWaveform(kickerwave)
if hkick11_2 is not None:
	hkick11_2.setParam("kx", strength_hkicker11)
	hkick11_2.setWaveform(kickerwave)
if vkick11_2 is not None:
	vkick11_2.setParam("ky", strength_vkicker11)
	vkick11_2.setWaveform(kickerwave)
if vkick12_2 is not None:
	vkick12_2.setParam("ky", strength_vkicker12)
	vkick12_2.setWaveform(kickerwave)
if hkick12_2 is not None:
	hkick12_2.setParam("kx", strength_hkicker12)
	hkick12_2.setWaveform(kickerwave)
if vkick13_2 is not None:
	vkick13_2.setParam("ky", strength_vkicker13)
	vkick13_2.setWaveform(kickerwave)
if hkick13_2 is not None:
	hkick13_2.setParam("kx", strength_hkicker13)
	hkick13_2.setWaveform(kickerwave)	

teapot_latt_full.initialize()
teapot_latt_partial.initialize()


OL_teapot_latt_full =OptimizerLattice(teapot_latt_full)
OL_teapot_latt_partial =OptimizerLattice(teapot_latt_partial)

if doDipoleStrippersClosed:
	OL_teapot_latt_full.setDoDipoleStrippers(True)
	OL_teapot_latt_partial.setDoDipoleStrippers(True)

#scorer = MyScorer(0.005769,0.002069,0.001778,-0.000359,-0.003845,0.000000)
#scorer = MyScorer(0.004334,0.000192,0.001710,-0.000349,-0.004286,0.000000)
#scorer = MyScorer(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000)

scorer = MyScorer(OL_teapot_latt_full,OL_teapot_latt_partial,OL_inj_latt_start,OL_inj_latt_full,args.beamLatticeFile,args.optimizerConfigFile)
scorer.initScorer()	

#searchAlgorithm   = RandomSearchAlgorithm()
searchAlgorithm = SimplexSearchAlgorithm()

#max_time = 0.05
#max_time = 100
max_time = 200
min_score=1E-8
max_iterations=1000

iterationStopper=SolveStopperFactory.maxIterationStopper(max_iterations)
timeStopper = SolveStopperFactory.maxTimeStopper(max_time)
accuracyStopper = SolveStopperFactory.minScoreStopper(min_score)
solverStopper=SolveStopperFactory.comboStopper()
solverStopper.addStopper(iterationStopper)
solverStopper.addStopper(timeStopper)
solverStopper.addStopper(accuracyStopper)
#solverStopper = SolveStopperFactory.maxAccuracyStopper(max_accuracy)

solver = Solver()
solver.setAlgorithm(searchAlgorithm)
solver.setStopper(solverStopper)

useChicaneScaleFile=False
if (optimizerSettingsDictionary.hasKey("readChicaneFile") and optimizerSettingsDictionary.getValue("readChicaneFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale10FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale10FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale11FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale11FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale12FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale12FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale13FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale13FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPXInjectionFromFile") and optimizerSettingsDictionary.getValue("readInitialPXInjectionFromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPYInjectionFromFile") and optimizerSettingsDictionary.getValue("readInitialPYInjectionFromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPXClosedFromFile") and optimizerSettingsDictionary.getValue("readInitialPXClosedFromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPYClosedFromFile") and optimizerSettingsDictionary.getValue("readInitialPYClosedFromFile")=="True"):
	useChicaneScaleFile=True

chicaneScale10=1.
chicaneScale11=1.
chicaneScale12=1.
chicaneScale13=1.
initialPXInjection=float(beamLatticeDictionary.getValue("pxOffsetInjection"))
initialPYInjection=float(beamLatticeDictionary.getValue("pyOffsetInjection"))
initialPXClosed=float(beamLatticeDictionary.getValue("pxOffsetClosed"))
initialPYClosed=float(beamLatticeDictionary.getValue("pyOffsetClosed"))	
initialLengthScaleAngleOffset=[]
initialPositionOffset=[]
for index in range(scorer.getNStripper()):
	pair=[]
	pair.append(1.)
	pair.append(0.)
	initialLengthScaleAngleOffset.append(pair)
	initialPositionOffset.append(0.)		
	
if useChicaneScaleFile:
	openedFile=open("%s/ChicaneScales.txt"%(inputDirectoryChicaneScales),'r')
	line=openedFile.readline()
	print line
	theScales=line.split(",")
	if optimizerSettingsDictionary.hasKey("readChicaneScale10FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale10FromFile")=="True":
		chicaneScale10=float(theScales[0].strip())
	if optimizerSettingsDictionary.hasKey("readChicaneScale11FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale11FromFile")=="True":
		chicaneScale11=float(theScales[1].strip())
	if optimizerSettingsDictionary.hasKey("readChicaneScale12FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale12FromFile")=="True":
		chicaneScale12=float(theScales[2].strip())
	if optimizerSettingsDictionary.hasKey("readChicaneScale13FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale13FromFile")=="True":
		chicaneScale13=float(theScales[3].strip())
	if optimizerSettingsDictionary.hasKey("readInitialPXInjectionFromFile") and optimizerSettingsDictionary.getValue("readInitialPXInjectionFromFile")=="True":
		initialPXInjection=float(theScales[4].strip())		
	if optimizerSettingsDictionary.hasKey("readInitialPYInjectionFromFile") and optimizerSettingsDictionary.getValue("readInitialPYInjectionFromFile")=="True":
		initialPYInjection=float(theScales[5].strip())		
	if optimizerSettingsDictionary.hasKey("readInitialPXClosedFromFile") and optimizerSettingsDictionary.getValue("readInitialPXClosedFromFile")=="True":
		initialPXClosed=float(theScales[6].strip())	
	if optimizerSettingsDictionary.hasKey("readInitialPYClosedFromFile") and optimizerSettingsDictionary.getValue("readInitialPYClosedFromFile")=="True":
		initialPYClosed=float(theScales[7].strip())
	for index in range(scorer.getNStripper()):
		if optimizerSettingsDictionary.hasKey("chicaneScaleFile_useStripperLength%d"%(index+1)) and optimizerSettingsDictionary.getValue("chicaneScaleFile_useStripperLength%d"%(index+1))=="True":
			initialLengthScaleAngleOffset[index][0]=float(theScales[startingNumber+index*2].strip())
		if optimizerSettingsDictionary.hasKey("chicaneScaleFile_useStripperAngle%d"%(index+1)) and optimizerSettingsDictionary.getValue("chicaneScaleFile_useStripperAngle%d"%(index+1))=="True":
			initialLengthScaleAngleOffset[index][1]=float(theScales[startingNumber+1+index*2].strip())	
		if optimizerSettingsDictionary.hasKey("chicaneScaleFile_useStripperOffset%d"%(index+1)) and optimizerSettingsDictionary.getValue("chicaneScaleFile_useStripperOffset%d"%(index+1))=="True":
			initialPositionOffset[index]=float(theScales[startingNumber+scorer.getNStripper()*2+index].strip())				
	openedFile.close()

trialPoint = TrialPoint()
#x0-x3 are chicane10-13 scales
#x4 is px inject offset
#x5 is py inject offset
#x6 is px closed offset
#x7 is py closed offset
#x8 is length of stripper scale
#x9 is angle stripper offset

trialPoint.addVariableProxy(VariableProxy(name = "x0", value = chicaneScale10, step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x1", value = chicaneScale11, step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x2", value = chicaneScale12, step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x3", value = chicaneScale13, step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x4", value = initialPXInjection, step = 0.01))
trialPoint.addVariableProxy(VariableProxy(name = "x5", value = initialPYInjection, step = 0.01))
trialPoint.addVariableProxy(VariableProxy(name = "x6", value = initialPXClosed, step = 0.01))
trialPoint.addVariableProxy(VariableProxy(name = "x7", value = initialPYClosed, step = 0.01))

for index in range(scorer.getNStripper()):
	trialPoint.addVariableProxy(VariableProxy(name = "x%d"%(startingNumber+index*2), value = initialLengthScaleAngleOffset[index][0], step = 0.01))
	trialPoint.addVariableProxy(VariableProxy(name = "x%d"%(startingNumber+1+index*2), value = initialLengthScaleAngleOffset[index][1], step = 0.01))
for index in range(scorer.getNStripper()):	
	trialPoint.addVariableProxy(VariableProxy(name = "x%d"%(startingNumber+scorer.getNStripper()*2+index), value = initialPositionOffset[index], step = 0.01))
x0 = trialPoint.getVariableProxyArr()[0]
x1 = trialPoint.getVariableProxyArr()[1]
x2 = trialPoint.getVariableProxyArr()[2]
x3 = trialPoint.getVariableProxyArr()[3]
x4 = trialPoint.getVariableProxyArr()[4]
x5 = trialPoint.getVariableProxyArr()[5]
x6 = trialPoint.getVariableProxyArr()[6]
x7 = trialPoint.getVariableProxyArr()[7]
for index in range(scorer.getNStripper()):
	x8 = trialPoint.getVariableProxyArr()[(startingNumber+index*2)]
	x9 = trialPoint.getVariableProxyArr()[(startingNumber+1+index*2)]
	x10 = trialPoint.getVariableProxyArr()[(startingNumber+scorer.getNStripper()*2+index)]
	x8.setUseInSolver(False)
	x9.setUseInSolver(False)
	x10.setUseInSolver(False)
	if optimizerSettingsDictionary.hasKey("floatStripperLength%d"%(index+1)) and optimizerSettingsDictionary.getValue("floatStripperLength%d"%(index+1))=="True":
		x8.setUseInSolver(True)
		x8.setLowerLimit(.01)
	if optimizerSettingsDictionary.hasKey("floatStripperAngle%d"%(index+1)) and optimizerSettingsDictionary.getValue("floatStripperAngle%d"%(index+1))=="True":
		x9.setUseInSolver(True)		
	if optimizerSettingsDictionary.hasKey("floatStripperOffset%d"%(index+1)) and optimizerSettingsDictionary.getValue("floatStripperOffset%d"%(index+1))=="True":
		x10.setUseInSolver(True)			

if optimizerSettingsDictionary.hasKey("fixChicaneScale10") and optimizerSettingsDictionary.getValue("fixChicaneScale10")=="True":
	x0.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixChicaneScale11") and optimizerSettingsDictionary.getValue("fixChicaneScale11")=="True":
	x1.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixChicaneScale12") and optimizerSettingsDictionary.getValue("fixChicaneScale12")=="True":
	x2.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixChicaneScale13") and optimizerSettingsDictionary.getValue("fixChicaneScale13")=="True":
	x3.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixInitialPXInjection") and optimizerSettingsDictionary.getValue("fixInitialPXInjection")=="True":
	x4.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixInitialPYInjection") and optimizerSettingsDictionary.getValue("fixInitialPYInjection")=="True":
	x5.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixInitialPXClosed") and optimizerSettingsDictionary.getValue("fixInitialPXClosed")=="True":
	x6.setUseInSolver(False)
if optimizerSettingsDictionary.hasKey("fixInitialPYClosed") and optimizerSettingsDictionary.getValue("fixInitialPYClosed")=="True":
	x7.setUseInSolver(False)


solver.solve(scorer,trialPoint)

print "===== best score ========== fitting time = ", solver.getScoreboard().getRunTime()

bestScore = solver.getScoreboard().getBestScore()	
print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()

trialPoint = solver.getScoreboard().getBestTrialPoint()

print trialPoint.textDesciption()

print "(%f,%f,%f,%f)"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3])
#outputDirectory="WasteBeamClosed"
fileOut=open("%s/ChicaneScales.txt"%(outputDirectory),'w')
fileOut.write("%f,%f,%f,%f,%f,%f,%f,%f"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3],trialPoint.getVariableProxyValuesArr()[4],trialPoint.getVariableProxyValuesArr()[5],trialPoint.getVariableProxyValuesArr()[6],trialPoint.getVariableProxyValuesArr()[7]))
for index in range(scorer.getNStripper()):
	fileOut.write(",%f,%f"%(trialPoint.getVariableProxyValuesArr()[startingNumber+index*2],trialPoint.getVariableProxyValuesArr()[startingNumber+1+index*2]))
for index in range(scorer.getNStripper()):
	fileOut.write(",%f"%(trialPoint.getVariableProxyValuesArr()[startingNumber+scorer.getNStripper()*2+index]))	
fileOut.write("\n")
fileOut.flush() 
fileOut.close() 

