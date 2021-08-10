##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################


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
from ConfigureFileClass import ConfigureFileReader

class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self,OL_teapot_latt_full,OL_teapot_latt_partial,OL_inject_start,beamLatticeFileName,optimizerSettingsFileName,magneticFieldFileName,currentPart=0,nPartsChicane=0,currentPart2=0,nPartsChicane2=0):
		Scorer.__init__(self)
		self.OL_teapot_latt=OL_teapot_latt_full
		self.OL_teapot_latt_full=OL_teapot_latt_full
		self.OL_teapot_latt_partial=OL_teapot_latt_partial
		self.OL_inject_start=OL_inject_start
		self.currentPart=currentPart
		self.nPartsChicane=nPartsChicane
		self.currentPart2=currentPart2
		self.nPartsChicane2=nPartsChicane2		
		#self.doDipoleKickers=doDipoleKickers
		self.addChicaneFieldToStripper=True
		#forDebuggingShouldAlwaysBeTrue
		self.rescaleChicaneFieldInStripper=True

		self.beamLatticeDictionary=ConfigureFileReader(beamLatticeFileName)
		self.beamLatticeDictionary.printDictionary()
		self.optimizerSettingsDictionary=ConfigureFileReader(optimizerSettingsFileName)
		self.optimizerSettingsDictionary.printDictionary()	
		self.magneticFieldDictionary=ConfigureFileReader(magneticFieldFileName)
		self.magneticFieldDictionary.printDictionary()		
		
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
		
		self.theEffLength1=0.03*2
		#theEffLength=0.01
		self.fieldStrength1=1.3
		self.fieldStrengthMin1=.2
		self.cutLength1=0.03
		self.fieldDirection1=math.pi/2.
		self.theEffLength2=0.03*2
		self.fieldStrength2=1.3
		self.fieldStrengthMin2=.2
		self.cutLength2=0.03		
		self.fieldDirection2=math.pi/2.
		
		self.theEffLength1=float(self.magneticFieldDictionary.getValue("stripperLength1"))
		self.fieldStrength1=float(self.magneticFieldDictionary.getValue("stripperStrengthMax1"))
		self.fieldStrengthMin1=float(self.magneticFieldDictionary.getValue("stripperStrengthMin1"))
		self.cutLength1=float(self.magneticFieldDictionary.getValue("cutLength1"))
		if self.magneticFieldDictionary.getValue("fieldDirection1").lower()=="up":
			self.fieldDirection1=math.pi/2.
		elif self.magneticFieldDictionary.getValue("fieldDirection1").lower()=="down":
			self.fieldDirection1=-math.pi/2.
		elif self.magneticFieldDictionary.getValue("fieldDirection1").lower()=="left":
			self.fieldDirection1=0
		elif self.magneticFieldDictionary.getValue("fieldDirection1").lower()=="right":
			self.fieldDirection1=math.pi
		else:
			self.fieldDirection1=float(self.magneticFieldDictionary.getValue("fieldDirection1"))	
			
		self.theEffLength2=float(self.magneticFieldDictionary.getValue("stripperLength2"))
		self.fieldStrength2=float(self.magneticFieldDictionary.getValue("stripperStrengthMax2"))
		self.fieldStrengthMin2=float(self.magneticFieldDictionary.getValue("stripperStrengthMin2"))
		self.cutLength2=float(self.magneticFieldDictionary.getValue("cutLength2"))
		if self.magneticFieldDictionary.getValue("fieldDirection2").lower()=="up":
			self.fieldDirection2=math.pi/2.
		elif self.magneticFieldDictionary.getValue("fieldDirection2").lower()=="down":
			self.fieldDirection2=-math.pi/2.
		elif self.magneticFieldDictionary.getValue("fieldDirection2").lower()=="left":
			self.fieldDirection2=0
		elif self.magneticFieldDictionary.getValue("fieldDirection2").lower()=="right":
			self.fieldDirection2=math.pi
		else:
			self.fieldDirection2=float(self.magneticFieldDictionary.getValue("fieldDirection2"))	
			
		self.n1=1000
		self.maxValue1=self.theEffLength1
		self.step1=self.maxValue1/self.n1

		self.n2_Start=1000
		self.maxValue2_Start=self.theEffLength2
		self.step2_Start=self.maxValue2_Start/self.n2_Start
		self.theEffLength2_Start=self.theEffLength2
		self.fieldDirection2_Start=self.fieldDirection2
		
		self.n2=1000
		self.maxValue2=self.theEffLength2
		self.step2=self.maxValue2/self.n2
		
		self.magneticFieldx= Function()
		self.magneticFieldy= Function()	
		self.magneticFieldx2= Function()
		self.magneticFieldy2= Function()	
		
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
	def pieceWiseField1(self,x):
		if x<self.cutLength1 :
			return (self.fieldStrength1-self.fieldStrengthMin1)/self.cutLength1*x+self.fieldStrengthMin1
		elif x>=self.cutLength1:
			return self.fieldStrength1
		pass			
	def pieceWiseField2(self,x):
		if x<self.cutLength2 :
			return (self.fieldStrength2-self.fieldStrengthMin2)/self.cutLength2*x+self.fieldStrengthMin2
		elif x>=self.cutLength2:
			return self.fieldStrength2
		pass		
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
	def setSecondStripperLength(self,scale):
		self.n2=int(self.n2_Start*scale)
		self.theEffLength2=self.theEffLength2_Start*scale
		self.maxValue2=self.theEffLength2
		self.step2=self.maxValue2/self.n2
	def makeNewInjectLattice(self):
		inj_latt_start = teapot.TEAPOT_Ring()
		print "Read MAD."
		#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
		inj_latt_start.readMAD(self.beamLatticeDictionary.getValue("latticeInjection"),"RING")		
		self.OL_inject_start =OptimizerLattice(inj_latt_start)
		if (self.beamLatticeDictionary.hasKey("doDipoleStrippersInjection") and self.beamLatticeDictionary.getValue("doDipoleStrippersInjection")=="True"):
			self.OL_inject_start.setDoDipoleStrippers(True)
			if self.beamLatticeDictionary.hasKey("firstStripperIsStripperInjection") and self.beamLatticeDictionary.getValue("firstStripperIsStripperInjection")=="True":
				self.OL_inject_start.setFirstDipoleIsStripper(True)
			if self.beamLatticeDictionary.hasKey("secondStripperIsStripperInjection") and self.beamLatticeDictionary.getValue("secondStripperIsStripperInjection")=="True":
				self.OL_inject_start.setSecondDipoleIsStripper(True)	
			if self.beamLatticeDictionary.hasKey("secondStrippingLengthInjection"):
				self.OL_inject_start.setSecondDipoleFixedStripLength(float(self.beamLatticeDictionary.getValue("secondStrippingLengthInjection")))
		print "boom1"
		self.changeLattice(self.OL_inject_start)
		self.initChicanes()
		print "boom2"
		self.OL_teapot_latt=self.OL_inject_start
	def rotateSecondStripperField(self,scale):
		self.fieldDirection2=self.fieldDirection2_Start*scale
		nodes=self.OL_teapot_latt.getTeapotLattice().getNodes()
		
		self.magneticFieldx2= Function()
		self.magneticFieldy2= Function()	
		for i in range(self.n2):
			x = self.step2*i;
			#y = constantField(x)
			y = self.pieceWiseField2(x)
			self.magneticFieldx2.add(x,y*math.cos(self.fieldDirection2))
			self.magneticFieldy2.add(x,y*math.sin(self.fieldDirection2))	
		nodes[self.OL_teapot_latt.getSecondDipoleNode()].setFunctionMagneticFieldx(self.magneticFieldx2)
		nodes[self.OL_teapot_latt.getSecondDipoleNode()].setFunctionMagneticFieldy(self.magneticFieldy2)
		nodes[self.OL_teapot_latt.getSecondDipoleNode()].computeFunctions()			
	def setScaleChicane(self,chicaneToScale,scale):
		nodes=self.OL_teapot_latt.getTeapotLattice().getNodes()
		for i in range(len(self.OL_teapot_latt.getChicaneNodes()[chicaneToScale])):
			if self.OL_teapot_latt.getChicaneNodes()[chicaneToScale][i] >=0:
				nodes[self.OL_teapot_latt.getChicaneNodes()[chicaneToScale][i]].setParam("kx", self.OL_teapot_latt.getChicaneNodeStrength()[chicaneToScale][i]*scale)
			
		#need to recompute dipole field if stripper is in chicane field because chicane field has been changed
		if chicaneToScale==1 and self.OL_teapot_latt.getFirstDipoleInChicane() and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
			print "we are recomputing stripper1"
			self.magneticFieldx= Function()
			self.magneticFieldy= Function()				
			xkickerField=nodes[self.OL_teapot_latt.getFirstDipoleNode()].getChicaneFieldx()
			ykickerField=nodes[self.OL_teapot_latt.getFirstDipoleNode()].getChicaneFieldy()
			for i in range(self.n1):
				x = self.step1*i;
				#y = constantField(x)
				y = self.pieceWiseField1(x)
				self.magneticFieldx.add(x,y*math.cos(self.fieldDirection1)+xkickerField*scale)
				self.magneticFieldy.add(x,y*math.sin(self.fieldDirection1)+ykickerField*scale)	
				
			nodes[self.OL_teapot_latt.getFirstDipoleNode()].setFunctionMagneticFieldx(self.magneticFieldx)
			nodes[self.OL_teapot_latt.getFirstDipoleNode()].setFunctionMagneticFieldy(self.magneticFieldy)
			nodes[self.OL_teapot_latt.getFirstDipoleNode()].computeFunctions()
			
		elif chicaneToScale==2 and self.OL_teapot_latt.getSecondDipoleInChicane() and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
			print "we are recomputing stripper2"
			self.magneticFieldx2= Function()
			self.magneticFieldy2= Function()				
			xkickerField=nodes[self.OL_teapot_latt.getSecondDipoleNode()].getChicaneFieldx()
			ykickerField=nodes[self.OL_teapot_latt.getSecondDipoleNode()].getChicaneFieldy()
			for i in range(self.n2):
				x = self.step2*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx2.add(x,y*math.cos(self.fieldDirection2)+xkickerField*scale)
				self.magneticFieldy2.add(x,y*math.sin(self.fieldDirection2)+ykickerField*scale)	
				
			nodes[self.OL_teapot_latt.getSecondDipoleNode()].setFunctionMagneticFieldx(self.magneticFieldx2)
			nodes[self.OL_teapot_latt.getSecondDipoleNode()].setFunctionMagneticFieldy(self.magneticFieldy2)
			nodes[self.OL_teapot_latt.getSecondDipoleNode()].computeFunctions()			

	#find location of stripper nodes
	def findStrippers(self):
		index=0
		for node in self.OL_teapot_latt.getTeapotLattice().getNodes():
			if node.getName().strip() == "Dipole_DH_A11":
				self.OL_teapot_latt.setFirstDipoleNode(index)
			elif node.getName().strip() == "Dipole_DH_A12":
				self.OL_teapot_latt.setSecondDipoleNode(index)			
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
		
	#initialize the chicanes in the teapot lattice.
	def initChicanes(self):
		self.findChicanes()
		print self.OL_teapot_latt.getChicaneNodes()
		chicanewave = flatTopWaveform(1.0)
		nodes = self.OL_teapot_latt.getTeapotLattice().getNodes()
		for chicane in range(len(self.OL_teapot_latt.getChicaneNodes())):
			for i in range(len(self.OL_teapot_latt.getChicaneNodes()[chicane])):
				if self.OL_teapot_latt.getChicaneNodes()[chicane][i] >=0:
					nodes[self.OL_teapot_latt.getChicaneNodes()[chicane][i]].setWaveform(chicanewave)
					print chicane, "i= ",i
					print self.OL_teapot_latt.getChicaneNodeStrength()[chicane][i]
					nodes[self.OL_teapot_latt.getChicaneNodes()[chicane][i]].setParam("kx",self.OL_teapot_latt.getChicaneNodeStrength()[chicane][i])	
		self.addFirstStripperDipole()
		self.addSecondStripperDipole()
		self.findChicanes()
		self.findChicaneStrength()
		self.findStrippers()				
	#change the lattice
	def changeLattice(self,OL_lattice):
		self.OL_teapot_latt=OL_lattice
	#addFirstStripper
	def addFirstStripperDipole(self):
		if self.OL_teapot_latt.getDoDipoleStrippers() and self.OL_teapot_latt.getChicaneNodes()[1][0] >=0:
			self.findChicanes()
			#calculate where to place 1st stripper dipole
			position=-100.
			if self.currentPart==-1:
				#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[1][0]]][0]-self.theEffLength
				#print self.OL_teapot_latt.getDriftNodes()
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getDriftNodes()[0][0]]][1]-self.theEffLength1
			elif self.currentPart==0:
				print "OL_teapot_latt.getChicaneNodes()[1][0]= ",self.OL_teapot_latt.getChicaneNodes()[1][0]
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[1][0]]][0]
			elif self.currentPart is self.nPartsChicane:
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[1][0]]][1]
			else :
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[1][0]]][0]+self.OL_teapot_latt.getTeapotLattice().getNodes()[OL_teapot_latt.getChicaneNodes()[1][0]].getLength()*self.currentPart/self.nPartsChicane
								
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + self.theEffLength1
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
				print position_start
				print position_stop
				sys.exit(0)
			nodeCurrent=self.OL_teapot_latt.getTeapotLattice().getNodes()[node_start_ind]
			xkickerField=0.
			ykickerField=0.
			#nothing to change because no field in drift
			if(isinstance(nodeCurrent,DriftTEAPOT)):
				self.OL_teapot_latt.setFirstDipoleInChicane(False)
				print "stripper dipole is in drift"
			#in kicker so add kicker field to stripping field
			elif (isinstance(nodeCurrent,KickTEAPOT)):
				print "stripper dipole is in kick node"
				if self.addChicaneFieldToStripper:
					self.OL_teapot_latt.setFirstDipoleInChicane(True)
					#compute field to create kick
					length=nodeCurrent.getLength()
					kx=nodeCurrent.getParam("kx")
					ykickerField=-kx*self.rigidity/length
					ky=nodeCurrent.getParam("ky")
					xkickerField=ky*self.rigidity/length
			#actually creates the field functions
		
			for i in range(self.n1):
				x = self.step1*i;
				#y = constantField(x)
				y = self.pieceWiseField1(x)
				self.magneticFieldx.add(x,y*math.cos(self.fieldDirection1)+xkickerField)
				self.magneticFieldy.add(x,y*math.sin(self.fieldDirection1)+ykickerField)
			
			if self.OL_teapot_latt.getFirstDipoleIsStripper():
				myDipole_DH_A11=GeneralDipoleStripSeperateField(self.magneticFieldx,self.magneticFieldy,self.n1,self.maxValue1,self.gamma,self.beta,"Dipole_DH_A11",self.OL_teapot_latt.getFirstDipoleFixedStripLength())
			else:
				myDipole_DH_A11=GeneralDipoleNoStripSeperateField(self.magneticFieldx,self.magneticFieldy,self.n1,self.maxValue1,self.gamma,self.beta,"Dipole_DH_A11")
			#print "xkickerField=",xkickerField
			myDipole_DH_A11.setChicaneFieldx(xkickerField)
			myDipole_DH_A11.setChicaneFieldy(ykickerField)
			addDipoleStripperNode(self.OL_teapot_latt.getTeapotLattice(),position,myDipole_DH_A11)	
			
	#addSecondStripper
	def addSecondStripperDipole(self):
		if self.OL_teapot_latt.getDoDipoleStrippers() and self.OL_teapot_latt.getChicaneNodes()[2][0] >=0:
			self.findChicanes()
			#calculate where to place 1st stripper dipole
			position=-100.
			#place second stripper 5/6 of the way into chicane3/12. temporary position for consistency
			#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]+self.teapot_latt.getNodes()[self.chicaneNodes[2][0]].getLength()*5./6.
			if self.currentPart2==-1:
				#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]-self.theEffLength
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getDriftNodes()[1][0]]][1]-self.theEffLength2
			elif self.currentPart2==0:
				print "self.chicaneNodes[1][0]= ",self.OL_teapot_latt.getChicaneNodes()[1][0]
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]]][0]
			elif self.currentPart2 is self.nPartsChicane2:
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]]][1]
			else :
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]]][0]+self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]].getLength()*self.currentPart2/self.nPartsChicane2			
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + self.theEffLength2
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
				print "something is going to be broken2"
				sys.exit(0)
			nodeCurrent=self.OL_teapot_latt.getTeapotLattice().getNodes()[node_start_ind]
			xkickerField=0.
			ykickerField=0.
			#nothing to change because no field in drift
			if(isinstance(nodeCurrent,DriftTEAPOT)):
				self.OL_teapot_latt.setSecondDipoleInChicane(False)
				print "stripper dipole is in drift"
			#in kicker so add kicker field to stripping field
			elif (isinstance(nodeCurrent,KickTEAPOT)):
				print "stripper dipole is in kick node"
				if self.addChicaneFieldToStripper:
					self.OL_teapot_latt.setSecondDipoleInChicane(True)
					#compute field to create kick
					length=nodeCurrent.getLength()
					kx=nodeCurrent.getParam("kx")
					ykickerField=-kx*self.rigidity/length
					ky=nodeCurrent.getParam("ky")
					xkickerField=ky*self.rigidity/length	
			#actually creates the field functions
		
			for i in range(self.n2):
				x = self.step2*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx2.add(x,y*math.cos(self.fieldDirection2)+xkickerField)
				self.magneticFieldy2.add(x,y*math.sin(self.fieldDirection2)+ykickerField)
			
			if self.OL_teapot_latt.getSecondDipoleIsStripper():
				myDipole_DH_A12=GeneralDipoleStripSeperateField(self.magneticFieldx2,self.magneticFieldy2,self.n2,self.maxValue2,self.gamma,self.beta,"Dipole_DH_A12",self.OL_teapot_latt.getSecondDipoleFixedStripLength())
			else:
				myDipole_DH_A12=GeneralDipoleNoStripSeperateField(self.magneticFieldx2,self.magneticFieldy2,self.n2,self.maxValue2,self.gamma,self.beta,"Dipole_DH_A12")
			myDipole_DH_A12.setChicaneFieldx(xkickerField)
			myDipole_DH_A12.setChicaneFieldy(ykickerField)			
			addDipoleStripperNode(self.OL_teapot_latt.getTeapotLattice(),position,myDipole_DH_A12)				
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
		x8 = trialPoint.getVariableProxyArr()[8].getValue()
		x9 = trialPoint.getVariableProxyArr()[9].getValue()
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
		#if we are changing length need to remake lattice
		if (optimizerSettingsDictionary.hasKey("floatSecondStripperLength") and optimizerSettingsDictionary.getValue("floatSecondStripperLength")=="True"):
			self.setSecondStripperLength(x8)
			self.makeNewInjectLattice()
		#self.changeLattice(self.OL_inject_start)
		#need to rotate field
		if (optimizerSettingsDictionary.hasKey("floatSecondStripperAngle") and optimizerSettingsDictionary.getValue("floatSecondStripperAngle")=="True"):
			self.rotateSecondStripperField(x9)
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)
			
		
		score=0
		self.OL_teapot_latt_full.getTeapotLattice().trackBunch(self.b, self.paramsDict)
		if optimizerSettingsDictionary.hasKey("usedClosedScore") and optimizerSettingsDictionary.getValue("usedClosedScore")=="True":
			score = score +(self.b.x(0)-self.xTarget)**2 + (self.b.px(0)-self.pxTarget)**2 + (self.b.y(0)-self.yTarget)**2+(self.b.py(0)-self.pyTarget)**2
		self.resetBunch()
		for i in range(self.turns):
			self.OL_teapot_latt_partial.getTeapotLattice().trackBunch(self.b, self.paramsDict)
		self.OL_inject_start.getTeapotLattice().trackBunch(self.b2, self.paramsDict2)
		#self.OL_inject_end.getTeapotLattice().trackBunch(self.b2, self.paramsDict2)
		twiss_analysis = BunchTwissAnalysis()  
		twiss_analysis.analyzeBunch(self.b2)
		(xavg,xpavg,yavg,ypavg)=(twiss_analysis.getAverage(0),twiss_analysis.getAverage(1),twiss_analysis.getAverage(2),twiss_analysis.getAverage(3))
		if optimizerSettingsDictionary.hasKey("useParallelScore") and optimizerSettingsDictionary.getValue("useParallelScore")=="True":
			score = score+(self.b.px(0)-xpavg)**2 +(self.b.py(0)-ypavg)**2
		#print "self.b.px(0)=",self.b.px(0), " xpavg=",xpavg," score=",score, " x4=",x4
		#print "self.b.py(0)=",self.b.py(0), " ypavg=",ypavg," score=",score, " x5=",x5
		return score	
		
print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
#parser.add_argument("--addChicaneFieldToStripper",type=bool, dest='addChicaneFieldToStripper', default=True, help="Include the chicane fields in the stripper if stripper is inside chicane")
parser.add_argument("--outputDirectory", dest='outputDirectory', default="InjectBeam4_ReverseSecond_ChangeOffset", help="Where to put output")
#parser.add_argument("--inputDirectory", dest='chicaneScaleDirectory', default="InjectBeam4_ReverseSecond", help="Where to get chicane scales from")

parser.add_argument("--magneticFieldFile", dest='magneticFieldFile', default="MagneticFieldFiles/magneticFieldUpUp.txt", help="infoOnMagneticField")
parser.add_argument("--optimizerConfigFile", dest='optimizerConfigFile', default="OptimizerConfigFiles/DefaultSettings.txt", help="info on optimizer configing")
parser.add_argument("--beamLatticeFile", dest='beamLatticeFile', default="OptimizerConfigFiles/DefaultBeamLattice.txt", help="infoOnInitalBeamAndLattices")
parser.add_argument("--chicaneScaleDirectory", dest='chicaneScaleDirectory', default="InjectBeam3_ChangeOffset", help="Where to get chicane scales from")
args = parser.parse_args()

outputDirectory=args.outputDirectory
inputDirectoryChicaneScales=args.chicaneScaleDirectory

beamLatticeDictionary=ConfigureFileReader(args.beamLatticeFile)
beamLatticeDictionary.printDictionary()
optimizerSettingsDictionary=ConfigureFileReader(args.optimizerConfigFile)
optimizerSettingsDictionary.printDictionary()
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
#this sets how to divide up chicane2/11 in terms of where 1st stripper is placed.
nPartsChicane=1
nPartsChicane2=0
stripperPositionArray=["0"]
stripperPositionArray2=["0"]
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
	nPartsChicane=int(beamLatticeDictionary.getValue("firstStripperPositionMax"))
	nPartsChicane2=int(beamLatticeDictionary.getValue("secondStripperPositionMax"))	
	stripperPositionArray=beamLatticeDictionary.getArray("firstStripperPositionArray")
	stripperPositionArray2=beamLatticeDictionary.getArray("secondStripperPositionArray")

latticeInjectionName="none.txt"
latticeClosedCompareToInjectionName="none.txt"
latticeClosedName="none.txt"
if beamLatticeDictionary.hasKey("latticeInjection"):
	latticeInjectionName=beamLatticeDictionary.getValue("latticeInjection")
if beamLatticeDictionary.hasKey("latticeClosedCompareToInjection"):
	latticeClosedCompareToInjectionName=beamLatticeDictionary.getValue("latticeClosedCompareToInjection")
if beamLatticeDictionary.hasKey("latticeClosed"):
	latticeClosedName=beamLatticeDictionary.getValue("latticeClosed")
	
useChicaneScaleFile=False
if (optimizerSettingsDictionary.hasKey("readChicaneScale10FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale10FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale11FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale11FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale12FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale12FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readChicaneScale13FromFile") and optimizerSettingsDictionary.getValue("readChicaneScale13FromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPXInjectionFromFile") and optimizerSettingsDictionary.getValue("readInitialPXInjectionFromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPYInjectionFromFile") and optimizerSettingsDictionary.getValue("readInitialPYInjectionFromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPXClosedFromFile") and optimizerSettingsDictionary.getValue("readInitialPXClosedFromFile")=="True") or (optimizerSettingsDictionary.hasKey("readInitialPYClosedFromFile") and optimizerSettingsDictionary.getValue("readInitialPYClosedFromFile")=="True"):
	useChicaneScaleFile=True
#for currentPart in range(-1,nPartsChicane+1):
for currentPart in stripperPositionArray:
	currentPart=int(currentPart)
	for currentPart2 in stripperPositionArray2:
		currentPart2=int(currentPart2)
		chicaneScale10=1.
		chicaneScale11=1.
		chicaneScale12=1.
		chicaneScale13=1.
		initialPXInjection=float(beamLatticeDictionary.getValue("pxOffsetInjection"))
		initialPYInjection=float(beamLatticeDictionary.getValue("pyOffsetInjection"))
		initialPXClosed=float(beamLatticeDictionary.getValue("pxOffsetClosed"))
		initialPYClosed=float(beamLatticeDictionary.getValue("pyOffsetClosed"))		
		if useChicaneScaleFile:
			openedFile=open("%s/ChicaneScales_%d_%d_%d_%d.txt"%(inputDirectoryChicaneScales,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'r')
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
			openedFile.close()
			
		inj_latt_start = teapot.TEAPOT_Ring()
		print "Read MAD."
		#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
		inj_latt_start.readMAD(latticeInjectionName,"RING")
		print "Lattice=",inj_latt_start.getName()," length [m] =",inj_latt_start.getLength()," nodes=",len(inj_latt_start.getNodes())
		
		OL_inj_latt_start =OptimizerLattice(inj_latt_start)
		if doDipoleStrippersInjection:
			OL_inj_latt_start.setDoDipoleStrippers(True)
			if beamLatticeDictionary.hasKey("firstStripperIsStripperInjection") and beamLatticeDictionary.getValue("firstStripperIsStripperInjection")=="True":
				OL_inj_latt_start.setFirstDipoleIsStripper(True)
			if beamLatticeDictionary.hasKey("secondStripperIsStripperInjection") and beamLatticeDictionary.getValue("secondStripperIsStripperInjection")=="True":
				OL_inj_latt_start.setSecondDipoleIsStripper(True)	
			if beamLatticeDictionary.hasKey("secondStrippingLengthInjection"):
				OL_inj_latt_start.setSecondDipoleFixedStripLength(float(beamLatticeDictionary.getValue("secondStrippingLengthInjection")))
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
		
		scorer = MyScorer(OL_teapot_latt_full,OL_teapot_latt_partial,OL_inj_latt_start,args.beamLatticeFile,args.optimizerConfigFile,args.magneticFieldFile,currentPart,nPartsChicane,currentPart2,nPartsChicane2)
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
		
		trialPoint = TrialPoint()
		#x0-x3 are chicane10-13 scales
		#x4 is px inject offset
		#x5 is py inject offset
		#x6 is px closed offset
		#x7 is py closed offset
		#x8 is length of stripper scale
		#x9 is angle stripper scale
		
		trialPoint.addVariableProxy(VariableProxy(name = "x0", value = chicaneScale10, step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x1", value = chicaneScale11, step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x2", value = chicaneScale12, step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x3", value = chicaneScale13, step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x4", value = initialPXInjection, step = 0.01))
		trialPoint.addVariableProxy(VariableProxy(name = "x5", value = initialPYInjection, step = 0.01))
		trialPoint.addVariableProxy(VariableProxy(name = "x6", value = initialPXClosed, step = 0.01))
		trialPoint.addVariableProxy(VariableProxy(name = "x7", value = initialPYClosed, step = 0.01))
		trialPoint.addVariableProxy(VariableProxy(name = "x8", value = 1., step = 0.01))
		trialPoint.addVariableProxy(VariableProxy(name = "x9", value = 1., step = 0.01))		
		x0 = trialPoint.getVariableProxyArr()[0]
		x1 = trialPoint.getVariableProxyArr()[1]
		x2 = trialPoint.getVariableProxyArr()[2]
		x3 = trialPoint.getVariableProxyArr()[3]
		x4 = trialPoint.getVariableProxyArr()[4]
		x5 = trialPoint.getVariableProxyArr()[5]
		x6 = trialPoint.getVariableProxyArr()[6]
		x7 = trialPoint.getVariableProxyArr()[7]
		x8 = trialPoint.getVariableProxyArr()[8]
		x9 = trialPoint.getVariableProxyArr()[9]
		
		x8.setUseInSolver(False)
		x9.setUseInSolver(False)
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
		if optimizerSettingsDictionary.hasKey("floatSecondStripperLength") and optimizerSettingsDictionary.getValue("floatSecondStripperLength")=="True":
			x8.setUseInSolver(True)
		if optimizerSettingsDictionary.hasKey("floatSecondStripperAngle") and optimizerSettingsDictionary.getValue("floatSecondStripperAngle")=="True":
			x9.setUseInSolver(True)			
		solver.solve(scorer,trialPoint)
		
		print "===== best score ========== fitting time = ", solver.getScoreboard().getRunTime()
		
		bestScore = solver.getScoreboard().getBestScore()	
		print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
		
		trialPoint = solver.getScoreboard().getBestTrialPoint()
		
		print trialPoint.textDesciption()
		
		print "(%f,%f,%f,%f)"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3])
		#outputDirectory="WasteBeamClosed"
		fileOut=open("%s/ChicaneScales_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.write("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3],trialPoint.getVariableProxyValuesArr()[4],trialPoint.getVariableProxyValuesArr()[5],trialPoint.getVariableProxyValuesArr()[6],trialPoint.getVariableProxyValuesArr()[7],trialPoint.getVariableProxyValuesArr()[8],trialPoint.getVariableProxyValuesArr()[9]) +"\n")
		fileOut.flush() 
		fileOut.close() 

