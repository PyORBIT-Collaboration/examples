##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################

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
from orbit.teapot import GeneralDipole
from orbit.teapot import YDipole
from orbit.teapot import XDipole
from orbit.teapot import GeneralDipoleNoStrip
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

class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self,OL_teapot_latt_full,OL_teapot_latt_partial,OL_inject_start,OL_inject_end,doDipoleKickers,currentPart=0,nPartsChicane=0,currentPart2=0,nPartsChicane2=0):
		Scorer.__init__(self)
		self.OL_teapot_latt=OL_teapot_latt_full
		self.OL_teapot_latt_full=OL_teapot_latt_full
		self.OL_teapot_latt_partial=OL_teapot_latt_partial
		self.OL_inject_start=OL_inject_start
		self.OL_inject_end=OL_inject_end
		self.currentPart=currentPart
		self.nPartsChicane=nPartsChicane
		self.currentPart2=currentPart2
		self.nPartsChicane2=nPartsChicane2		
		self.doDipoleKickers=doDipoleKickers
		self.addChicaneFieldToStripper=True
		#forDebuggingShouldAlwaysBeTrue
		self.rescaleChicaneFieldInStripper=True
		
		self.xTarget=0
		self.pxTarget=0
		self.yTarget=0
		self.pyTarget=0
		self.zTarget=0
		self.dETarget=0
		
		self.xInit=0
		self.pxInit=0
		self.yInit=0
		self.pyInit=0
		self.zInit=0
		self.dEInit=0	
		
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
		self.b.mass(0.93827231)
		self.b.macroSize(macrosize)
		#energy = 1.0 #Gev
		energy = 1.3 #Gev
		self.b.getSyncParticle().kinEnergy(energy)
		self.b.addParticle(0,0,0,0,0,0)

		#this is the injected bunch
		self.b2 = Bunch()

		sp = self.b.getSyncParticle()
		self.beta= sp.beta()
		self.gamma=sp.gamma()
		c=299792458
		self.rigidity= sp.momentum()/(c/math.pow(10.,9))
		
		self.paramsDict = {}
		self.paramsDict["bunch"]= self.b

		self.paramsDict2 = {}
		self.paramsDict2["bunch"]= self.b2
		
		self.theEffLength=0.03*2
		#theEffLength=0.01
		self.fieldStrength=1.3
		self.fieldStrengthMin=.2
		self.cutLength=0.03
		self.fieldDirection1=-math.pi/2.
		self.fieldDirection2=-math.pi/2.
		#self.fieldDirection1=0
		#self.fieldDirection2=math.pi	
		self.n=1000
		self.maxValue=self.theEffLength
		self.step=self.maxValue/self.n
		
		self.magneticFieldx= Function()
		self.magneticFieldy= Function()	
		self.magneticFieldx2= Function()
		self.magneticFieldy2= Function()	
		
		self.pxOffset=-.042
		self.pyOffset=0

	def getpxOffset(self):
		return self.pxOffset
	def getpyOffset(self):
		return self.pyOffset	
	def setpxOffset(self, pxOffset):
		self.pxOffset=pxOffset
	def setpyOffset(self,pyOffset):
		self.pyOffset=pyOffset			
	def pieceWiseField2(self,x):
		if x<self.cutLength :
			return (self.fieldStrength-self.fieldStrengthMin)/self.cutLength*x+self.fieldStrengthMin
		elif x>=self.cutLength:
			return self.fieldStrength
		pass		
	def resetBunch2(self):
		nParts=10000
		intensity = 7.8e13
		turns=1
		e_kin_ini = 1.3 # in [GeV]
		mass =  0.93827231 #0.939294    # in [GeV]
		gamma = (mass + e_kin_ini)/mass
		beta = math.sqrt(gamma*gamma - 1.0)/gamma
		macrosize = intensity/turns/nParts
		#print "relat. gamma=",gamma
		#print "relat.  beta=",beta
		
		
		#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
		(alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)
		
		(alphaX,betaX,emittX) = (.224, 10.5, 1.445)
		(alphaY,betaY,emittY) = (.224, 10.5, 1.445)
		
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
		bunch_gen = SNS_Linac_BunchGenerator(twissX,twissY,twissZ)
		
		self.b2 = Bunch()

		#generate initial bunch

		self.b2 = bunch_gen.getBunch(nParticles = nParts, distributorClass = WaterBagDist3D)
		#bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
		#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)
	
		xOffset=0.25671
		#pxOffset=-.042
		pxOffset=self.getpxOffset()
		yOffset=0.046
		pyOffset=self.getpyOffset()
		#if reading bunch from file the offset should already have been added
		
		for i in range(self.b2.getSize()):
			self.b2.x(i,self.b2.x(i)+xOffset)
			self.b2.px(i,self.b2.px(i)+pxOffset)
			self.b2.y(i,self.b2.y(i)+yOffset)
			self.b2.py(i,self.b2.py(i)+pyOffset)			
		self.b2.mass(mass) #mass
		self.b2.macroSize(macrosize)
		energy = e_kin_ini # 1.0 #Gev
		self.b2.getSyncParticle().kinEnergy(energy)
		self.b2.charge(-1)
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

	def setScaleChicane(self,chicaneToScale,scale):
		nodes=self.OL_teapot_latt.getTeapotLattice().getNodes()
		for i in range(len(self.OL_teapot_latt.getChicaneNodes()[chicaneToScale])):
			if self.OL_teapot_latt.getChicaneNodes()[chicaneToScale][i] >=0:
				nodes[self.OL_teapot_latt.getChicaneNodes()[chicaneToScale][i]].setParam("kx", self.OL_teapot_latt.getChicaneNodeStrength()[chicaneToScale][i]*scale)
			
		#need to recompute dipole field if stripper is in chicane field because chicane field has been changed
		if chicaneToScale==1 and self.OL_teapot_latt.getFirstDipoleInChicane() and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
			self.magneticFieldx= Function()
			self.magneticFieldy= Function()				
			xkickerField=nodes[self.OL_teapot_latt.getFirstDipoleNode()].getChicaneFieldx()
			ykickerField=nodes[self.OL_teapot_latt.getFirstDipoleNode()].getChicaneFieldy()
			for i in range(self.n):
				x = self.step*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx.add(x,y*math.cos(self.fieldDirection1)+xkickerField*scale)
				self.magneticFieldy.add(x,y*math.sin(self.fieldDirection1)+ykickerField*scale)	
				
			nodes[self.OL_teapot_latt.getFirstDipoleNode()].setFunctionMagneticFieldx(self.magneticFieldx)
			nodes[self.OL_teapot_latt.getFirstDipoleNode()].setFunctionMagneticFieldy(self.magneticFieldy)
			nodes[self.OL_teapot_latt.getFirstDipoleNode()].computeFunctions()
			
		elif chicaneToScale==2 and self.OL_teapot_latt.getSecondDipoleInChicane() and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
			self.magneticFieldx2= Function()
			self.magneticFieldy2= Function()				
			xkickerField=nodes[self.OL_teapot_latt.getSecondDipoleNode()].getChicaneFieldx()
			ykickerField=nodes[self.OL_teapot_latt.getSecondDipoleNode()].getChicaneFieldy()
			for i in range(self.n):
				x = self.step*i;
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
		self.changeLattice(self.OL_inject_end)
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
		if self.doDipoleKickers and self.OL_teapot_latt.getChicaneNodes()[1][0] >=0:
			self.findChicanes()
			#calculate where to place 1st stripper dipole
			position=-100.
			if self.currentPart==-1:
				#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[1][0]]][0]-self.theEffLength
				#print self.OL_teapot_latt.getDriftNodes()
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getDriftNodes()[0][0]]][1]-self.theEffLength
			elif self.currentPart==0:
				print "OL_teapot_latt.getChicaneNodes()[1][0]= ",self.OL_teapot_latt.getChicaneNodes()[1][0]
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[1][0]]][0]
			elif self.currentPart is self.nPartsChicane:
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[1][0]]][1]
			else :
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[1][0]]][0]+self.OL_teapot_latt.getTeapotLattice().getNodes()[OL_teapot_latt.getChicaneNodes()[1][0]].getLength()*self.currentPart/self.nPartsChicane
								
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + self.theEffLength
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
		
			for i in range(self.n):
				x = self.step*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx.add(x,y*math.cos(self.fieldDirection1)+xkickerField)
				self.magneticFieldy.add(x,y*math.sin(self.fieldDirection1)+ykickerField)
			
			if self.OL_teapot_latt.getFirstDipoleIsStripper():
				myDipole_DH_A11=GeneralDipoleStripSeperateField(self.magneticFieldx,self.magneticFieldy,self.n,self.maxValue,self.gamma,self.beta,"Dipole_DH_A11",self.OL_teapot_latt.getFirstDipoleFixedStripLength())
			else:
				myDipole_DH_A11=GeneralDipoleNoStripSeperateField(self.magneticFieldx,self.magneticFieldy,self.n,self.maxValue,self.gamma,self.beta,"Dipole_DH_A11")
			#print "xkickerField=",xkickerField
			myDipole_DH_A11.setChicaneFieldx(xkickerField)
			myDipole_DH_A11.setChicaneFieldy(ykickerField)
			addDipoleStripperNode(self.OL_teapot_latt.getTeapotLattice(),position,myDipole_DH_A11)	
			
	#addSecondStripper
	def addSecondStripperDipole(self):
		if self.doDipoleKickers and self.OL_teapot_latt.getChicaneNodes()[2][0] >=0:
			self.findChicanes()
			#calculate where to place 1st stripper dipole
			position=-100.
			#place second stripper 5/6 of the way into chicane3/12. temporary position for consistency
			#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]+self.teapot_latt.getNodes()[self.chicaneNodes[2][0]].getLength()*5./6.
			if self.currentPart2==-1:
				#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]-self.theEffLength
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getDriftNodes()[1][0]]][1]-self.theEffLength
			elif self.currentPart2==0:
				print "self.chicaneNodes[1][0]= ",self.OL_teapot_latt.getChicaneNodes()[1][0]
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]]][0]
			elif self.currentPart2 is self.nPartsChicane2:
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]]][1]
			else :
				position =self.OL_teapot_latt.getTeapotLattice().getNodePositionsDict()[self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]]][0]+self.OL_teapot_latt.getTeapotLattice().getNodes()[self.OL_teapot_latt.getChicaneNodes()[2][0]].getLength()*self.currentPart2/self.nPartsChicane2			
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + self.theEffLength
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
		
			for i in range(self.n):
				x = self.step*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx2.add(x,y*math.cos(self.fieldDirection2)+xkickerField)
				self.magneticFieldy2.add(x,y*math.sin(self.fieldDirection2)+ykickerField)
			
			if self.OL_teapot_latt.getSecondDipoleIsStripper():
				myDipole_DH_A12=GeneralDipoleStripSeperateField(self.magneticFieldx2,self.magneticFieldy2,self.n,self.maxValue,self.gamma,self.beta,"Dipole_DH_A12",self.OL_teapot_latt.getSecondDipoleFixedStripLength())
			else:
				myDipole_DH_A12=GeneralDipoleNoStripSeperateField(self.magneticFieldx2,self.magneticFieldy2,self.n,self.maxValue,self.gamma,self.beta,"Dipole_DH_A12")
			myDipole_DH_A12.setChicaneFieldx(xkickerField)
			myDipole_DH_A12.setChicaneFieldy(ykickerField)			
			addDipoleStripperNode(self.OL_teapot_latt.getTeapotLattice(),position,myDipole_DH_A12)				
	def getScore(self,trialPoint):
		self.resetBunch()
		#self.resetBunch2()
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()
		x3 = trialPoint.getVariableProxyArr()[3].getValue()
		x4 = trialPoint.getVariableProxyArr()[4].getValue()
		self.setpxOffset(x4)
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
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)
		
		self.changeLattice(self.OL_inject_end)
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)		
		for i in range(self.turns):
			self.OL_teapot_latt_full.getTeapotLattice().trackBunch(self.b, self.paramsDict)
			
		#score = (b.x(0)-self.xTarget)**2 + (b.px(0)-self.pxTarget)**2 + (b.y(0)-self.yTarget)**2+(b.py(0)-self.pyTarget)**2+(b.z(0)-self.zTarget)**2+(b.pz(0)-self.dETarget)**2
		score = (self.b.x(0)-self.xTarget)**2 + (self.b.px(0)-self.pxTarget)**2 + (self.b.y(0)-self.yTarget)**2+(self.b.py(0)-self.pyTarget)**2

		#print "score= ",score, " self.b.x(0)=",self.b.x(0)
		#print "score= ",score, " self.b.px(0)=",self.b.px(0)
		self.resetBunch()
		#print "score= ",score, " self.b.x(0)=",self.b.x(0)
		for i in range(self.turns):
			self.OL_teapot_latt_partial.getTeapotLattice().trackBunch(self.b, self.paramsDict)
		self.OL_inject_start.getTeapotLattice().trackBunch(self.b2, self.paramsDict2)
		#self.OL_inject_end.getTeapotLattice().trackBunch(self.b2, self.paramsDict2)
		twiss_analysis = BunchTwissAnalysis()  
		twiss_analysis.analyzeBunch(self.b2)
		(xavg,xpavg,yavg,ypavg)=(twiss_analysis.getAverage(0),twiss_analysis.getAverage(1),twiss_analysis.getAverage(2),twiss_analysis.getAverage(3))
		score =score +(self.b.px(0)-xpavg)**2 +(self.b.py(0)-ypavg)**2
		#print "score= ",score, " self.b.x(0)=",self.b.x(0)
		#print "self.b.px(0)=",self.b.px(0), " xpavg=",xpavg, " x4=",x4
		return score	
		
print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
parser.add_argument("--addChicaneFieldToStripper",type=bool, dest='addChicaneFieldToStripper', default=True, help="Include the chicane fields in the stripper if stripper is inside chicane")
parser.add_argument("--outputDirectory", dest='outputDirectory', default="ClosureTest", help="Where to put output")
args = parser.parse_args()
doDipoleKickers=args.doDipoleKickers
outputDirectory=args.outputDirectory
addChicaneFieldToStripper=args.addChicaneFieldToStripper
if not os.path.exists(outputDirectory):
	os.mkdir(outputDirectory)
#=====Make a Teapot style lattice======


#the default chicane kick strength array
chicaneStrengthArray=[-0.041456,0.052434,0.0298523,-0.0398609]
chicaneNodes=[29,31,34,36]
#this sets how to divide up chicane2/11 in terms of where 1st stripper is placed.
nPartsChicane=1
nPartsChicane2=0

#for currentPart in range(-1,nPartsChicane+1):
for currentPart in range(1,nPartsChicane+1):
	for currentPart2 in range(-1,nPartsChicane2):
		inj_latt_start = teapot.TEAPOT_Ring()
		print "Read MAD."
		#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
		inj_latt_start.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers_onlyChicane2.LAT","RING")
		#print "Lattice=",inj_latt_start.getName()," length [m] =",inj_latt_start.getLength()," nodes=",len(inj_latt_start.getNodes())
		
		OL_inj_latt_start =OptimizerLattice(inj_latt_start)   
		OL_inj_latt_start.setFirstDipoleIsStripper(True)
		OL_inj_latt_start.setSecondDipoleIsStripper(True)
		OL_inj_latt_start.setSecondDipoleFixedStripLength(.02)
		
		inj_latt_end = teapot.TEAPOT_Ring()
		print "Read MAD."
		#this lattice contains chicane 4 and the drift leading to the waste septum
		inj_latt_end.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers_onlyChicane3.LAT","RING")
		#print "Lattice=",inj_latt_end.getName()," length [m] =",inj_latt_end.getLength()," nodes=",len(inj_latt_end.getNodes())
		
		OL_inj_latt_end =OptimizerLattice(inj_latt_end)
		teapot_latt_full = teapot.TEAPOT_Ring()
		teapot_latt_full.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_KickersJustBeforeQuadAfterChicane4.LAT","RING")
		print "Lattice=",teapot_latt_full.getName()," length [m] =",teapot_latt_full.getLength()," nodes=",len(teapot_latt_full.getNodes())
		
		
		
		teapot_latt_partial = teapot.TEAPOT_Ring()
		teapot_latt_partial.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers_Start_To_JustBeforeChicane4.LAT","RING")
		print "Lattice=",teapot_latt_partial.getName()," length [m] =",teapot_latt_partial.getLength()," nodes=",len(teapot_latt_partial.getNodes())
		
		#Turn off injection kickers
		
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
		hkick10 = nodes[10]
		vkick10 = nodes[12]
		hkick11	= nodes[14]
		vkick11 = nodes[16]
		vkick12 = nodes[49]
		hkick12 = nodes[51]
		vkick13 = nodes[53]
		hkick13	= nodes[55]

		hkick10_2 = nodes2[10]
		vkick10_2 = nodes2[12]
		hkick11_2 = nodes2[14]
		vkick11_2 = nodes2[16]
		
		vkick10.setParam("ky", strength_vkicker10)
		hkick10.setParam("kx", strength_hkicker10)
		vkick11.setParam("ky", strength_vkicker11)
		hkick11.setParam("kx", strength_hkicker11)
		vkick12.setParam("ky", strength_vkicker12)
		hkick12.setParam("kx", strength_hkicker12)
		vkick13.setParam("ky", strength_vkicker13)
		hkick13.setParam("kx", strength_hkicker13)

		vkick10_2.setParam("ky", strength_vkicker10)
		hkick10_2.setParam("kx", strength_hkicker10)
		vkick11_2.setParam("ky", strength_vkicker11)
		hkick11_2.setParam("kx", strength_hkicker11)
		
		vkick10.setWaveform(kickerwave)
		hkick10.setWaveform(kickerwave)
		vkick11.setWaveform(kickerwave)
		hkick11.setWaveform(kickerwave)
		vkick12.setWaveform(kickerwave)
		hkick12.setWaveform(kickerwave)
		vkick13.setWaveform(kickerwave)
		hkick13.setWaveform(kickerwave)

		vkick10_2.setWaveform(kickerwave)
		hkick10_2.setWaveform(kickerwave)
		vkick11_2.setWaveform(kickerwave)
		hkick11_2.setWaveform(kickerwave)
		
		teapot_latt_full.initialize()
		teapot_latt_partial.initialize()
		
		
		OL_teapot_latt_full =OptimizerLattice(teapot_latt_full)
		OL_teapot_latt_partial =OptimizerLattice(teapot_latt_partial)
		
		
		
		#scorer = MyScorer(0.005769,0.002069,0.001778,-0.000359,-0.003845,0.000000)
		#scorer = MyScorer(0.004334,0.000192,0.001710,-0.000349,-0.004286,0.000000)
		#scorer = MyScorer(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000)
		
		scorer = MyScorer(OL_teapot_latt_full,OL_teapot_latt_partial,OL_inj_latt_start,OL_inj_latt_end,doDipoleKickers,currentPart,nPartsChicane,currentPart2,nPartsChicane2)
		scorer.initScorer()	
		
		#searchAlgorithm   = RandomSearchAlgorithm()
		searchAlgorithm = SimplexSearchAlgorithm()
		
		#max_time = 0.05
		max_time = 100
		max_accuracy=1E-10
		max_iterations=1000
		
		iterationStopper=SolveStopperFactory.maxIterationStopper(max_iterations)
		timeStopper = SolveStopperFactory.maxTimeStopper(max_time)
		solverStopper=SolveStopperFactory.comboStopper()
		solverStopper.addStopper(iterationStopper)
		solverStopper.addStopper(timeStopper)
		#solverStopper = SolveStopperFactory.maxAccuracyStopper(max_accuracy)
		
		solver = Solver()
		solver.setAlgorithm(searchAlgorithm)
		solver.setStopper(solverStopper)
		
		trialPoint = TrialPoint()
		#trialPoint.addVariableProxy(VariableProxy(name = "x0", value = 1., step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x0", value = 1., step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x1", value = 1., step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x2", value = 1., step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x3", value = 1., step = 0.1))
		trialPoint.addVariableProxy(VariableProxy(name = "x4", value = -.042, step = 0.01))
		x0 = trialPoint.getVariableProxyArr()[0]
		x1 = trialPoint.getVariableProxyArr()[1]
		x2 = trialPoint.getVariableProxyArr()[2]
		x3 = trialPoint.getVariableProxyArr()[3]
		x4 = trialPoint.getVariableProxyArr()[4]
		
		x4.setUseInSolver(False)
		solver.solve(scorer,trialPoint)
		
		print "===== best score ========== fitting time = ", solver.getScoreboard().getRunTime()
		
		bestScore = solver.getScoreboard().getBestScore()	
		print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
		
		trialPoint = solver.getScoreboard().getBestTrialPoint()
		
		print trialPoint.textDesciption()
		
		print "(%f,%f,%f,%f)"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3])
		#outputDirectory="WasteBeamClosed"
		fileOut=open("%s/ChicaneScales_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.write("%f,%f,%f,%f,%f"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3],trialPoint.getVariableProxyValuesArr()[4]) +"\n")
		fileOut.flush() 
		fileOut.close() 
		
		#scorer.getScore(trialPoint)

