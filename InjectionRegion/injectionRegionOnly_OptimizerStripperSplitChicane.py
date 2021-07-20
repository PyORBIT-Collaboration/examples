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
from orbit.teapot import addDipoleStripperNode

from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self,teapot_latt,doDipoleKickers,currentPart=0,nPartsChicane=0,currentPart2=0,nPartsChicane2=0):
		Scorer.__init__(self)
		self.teapot_latt=teapot_latt
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
		
		self.b = Bunch()
		self.b.mass(0.93827231)
		self.b.macroSize(macrosize)
		#energy = 1.0 #Gev
		energy = 1.3 #Gev
		self.b.getSyncParticle().kinEnergy(energy)
		self.b.addParticle(0,0,0,0,0,0)

		sp = self.b.getSyncParticle()
		self.beta= sp.beta()
		self.gamma=sp.gamma()
		c=299792458
		self.rigidity= sp.momentum()/(c/math.pow(10.,9))
		
		self.paramsDict = {}
		self.paramsDict["bunch"]= self.b
		
		self.theEffLength=0.03*2
		#theEffLength=0.01
		self.fieldStrength=1.3
		self.fieldStrengthMin=.2
		self.cutLength=0.03
		self.fieldDirection1=math.pi/2.
		self.fieldDirection2=-math.pi/2.
		#self.fieldDirection1=0
		#self.fieldDirection2=math.pi	
		self.n=1000
		self.maxValue=self.theEffLength
		self.step=self.maxValue/self.n
		
		self.firstDipoleInChicane=False
		self.secondDipoleInChicane=False
		self.firstDipoleNode=-1
		self.secondDipoleNode=-1
		
		self.magneticFieldx= Function()
		self.magneticFieldy= Function()	
		self.magneticFieldx2= Function()
		self.magneticFieldy2= Function()	
		#[0-3] are chicane10-13, 
		self.chicaneNodes=[]
		#[0] is drift DB12 and [1] is drift DB23
		self.driftNodes=[]

		self.chicaneNodeStrength=[[0.041456],[-0.052434],[-0.0298523],[0.0398609]]
	def pieceWiseField2(self,x):
		if x<self.cutLength :
			return (self.fieldStrength-self.fieldStrengthMin)/self.cutLength*x+self.fieldStrengthMin
		elif x>=self.cutLength:
			return self.fieldStrength
		pass		
	#currently assumes bunch has just one macroparticle
	def resetBunch(self):
		self.b.x(0,self.xInit)
		self.b.px(0,self.pxInit)
		self.b.y(0,self.yInit)
		self.b.py(0,self.pyInit)
		self.b.z(0,self.zInit)
		self.b.dE(0,self.dEInit)

	def setScaleChicane(self,chicaneToScale,scale):
		nodes=self.teapot_latt.getNodes()
		for i in range(len(self.chicaneNodes[chicaneToScale])):
			nodes[self.chicaneNodes[chicaneToScale][i]].setParam("kx", self.chicaneNodeStrength[chicaneToScale][i]*scale)
			
		#need to recompute dipole field if stripper is in chicane field because chicane field has been changed
		if chicaneToScale==1 and self.firstDipoleInChicane and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
			self.magneticFieldx= Function()
			self.magneticFieldy= Function()				
			xkickerField=nodes[self.firstDipoleNode].getChicaneFieldx()
			ykickerField=nodes[self.firstDipoleNode].getChicaneFieldy()
			for i in range(self.n):
				x = self.step*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx.add(x,y*math.cos(self.fieldDirection1)+xkickerField*scale)
				self.magneticFieldy.add(x,y*math.sin(self.fieldDirection1)+ykickerField*scale)	
				
			nodes[self.firstDipoleNode].setFunctionMagneticFieldx(self.magneticFieldx)
			nodes[self.firstDipoleNode].setFunctionMagneticFieldy(self.magneticFieldy)
			nodes[self.firstDipoleNode].computeFunctions()
			
		elif chicaneToScale==2 and self.secondDipoleInChicane and self.addChicaneFieldToStripper and self.rescaleChicaneFieldInStripper:
			self.magneticFieldx2= Function()
			self.magneticFieldy2= Function()				
			xkickerField=nodes[self.secondDipoleNode].getChicaneFieldx()
			ykickerField=nodes[self.secondDipoleNode].getChicaneFieldy()
			for i in range(self.n):
				x = self.step*i;
				#y = constantField(x)
				y = self.pieceWiseField2(x)
				self.magneticFieldx2.add(x,y*math.cos(self.fieldDirection2)+xkickerField*scale)
				self.magneticFieldy2.add(x,y*math.sin(self.fieldDirection2)+ykickerField*scale)	
				
			nodes[self.secondDipoleNode].setFunctionMagneticFieldx(self.magneticFieldx2)
			nodes[self.secondDipoleNode].setFunctionMagneticFieldy(self.magneticFieldy2)
			nodes[self.secondDipoleNode].computeFunctions()			
	def getScaleChicane(self,chicaneToScale):
		nodeIndex=self.chicaneNodes[chicaneToScale]
		chicaneNode=self.teapot_latt.getNodes()[nodeIndex]
		return chicaneNode.getParam("kx")
	#used to change the location of the chicane node if different from original location
	def changeChicaneNode(self, nodeToChange, newNodeLocation):
		self.chicaneNodes[nodeToChange]=newNodeLocation
	#changes the starting strength of the chicane
	def setStartStrengthChicane(self,chicaneToChange,strength):
		self.chicaneNodeStrength[chicaneToChange]=strength
	#gets the starting strength of the chicane
	def getStartStrengthChicane(self,chicaneToChange):
		return self.chicaneNodeStrength[chicaneToChange]
	#find location of stripper nodes
	def findStrippers(self):
		index=0
		for node in teapot_latt.getNodes():
			if node.getName().strip() == "Dipole_DH_A11":
				self.firstDipoleNode=index
			elif node.getName().strip() == "Dipole_DH_A12":
				self.secondDipoleNode=index
			index+=1		
	#find location of chicane nodes
	def findChicanes(self):
		#[0-3] are chicane10-13,
		self.chicaneNodes=[]
		# [0] is drift DB12 and [1] is drift DB23
		self.driftNodes=[]
		chicane10Nodes=[]
		chicane11Nodes=[]
		chicane12Nodes=[]
		chicane13Nodes=[]
		db12Nodes=[]
		db23Nodes=[]
		index=0
		for node in teapot_latt.getNodes():
			if node.getName().strip() == "DH_A10":
				chicane10Nodes.append(index)
			elif node.getName().strip() == "DH_A11":
				chicane11Nodes.append(index)
			elif node.getName().strip() == "DH_A12":
				chicane12Nodes.append(index)
			elif node.getName().strip() == "DH_A13":
				chicane13Nodes.append(index)
			elif node.getName().strip() == "DB12":
				db12Nodes.append(index)
			elif node.getName().strip() == "DB23":
				db23Nodes.append(index)				
			index+=1
		self.chicaneNodes.append(chicane10Nodes)
		self.chicaneNodes.append(chicane11Nodes)
		self.chicaneNodes.append(chicane12Nodes)
		self.chicaneNodes.append(chicane13Nodes)
		self.driftNodes.append(db12Nodes)
		self.driftNodes.append(db23Nodes)
	def findChicaneStrength(self):	
		self.findChicanes()
		self.chicaneNodeStrength=[]
		#chicane10Str=[]
		#chicane11Str=[]
		#chicane12Str=[]
		#chicane13Str=[]
		for chicane in self.chicaneNodes:
			chicaneStr=[]
			for i in chicane:
				chicaneStr.append(self.teapot_latt.getNodes()[i].getParam("kx"))
			self.chicaneNodeStrength.append(chicaneStr)
	#initialize the chicanes in the teapot lattice.
	def initChicanes(self):
		self.findChicanes()
		chicanewave = flatTopWaveform(1.0)
		nodes = self.teapot_latt.getNodes()
		for chicane in range(len(self.chicaneNodes)):
			for i in range(len(self.chicaneNodes[chicane])):
				nodes[self.chicaneNodes[chicane][i]].setWaveform(chicanewave)
				print chicane, "i= ",i
				print self.chicaneNodeStrength[chicane][i]
				nodes[self.chicaneNodes[chicane][i]].setParam("kx",self.chicaneNodeStrength[chicane][i])	
		self.addFirstStripperDipole()
		self.addSecondStripperDipole()
		self.findChicanes()
		self.findChicaneStrength()
		self.findStrippers()				
	#change the lattice
	def changeLattice(self,lattice):
		self.teapot_latt=lattice
	#addFirstStripper
	def addFirstStripperDipole(self):
		if self.doDipoleKickers:
			self.findChicanes()
			#calculate where to place 1st stripper dipole
			position=-100.
			if self.currentPart==-1:
				#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[1][0]]][0]-self.theEffLength
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.driftNodes[0][0]]][1]-self.theEffLength
			elif self.currentPart==0:
				print "self.chicaneNodes[1][0]= ",self.chicaneNodes[1][0]
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[1][0]]][0]
			elif self.currentPart is self.nPartsChicane:
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[1][0]]][1]
			else :
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[1][0]]][0]+self.teapot_latt.getNodes()[self.chicaneNodes[1][0]].getLength()*self.currentPart/self.nPartsChicane
								
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + self.theEffLength
			(node_start_ind,node_stop_ind,z,ind) = (-1,-1, 0., 0)
			for nodeCurrent in self.teapot_latt.getNodes():
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
			nodeCurrent=self.teapot_latt.getNodes()[node_start_ind]
			xkickerField=0.
			ykickerField=0.
			#nothing to change because no field in drift
			if(isinstance(nodeCurrent,DriftTEAPOT)):
				self.firstDipoleInChicane=False
				print "stripper dipole is in drift"
			#in kicker so add kicker field to stripping field
			elif (isinstance(nodeCurrent,KickTEAPOT)):
				print "stripper dipole is in kick node"
				if self.addChicaneFieldToStripper:
					self.firstDipoleInChicane=True
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
			
			myDipole_DH_A11=GeneralDipoleNoStripSeperateField(self.magneticFieldx,self.magneticFieldy,self.n,self.maxValue,self.gamma,self.beta,"Dipole_DH_A11")
			#print "xkickerField=",xkickerField
			myDipole_DH_A11.setChicaneFieldx(xkickerField)
			myDipole_DH_A11.setChicaneFieldy(ykickerField)
			addDipoleStripperNode(self.teapot_latt,position,myDipole_DH_A11)	
			
	#addSecondStripper
	def addSecondStripperDipole(self):
		if self.doDipoleKickers:
			self.findChicanes()
			#calculate where to place 1st stripper dipole
			position=-100.
			#place second stripper 5/6 of the way into chicane3/12. temporary position for consistency
			#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]+self.teapot_latt.getNodes()[self.chicaneNodes[2][0]].getLength()*5./6.
			if self.currentPart2==-1:
				#position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]-self.theEffLength
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.driftNodes[1][0]]][1]-self.theEffLength
			elif self.currentPart2==0:
				print "self.chicaneNodes[1][0]= ",self.chicaneNodes[1][0]
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]
			elif self.currentPart2 is self.nPartsChicane2:
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][1]
			else :
				position =self.teapot_latt.getNodePositionsDict()[self.teapot_latt.getNodes()[self.chicaneNodes[2][0]]][0]+self.teapot_latt.getNodes()[self.chicaneNodes[2][0]].getLength()*self.currentPart2/self.nPartsChicane2			
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + self.theEffLength
			(node_start_ind,node_stop_ind,z,ind) = (-1,-1, 0., 0)
			for nodeCurrent in self.teapot_latt.getNodes():
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
			nodeCurrent=self.teapot_latt.getNodes()[node_start_ind]
			xkickerField=0.
			ykickerField=0.
			#nothing to change because no field in drift
			if(isinstance(nodeCurrent,DriftTEAPOT)):
				self.secondDipoleInChicane=False
				print "stripper dipole is in drift"
			#in kicker so add kicker field to stripping field
			elif (isinstance(nodeCurrent,KickTEAPOT)):
				print "stripper dipole is in kick node"
				if self.addChicaneFieldToStripper:
					self.secondDipoleInChicane=True
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
			
			myDipole_DH_A12=GeneralDipoleNoStripSeperateField(self.magneticFieldx2,self.magneticFieldy2,self.n,self.maxValue,self.gamma,self.beta,"Dipole_DH_A12")
			myDipole_DH_A12.setChicaneFieldx(xkickerField)
			myDipole_DH_A12.setChicaneFieldy(ykickerField)			
			addDipoleStripperNode(self.teapot_latt,position,myDipole_DH_A12)				
	def getScore(self,trialPoint):
		self.resetBunch()
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()
		x3 = trialPoint.getVariableProxyArr()[3].getValue()
		self.setScaleChicane(0,x0)
		self.setScaleChicane(1,x1)
		self.setScaleChicane(2,x2)
		self.setScaleChicane(3,x3)

		for i in range(self.turns):
			self.teapot_latt.trackBunch(self.b, self.paramsDict)
			
		#score = (b.x(0)-self.xTarget)**2 + (b.px(0)-self.pxTarget)**2 + (b.y(0)-self.yTarget)**2+(b.py(0)-self.pyTarget)**2+(b.z(0)-self.zTarget)**2+(b.pz(0)-self.dETarget)**2
		score = (self.b.x(0)-self.xTarget)**2 + (self.b.px(0)-self.pxTarget)**2 + (self.b.y(0)-self.yTarget)**2+(self.b.py(0)-self.pyTarget)**2

		return score	
		
print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
parser.add_argument("--addChicaneFieldToStripper",type=bool, dest='addChicaneFieldToStripper', default=True, help="Include the chicane fields in the stripper if stripper is inside chicane")
parser.add_argument("--outputDirectory", dest='outputDirectory', default="InjectBeam4_NoDipoles", help="Where to put output")
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

for currentPart in range(1,nPartsChicane+1):
	for currentPart2 in range(-1,nPartsChicane2):
		teapot_latt = teapot.TEAPOT_Ring()
		teapot_latt.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_KickersJustBeforeQuadAfterChicane4.LAT","RING")
		#print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())
		nodes2 = teapot_latt.getNodes()
			
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
		
		nodes = teapot_latt.getNodes()
		hkick10 = nodes[10]
		vkick10 = nodes[12]
		hkick11	= nodes[14]
		vkick11 = nodes[16]
		vkick12 = nodes[49]
		hkick12 = nodes[51]
		vkick13 = nodes[53]
		hkick13	= nodes[55]
		
		vkick10.setParam("ky", strength_vkicker10)
		hkick10.setParam("kx", strength_hkicker10)
		vkick11.setParam("ky", strength_vkicker11)
		hkick11.setParam("kx", strength_hkicker11)
		vkick12.setParam("ky", strength_vkicker12)
		hkick12.setParam("kx", strength_hkicker12)
		vkick13.setParam("ky", strength_vkicker13)
		hkick13.setParam("kx", strength_hkicker13)
		
		vkick10.setWaveform(kickerwave)
		hkick10.setWaveform(kickerwave)
		vkick11.setWaveform(kickerwave)
		hkick11.setWaveform(kickerwave)
		vkick12.setWaveform(kickerwave)
		hkick12.setWaveform(kickerwave)
		vkick13.setWaveform(kickerwave)
		hkick13.setWaveform(kickerwave)
			
		teapot_latt.initialize()
		
		
		
		
		
		
		
		#scorer = MyScorer(0.005769,0.002069,0.001778,-0.000359,-0.003845,0.000000)
		#scorer = MyScorer(0.004334,0.000192,0.001710,-0.000349,-0.004286,0.000000)
		#scorer = MyScorer(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000)
		
		scorer = MyScorer(teapot_latt,doDipoleKickers,currentPart,nPartsChicane,currentPart2,nPartsChicane2)
		scorer.initChicanes()
		
		#searchAlgorithm   = RandomSearchAlgorithm()
		searchAlgorithm = SimplexSearchAlgorithm()
		
		#max_time = 0.05
		max_time = 10
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
		x0 = trialPoint.getVariableProxyArr()[0]
		x1 = trialPoint.getVariableProxyArr()[1]
		x2 = trialPoint.getVariableProxyArr()[2]
		x3 = trialPoint.getVariableProxyArr()[3]
		
		#x1.setUseInSolver(False)
		solver.solve(scorer,trialPoint)
		
		print "===== best score ========== fitting time = ", solver.getScoreboard().getRunTime()
		
		bestScore = solver.getScoreboard().getBestScore()	
		print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
		
		trialPoint = solver.getScoreboard().getBestTrialPoint()
		
		print trialPoint.textDesciption()
		
		print "(%f,%f,%f,%f)"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3])
		#outputDirectory="WasteBeamClosed"
		fileOut=open("%s/ChicaneScales_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.write("%f,%f,%f,%f"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3]) +"\n")
		fileOut.flush() 
		fileOut.close() 

