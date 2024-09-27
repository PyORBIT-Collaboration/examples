##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
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
from KevinPython.function_stripping import probabilityStripping

from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self,teapot_latt):
		Scorer.__init__(self)
		self.teapot_latt=teapot_latt
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
		
		self.paramsDict = {}
		self.paramsDict["bunch"]= self.b
		self.chicaneNodes=[29,31,34,36]
		self.chicaneNodeStrength=[-0.041456,0.052434,0.0298523,-0.0398609]
	#currently assumes bunch has just one macroparticle
	def resetBunch(self):
		self.b.x(0,self.xInit)
		self.b.px(0,self.pxInit)
		self.b.y(0,self.yInit)
		self.b.py(0,self.pyInit)
		self.b.z(0,self.zInit)
		self.b.dE(0,self.dEInit)

	def setScaleChicane(self,chicaneToScale,scale):
		nodeIndex=self.chicaneNodes[chicaneToScale]
		chicaneNode=self.teapot_latt.getNodes()[nodeIndex]
		chicaneNode.setParam("kx", self.chicaneNodeStrength[chicaneToScale]*scale)
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
	#initialize the chicanes in the teapot lattice.
	def initChicanes(self):
		chicanewave = flatTopWaveform(1.0)
		nodes = self.teapot_latt.getNodes()
		for i in range(4):
			nodes[self.chicaneNodes[i]].setWaveform(chicanewave)
			nodes[self.chicaneNodes[i]].setParam("kx",self.chicaneNodeStrength[i])
		pass		
	def getScore(self,trialPoint):
		self.resetBunch()
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()
		x3 = trialPoint.getVariableProxyArr()[3].getValue()
		self.setScaleChicane(0,-x0)
		self.setScaleChicane(1,-x1)
		self.setScaleChicane(2,-x2)
		self.setScaleChicane(3,-x3)

		for i in range(self.turns):
			self.teapot_latt.trackBunch(self.b, self.paramsDict)
			
		#score = (b.x(0)-self.xTarget)**2 + (b.px(0)-self.pxTarget)**2 + (b.y(0)-self.yTarget)**2+(b.py(0)-self.pyTarget)**2+(b.z(0)-self.zTarget)**2+(b.pz(0)-self.dETarget)**2
		score = (self.b.x(0)-self.xTarget)**2 + (self.b.px(0)-self.pxTarget)**2 + (self.b.y(0)-self.yTarget)**2+(self.b.py(0)-self.pyTarget)**2

		return score	
		
print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
args = parser.parse_args()
#=====Make a Teapot style lattice======
theEffLength=0.03*2
#theEffLength=0.01
fieldStrength=1.3
fieldStrengthMin=.2
cutLength=0.03
A1=2.47e-6
A2=4.49e9

#=====Main bunch parameters============
intensity = 7.8e13
#turns = 1000.0
turns=1
#macrosperturn = 260
macrosperturn = 1
macrosize = intensity/turns/macrosperturn
bunch_in = Bunch()
mass =  0.93827231 #0.939294    # in [GeV]
bunch_in.mass(mass) #mass
bunch_in.macroSize(macrosize)
e_kin_ini = 1.3
energy = e_kin_ini # 1.0 #Gev
bunch_in.getSyncParticle().kinEnergy(energy)
bunch_in.charge(+1)

sp = bunch_in.getSyncParticle()
beta= sp.beta()
gamma=sp.gamma()

print "beta= %f"%beta
print "gamma= %f"%gamma

def constantField(x):
	return fieldStrength
def pieceWiseField(x):
	if x<=0.0:
		return 0
	elif x<cutLength :
		return fieldStrength/cutLength*x
	elif x>=cutLength:
		return fieldStrength
	pass
def pieceWiseField2(x):
	if x<cutLength :
		return (fieldStrength-fieldStrengthMin)/cutLength*x+fieldStrengthMin
	elif x>=cutLength:
		return fieldStrength
	pass
magneticField= Function()
n=1000
maxValue=theEffLength
step=maxValue/n

for i in range(n):
	x = step*i;
	#y = constantField(x)
	y = pieceWiseField2(x)
	magneticField.add(x,y)
	
theStrippingFunctions=probabilityStripping(magneticField,n,maxValue,gamma,beta)
theStrippingFunctions.computeFunctions()
accumlatedSum=theStrippingFunctions.getaccumlatedSum()
CDF=theStrippingFunctions.getCDF()
deltaxp_rigidity=theStrippingFunctions.getdeltaxp_rigidity()
deltax_rigidity=theStrippingFunctions.getdeltax_rigidity()
deltaxp_m_rigidity=theStrippingFunctions.getdeltaxp_m_rigidity()
deltax_m_rigidity=theStrippingFunctions.getdeltax_m_rigidity()
InverseFunction=theStrippingFunctions.getInverseFunction()

nPartsChicane=4
for currentPart in range(nPartsChicane+1):
	teapot_latt = teapot.TEAPOT_Ring()
	teapot_latt.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers.LAT","RING")
	#print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())
	nodes2 = teapot_latt.getNodes()
	#numberOfCustomDipoles=2
	numberOfParts_DH_A11=nPartsChicane
	numberOfParts_DH_A12=4
	counterDH_A11=0
	counterDH_A12=0
	if args.doDipoleKickers:
		for node in nodes2:
			if (node.getName().strip()=="DH_A11"):
				#pass
				node.setnParts(numberOfParts_DH_A11)
				print "total length= ",node.getLength()
				#print "segment length= ",node.getLength(3)
				#myDipole_DH_A11=YDipole("Dipole_DH_A11")
				myDipole_DH_A11=GeneralDipoleNoStrip("Dipole_DH_A11")
				myDipole_DH_A11.setFunctionXPRigidity(deltaxp_rigidity)
				myDipole_DH_A11.setFunctionXRigidity(deltax_rigidity)
				myDipole_DH_A11.setFieldDirection(math.pi/2)
				myDipole_DH_A11.setEffLength(theEffLength)
				if currentPart is not nPartsChicane:
					node.addChildNode(myDipole_DH_A11,AccNode.BODY,currentPart)
				else:
					pass
					#node.addChildNode(myEmitNode_DH11_3pre,AccNode.EXIT,currentPart-1)
					#node.addChildNode(myPrintNode_DH11_3pre,AccNode.EXIT,currentPart-1)
					#node.addChildNode(myDipole_DH_A11,AccNode.EXIT,currentPart-1)
			if (node.getName().strip()=="DB23"):
				node.setnParts(2)
				if currentPart is nPartsChicane:
					myDipole_DH_A11=GeneralDipoleNoStrip("Dipole_DH_A11")
					myDipole_DH_A11.setFunctionXPRigidity(deltaxp_rigidity)
					myDipole_DH_A11.setFunctionXRigidity(deltax_rigidity)
					myDipole_DH_A11.setFieldDirection(math.pi/2)
					myDipole_DH_A11.setEffLength(theEffLength)					
					node.addChildNode(myDipole_DH_A11,AccNode.BODY,0)				
			if (node.getName().strip()=="DH_A12"):
				node.setnParts(numberOfParts_DH_A12)
				print "total length= ",node.getLength()
				print "segment length= ",node.getLength(3)
				#myDipole_DH_A12=YDipole("Dipole_DH_A12")
				myDipole_DH_A12=GeneralDipoleNoStrip("Dipole_DH_A12")
				myDipole_DH_A12.setFunctionXPRigidity(deltaxp_rigidity)
				myDipole_DH_A12.setFunctionXRigidity(deltax_rigidity)
				myDipole_DH_A12.setFieldDirection(math.pi/2)
				myDipole_DH_A12.setEffLength(theEffLength)			
				node.addChildNode(myDipole_DH_A12,AccNode.BODY,3)
		
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
	
	scorer = MyScorer(teapot_latt)
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
	
	#x0.setUseInSolver(False)
	solver.solve(scorer,trialPoint)
	
	print "===== best score ========== fitting time = ", solver.getScoreboard().getRunTime()
	
	bestScore = solver.getScoreboard().getBestScore()	
	print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
	
	trialPoint = solver.getScoreboard().getBestTrialPoint()
	
	print trialPoint.textDesciption()
	
	print "(%f,%f,%f,%f)"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3])
	outputDirectory="WasteBeamClosed"
	fileOut=open("%s/ChicaneScales_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.write("%f,%f,%f,%f"%(trialPoint.getVariableProxyValuesArr()[0],trialPoint.getVariableProxyValuesArr()[1],trialPoint.getVariableProxyValuesArr()[2],trialPoint.getVariableProxyValuesArr()[3]) +"\n")
	fileOut.flush() 
	fileOut.close() 

