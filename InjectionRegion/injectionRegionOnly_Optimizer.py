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
from KevinPython.printNode import Print_Node
import argparse

from orbit.utils.fitting import Solver
from orbit.utils.fitting import Scorer
from orbit.utils.fitting import SolveStopperFactory
from orbit.utils.fitting import VariableProxy
from orbit.utils.fitting import TrialPoint

from orbit.utils.fitting import SimplexSearchAlgorithm


class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self,xTarget,pxTarget,yTarget,pyTarget,zTarget,dETarget):
		Scorer.__init__(self)
		self.xTarget=xTarget
		self.pxTarget=pxTarget
		self.yTarget=yTarget
		self.pyTarget=pyTarget
		self.zTarget=zTarget
		self.dETarget=dETarget

	def getScore(self,trialPoint):
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()
		x3 = trialPoint.getVariableProxyArr()[3].getValue()
		print "Start."
		parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
		parser.add_argument("--fileName", dest='fileName', default="outputAddMagnetInjectionRegionAllNodesOG.txt", help="file to print node info into")
		parser.add_argument("--nParts",type=int, dest='nParts', default=1, help="number of particles")
		parser.add_argument("--turns",type=int, dest='turns', default=1, help="number of complete orbits")
		#parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=66, help="What node to monitor")
		parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=-1, help="What node to monitor")
		parser.add_argument("--printNodes",type=bool, dest='printNodes', default=True, help="print node list")
		parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=False, help="print node list")
		parser.add_argument("--x",type=float, dest='x', default=0, help="x position")
		parser.add_argument("--px",type=float, dest='px', default=0, help="xp rad")
		parser.add_argument("--y",type=float, dest='y', default=0, help="y position")
		parser.add_argument("--py",type=float, dest='py', default=0, help="yp rad")
		parser.add_argument("--z",type=float, dest='z', default=0, help="z position")
		parser.add_argument("--pz",type=float, dest='pz', default=0, help="dE position")
		parser.add_argument("--doNormalKickers",type=bool, dest='doNormalKickers', default=False, help="print node list")
		parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=1, help="scaleChicane10")
		parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=1, help="scaleChicane11")
		parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=1, help="scaleChicane12")
		parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=1, help="scaleChicane13")
		#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=0.8739071, help="scaleChicane12")
		#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=0.9329522, help="scaleChicane13")
		args = parser.parse_args()
		scaleChicane10=args.scaleChicane10
		scaleChicane11=args.scaleChicane11
		scaleChicane12=args.scaleChicane12
		scaleChicane13=args.scaleChicane13
		#=====Main bunch parameters============
		intensity = 7.8e13
		#turns = 1000.0
		turns = args.turns
		#macrosperturn = 260
		macrosperturn = args.nParts
		macrosize = intensity/turns/macrosperturn
		
		b = Bunch()
		b.mass(0.93827231)
		b.macroSize(macrosize)
		energy = 1.0 #Gev
		b.getSyncParticle().kinEnergy(energy)
		
		paramsDict = {}
		lostbunch = Bunch()
		paramsDict["lostbunch"]=lostbunch
		paramsDict["bunch"]= b
		lostbunch.addPartAttr("LostParticleAttributes") 
		
		#=====Make a Teapot style lattice======
		
		teapot_latt = teapot.TEAPOT_Ring()
		#print "Read MAD."
		teapot_latt.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers.LAT","RING")
		#print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())
		nodes2 = teapot_latt.getNodes()
		#numberOfCustomDipoles=2
		numberOfParts_DH_A11=4
		numberOfParts_DH_A12=4
		counterDH_A11=0
		counterDH_A12=0
		if args.doDipoleKickers:
			for node in nodes2:
				if (node.getName().strip()=="DH_A11"):
					#pass
					node.setnParts(numberOfParts_DH_A11)
					#print "total length= ",node.getLength()
					#print "segment length= ",node.getLength(3)
					#myDipole_DH_A11=YDipole("Dipole_DH_A11")
					myDipole_DH_A11=GeneralDipole("Dipole_DH_A11")
					myDipole_DH_A11.setMagneticFieldStrength(1.0)
					myDipole_DH_A11.setFieldDirection(math.pi/2)
					myDipole_DH_A11.setEffLength(.0254)
					node.addChildNode(myDipole_DH_A11,AccNode.BODY,3)
				if (node.getName().strip()=="DH_A12"):
					node.setnParts(numberOfParts_DH_A12)
					#print "total length= ",node.getLength()
					#print "segment length= ",node.getLength(3)
					#myDipole_DH_A12=YDipole("Dipole_DH_A12")
					myDipole_DH_A12=GeneralDipole("Dipole_DH_A12")
					myDipole_DH_A12.setMagneticFieldStrength(-1.0)
					myDipole_DH_A12.setFieldDirection(math.pi/2)
					myDipole_DH_A12.setEffLength(.0254)
					node.addChildNode(myDipole_DH_A12,AccNode.BODY,3)	
				
		#====Add the injection kickers======
		
		duration = 0.001
		startamp = 1.0
		endamp = 0.58
		deltapos = 0.001
		
		if args.doNormalKickers:
			strength_hkicker10 = 7.211536E-03
			strength_hkicker13 = strength_hkicker10
			strength_hkicker11 = -2.278306E-03
			strength_hkicker12 = strength_hkicker11
			strength_vkicker10 = 4.188402E-03
			strength_vkicker13 = strength_vkicker10
			strength_vkicker11 = -2.118213E-03
			strength_vkicker12 = strength_vkicker11
		else:
			strength_hkicker10 = 0
			strength_hkicker13 = strength_hkicker10
			strength_hkicker11 = 0
			strength_hkicker12 = strength_hkicker11
			strength_vkicker10 = 0
			strength_vkicker13 = strength_vkicker10
			strength_vkicker11 = 0
			strength_vkicker12 = strength_vkicker11	
		
		
		strength_chicane10 = -0.041456 *scaleChicane10
		strength_chicane11 = 0.052434*scaleChicane11
		strength_chicane12 = 0.0298523*scaleChicane12
		strength_chicane13 = -0.0398609*scaleChicane13
		
		lattlength = teapot_latt.getLength()
		sp = b.getSyncParticle()
		kickerwave = rootTWaveform(sp, lattlength, duration, startamp, endamp)
		chicanewave = flatTopWaveform(1.0)
		
		nodes = teapot_latt.getNodes()
		hkick10 = nodes[10]
		vkick10 = nodes[12]
		hkick11	= nodes[14]
		vkick11 = nodes[16]
		vkick12 = nodes[49]
		hkick12 = nodes[51]
		vkick13 = nodes[53]
		hkick13	= nodes[55]
		
		chicane10 = nodes[29]
		chicane11 = nodes[31]
		chicane12 = nodes[34]
		chicane13 = nodes[36]
		
		vkick10.setParam("ky", strength_vkicker10)
		hkick10.setParam("kx", strength_hkicker10)
		vkick11.setParam("ky", strength_vkicker11)
		hkick11.setParam("kx", strength_hkicker11)
		vkick12.setParam("ky", strength_vkicker12)
		hkick12.setParam("kx", strength_hkicker12)
		vkick13.setParam("ky", strength_vkicker13)
		hkick13.setParam("kx", strength_hkicker13)
		
		chicane10.setParam("kx", strength_chicane10*x0)
		chicane11.setParam("kx", strength_chicane11*x1)
		chicane12.setParam("kx", strength_chicane12*x2)
		chicane13.setParam("kx", strength_chicane13*x3)
		
		vkick10.setWaveform(kickerwave)
		hkick10.setWaveform(kickerwave)
		vkick11.setWaveform(kickerwave)
		hkick11.setWaveform(kickerwave)
		vkick12.setWaveform(kickerwave)
		hkick12.setWaveform(kickerwave)
		vkick13.setWaveform(kickerwave)
		hkick13.setWaveform(kickerwave)
		
		chicane10.setWaveform(chicanewave)
		chicane11.setWaveform(chicanewave)
		chicane12.setWaveform(chicanewave)
		chicane13.setWaveform(chicanewave)
		#print "chicane10= ",chicane10.getParam("kx")
		
		
		#for node in kickernode:
		#print "node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
		#print "There are ", node.getNumberOfBodyChildren()," child nodes."
		
		#------------------------------
		#Initial Distribution Functions
		#------------------------------
		
		sp = b.getSyncParticle()
		
		order = 3.
		alphax = 0.063
		betax = 10.209
		alphay = 0.063
		betay = 10.776
		emitlim = 0.152 * 2*(order + 1) * 1e-6
		#xcenterpos = 0.0468
		xcenterpos = 0.0
		xcentermom = 0.00
		#ycenterpos = 0.0492
		ycenterpos = 0.0
		ycentermom = 0.00
		

		
		#====Injection and foil aperature============
		
		xmin = xcenterpos - 0.0085
		xmax = xcenterpos + 0.0085
		ymin = ycenterpos - 0.0080
		ymax = ycenterpos + 0.100
		
		
		xFunc=notRandom(args.x,args.px)
		yFunc=notRandom(args.y,args.py)
		lFunc=notRandom(args.z,args.pz)
		
		
		
		#print xmin
		#print xmax
		#print ymin
		#print ymax
		#=================Add the injection node and foil node==  ==============
		
		nparts = macrosperturn
		injectparams = (xmin, xmax, ymin, ymax)
		#injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
		#addTeapotInjectionNode(teapot_latt, 0., injectnode) 
		
		#print "(x,px,y,py,z,pz)= (%f,%f,%f,%f,%f,%f) " %(b.getSyncParticle().x(),b.getSyncParticle().px(),b.getSyncParticle().y(),b.getSyncParticle().py(),b.getSyncParticle().z(),b.getSyncParticle().pz())
		#print "gamma= %f"%b.getSyncParticle().gamma()
		#print "beta= %f"%b.getSyncParticle().beta()
		#print "momentum= %f"%b.getSyncParticle().momentum()
		inject = InjectParts(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
		inject.addParticles()
		
		thick = 400.0
		foil = TeapotFoilNode(xmin, xmax, ymin, ymax, thick, "Foil 1")
		scatterchoice = 0
		foil.setScatterChoice(scatterchoice)
		#addTeapotFoilNode(teapot_latt,0.000001,foil)
		
		#----------------------------------------------
		# Add one black absorber collimator to act like
		# an aperture
		#----------------------------------------------
		colllength = 0.00001
		ma = 9
		density_fac = 1.0
		shape = 1
		radius = 0.110
		
		collimator = TeapotCollimatorNode(colllength, ma, density_fac, shape, radius, 0., 0., 0., 0., pos = 0., name = "Collimator 1")
		#addTeapotCollimatorNode(teapot_latt, 0.5, collimator)
		
		#-----------------------------
		# Add RF Node
		#-----------------------------
		
		teapot_latt.initialize()
		
		
		#----------------------------------------------
		#make 2.5D space charge calculator
		#----------------------------------------------
		#set boundary
		nboundarypoints = 128
		n_freespacemodes = 32
		r_boundary = 0.220
		boundary = Boundary2D(nboundarypoints,n_freespacemodes,"Circle",r_boundary,r_boundary)
		
		sizeX = 64   #number of grid points in horizontal direction
		sizeY = 64  #number of grid points in vertical direction
		sizeZ = 1     #number of longitudinal slices in the 2.5D space charge solver
		calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)
		sc_path_length_min = 0.00000001
		scLatticeModifications.setSC2p5DAccNodes(teapot_latt, sc_path_length_min,calc2p5d, boundary)
		

		
		#-------------------------------
		#  Lattice is ready
		#-------------------------------
		
		#fileOut=open(args.fileName,'w')
		#fileOut.close()
		#myPrintNode=Print_Node("MyPrintNode",True,args.fileName)
		
		#nodes = teapot_latt.getNodes()
		#if args.nodeMonitor >= 0:
			#nodes[args.nodeMonitor].addChildNode(myPrintNode,AccNode.EXIT)
		#else:
		#	for node in nodes:
		#		node.addChildNode(myPrintNode,AccNode.EXIT)	
		#i = 0
		#for node in nodes:
		#	pass
		#	if node.getName().strip() == "DH_A12":
		#		print node.getName().strip()
		#		print node.getnParts()
				#node.setnParts(10)
		#	if args.printNodes==True:
		#		print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
		#		print "There are ", node.getNumberOfBodyChildren()," child nodes."
		#		i=i+1
		
		#================Do some turns===========================================
		for i in range(args.turns):
			teapot_latt.trackBunch(b, paramsDict)
			
		score = (b.x(0)-self.xTarget)**2 + (b.px(0)-self.pxTarget)**2 + (b.y(0)-self.yTarget)**2+(b.py(0)-self.pyTarget)**2+(b.z(0)-self.zTarget)**2+(b.pz(0)-self.dETarget)**2

		return score	


#scorer = MyScorer(0.005769,0.002069,0.001778,-0.000359,-0.003845,0.000000)
#scorer = MyScorer(0.004334,0.000192,0.001710,-0.000349,-0.004286,0.000000)
scorer = MyScorer(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000)

#searchAlgorithm   = RandomSearchAlgorithm()
searchAlgorithm = SimplexSearchAlgorithm()

#max_time = 0.05
max_time = 100
max_accuracy=1E-10
solverStopper = SolveStopperFactory.maxTimeStopper(max_time)
#solverStopper = SolveStopperFactory.maxAccuracyStopper(max_accuracy)

solver = Solver()
solver.setAlgorithm(searchAlgorithm)
solver.setStopper(solverStopper)

trialPoint = TrialPoint()
#trialPoint.addVariableProxy(VariableProxy(name = "x0", value = 1., step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x0", value = 1., step = 0))
trialPoint.addVariableProxy(VariableProxy(name = "x1", value = 1., step = 0))
trialPoint.addVariableProxy(VariableProxy(name = "x2", value = 1., step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x3", value = 1., step = 0.1))

solver.solve(scorer,trialPoint)

print "===== best score ========== fitting time = ", solver.getScoreboard().getRunTime()

bestScore = solver.getScoreboard().getBestScore()	
print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()

trialPoint = solver.getScoreboard().getBestTrialPoint()

print trialPoint.textDesciption()


