#############################################################
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

from orbit.teapot import GeneralDipole
from orbit.teapot import YDipole
from orbit.teapot import XDipole
from orbit.teapot import GeneralDipoleStrip
from orbit.teapot import GeneralDipoleNoStrip
from KevinPython.printNode import Print_Node
from KevinPython.calculateEmit import Calc_Emit

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

from orbit_utils import Function
from KevinPython.function_stripping import probabilityStripping

import argparse

print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--fileName", dest='fileName', default="outputAddMagnet.txt", help="file to print node info into")
parser.add_argument("--fileName2", dest='fileName2', default="outputAddMagnetEmitNoSpaceCharge.txt", help="file to print node info into")
parser.add_argument("--nParts",type=int, dest='nParts', default=260, help="number of particles")
#parser.add_argument("--turns",type=int, dest='turns', default=100, help="number of complete orbits")
parser.add_argument("--turns",type=int, dest='turns', default=1, help="number of complete orbits")
parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=35, help="What node to monitor")
#parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=35, help="What node to monitor")#35 should be currently used
#parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=37, help="What node to monitor")
parser.add_argument("--printNodes",type=bool, dest='printNodes', default=True, help="print node list")
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
parser.add_argument("--doNormalKickers",type=bool, dest='doNormalKickers', default=False, help="print node list")
#With Dipoles Neg optimizer
parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-1.047686, help="scaleChicane10")
parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-0.9969892 , help="scaleChicane11")
parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-1.00911, help="scaleChicane12")
parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-0.9775873, help="scaleChicane13")

parser.add_argument("--useChicaneScaleFile",type=bool, dest='useChicaneScaleFile', default=True, help="whether or not to use chicane scales from file")
#With Dipoles
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-0.9804517, help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-0.9934917 , help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-1.009009, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-1.042836, help="scaleChicane13")
#No Dipoles
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-1.015811, help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-1.000257, help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-0.9997909, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-1.008058, help="scaleChicane13")
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-1., help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-1., help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-1., help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-1., help="scaleChicane13")
#No Offset
parser.add_argument("--xOffset",type=float, dest='xOffset', default=0.0, help="x injection offset")
parser.add_argument("--pxOffset",type=float, dest='pxOffset', default=0.0, help="px injection offset")
parser.add_argument("--yOffset",type=float, dest='yOffset', default=0.0, help="y injection offset")
parser.add_argument("--pyOffset",type=float, dest='pyOffset', default=0, help="py injection offset")
parser.add_argument("--initialDriftLength",type=float, dest='initialDriftLength', default=0, help="x injection offset")

#No Dipoles
#parser.add_argument("--xOffset",type=float, dest='xOffset', default=0.095824, help="x injection offset")
#parser.add_argument("--pxOffset",type=float, dest='pxOffset', default=-0.010336, help="px injection offset")
#With Dipoles Neg optimizer
#parser.add_argument("--xOffset",type=float, dest='xOffset', default=0.098971, help="x injection offset")
#parser.add_argument("--pxOffset",type=float, dest='pxOffset', default=-0.013333, help="px injection offset")
#With Dipoles
#parser.add_argument("--xOffset",type=float, dest='xOffset', default=0.091194, help="x injection offset")
#parser.add_argument("--pxOffset",type=float, dest='pxOffset', default=-0.015937, help="px injection offset")
args = parser.parse_args()
#=====Main bunch parameters============
intensity = 7.8e13
#turns = 1000.0
turns = args.turns
#macrosperturn = 260
macrosperturn = args.nParts
macrosize = intensity/turns/macrosperturn

if args.useChicaneScaleFile==False:
	chicaneScale10=args.scaleChicane10
	chicaneScale11=args.scaleChicane11
	chicaneScale12=args.scaleChicane12
	chicaneScale13=args.scaleChicane13
else:
	chicaneScale10=1.
	chicaneScale11=1.
	chicaneScale12=1.
	chicaneScale13=1.	
nPartsChicane=4
outputDirectory="WasteBeamClosed"
for currentPart in range(nPartsChicane+1):
	inj_latt_start = teapot.TEAPOT_Ring()
	print "Read MAD."
	inj_latt_start.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers_ClosedWaste.LAT","RING")
	print "Lattice=",inj_latt_start.getName()," length [m] =",inj_latt_start.getLength()," nodes=",len(inj_latt_start.getNodes())
	
	
	#------------------------------
	#Initial Distribution Functions
	#------------------------------
	
	e_kin_ini = 1.3 # in [GeV]
	mass =  0.93827231 #0.939294    # in [GeV]
	gamma = (mass + e_kin_ini)/mass
	beta = math.sqrt(gamma*gamma - 1.0)/gamma
	print "relat. gamma=",gamma
	print "relat.  beta=",beta
	
	
	#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
	#(alphaX,betaX,emittX) = (-1.9620, 0.1831, 0.21)
	#(alphaY,betaY,emittY) = ( 1.7681, 0.1620, 0.21)
	(alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)
	
	(alphaX,betaX,emittX) = (.224, 10.5, 1.445)
	#(alphaX,betaX,emittX) = (22.4, 100, 1.445)
	(alphaY,betaY,emittY) = ( .224, 10.5, 1.445)
	#(alphaY,betaY,emittY) = ( 1.7681, 0.1620, 1.445)
	
	#alphaZ = -alphaZ
	
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
	#bunch_gen.setKinEnergy(e_kin_ini)
	
	#set the beam peak current in mA
	#bunch_gen.setBeamCurrent(38.0)
	
	bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = WaterBagDist3D)
	#bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
	#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)
	
	
	bunch_in.mass(mass) #mass
	bunch_in.macroSize(macrosize)
	energy = e_kin_ini # 1.0 #Gev
	bunch_in.getSyncParticle().kinEnergy(energy)
	bunch_in.charge(+1)
	
	paramsDict = {}
	lostbunch = Bunch()
	paramsDict["lostbunch"]=lostbunch
	paramsDict["bunch"]= bunch_in
	lostbunch.addPartAttr("LostParticleAttributes") 
	
	theEffLength=0.03*2
	#theEffLength=0.01
	fieldStrength=1.3
	fieldStrengthMin=.2
	cutLength=0.03
	
	sp = bunch_in.getSyncParticle()
	beta= sp.beta()
	gamma=sp.gamma()
	
	print "beta= %f"%beta
	print "gamma= %f"%gamma
	
	A1=2.47e-6
	A2=4.49e9
	
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
	
	print "accumlatedSum->getY(maxValue)"
	print accumlatedSum.getY(maxValue)
	print "CDF->getY(maxValue)"
	print CDF.getY(maxValue)
	print "InverseFunction.getY(maxValue)"
	print InverseFunction.getY(.954)
	print deltaxp_rigidity.getY(maxValue)
	print deltax_rigidity.getY(maxValue)
	
	fileOut=open("%s/emmit_DH11_3pre_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_DH11_3pre_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	
	fileOut=open("%s/emmit_DH12_3pre_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_DH12_3pre_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	
	myEmitNode_DH11_3pre=Calc_Emit("MyEmitNode_DH11_3pre_%d"%(currentPart),True,"%s/emmit_DH11_3pre_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_DH11_3pre=Print_Node("MyPrintNode_DH11_3pre_%d"%(currentPart),True,"%s/print_DH11_3pre_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_DH12_3pre=Calc_Emit("MyEmitNode_DH12_3pre_%d"%(currentPart),True,"%s/emmit_DH12_3pre_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_DH12_3pre=Print_Node("MyPrintNode_DH12_3pre_%d"%(currentPart),True,"%s/print_DH12_3pre_%d.txt"%(outputDirectory,currentPart))
	
	nodes2 = inj_latt_start.getNodes()
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
					node.addChildNode(myEmitNode_DH11_3pre,AccNode.BODY,currentPart)
					node.addChildNode(myPrintNode_DH11_3pre,AccNode.BODY,currentPart)
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
					node.addChildNode(myEmitNode_DH11_3pre,AccNode.BODY,0)
					node.addChildNode(myPrintNode_DH11_3pre,AccNode.BODY,0)
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
				node.addChildNode(myEmitNode_DH12_3pre,AccNode.BODY,3)
				node.addChildNode(myPrintNode_DH12_3pre,AccNode.BODY,3)			
				node.addChildNode(myDipole_DH_A12,AccNode.BODY,3)
			#drift before chicane 2
			if (node.getName().strip()=="DB12"):
				if (args.initialDriftLength>0):
					node.setLength(args.initialDriftLength)
	print "counterDH_A11=",counterDH_A11
	print "counterDH_A12=",counterDH_A12
	
	
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
	
	if args.useChicaneScaleFile:
		openedFile=open("%s/ChicaneScales_%d.txt"%(outputDirectory,currentPart),'r')
		line=openedFile.readline()
		print line
		theScales=line.split(",")
		chicaneScale10=-float(theScales[0].strip())
		chicaneScale11=-float(theScales[1].strip())
		chicaneScale12=-float(theScales[2].strip())
		chicaneScale13=-float(theScales[3].strip())
	strength_chicane10 = -0.041456*chicaneScale10
	strength_chicane11 = 0.052434*chicaneScale11
	strength_chicane12 = 0.0298523*chicaneScale12
	strength_chicane13 = -0.0398609*chicaneScale13
	
	#strength_chicane12 = -0.041456*chicaneScale10
	#strength_chicane13 = 0.052434*chicaneScale11
	#strength_chicane10 = 0.0298523*chicaneScale12
	#strength_chicane11 = -0.0398609*chicaneScale13

		
	
	lattlength = inj_latt_start.getLength()
	sp = bunch_in.getSyncParticle()
	kickerwave = rootTWaveform(sp, lattlength, duration, startamp, endamp)
	chicanewave = flatTopWaveform(1.0)
	
	i = 0
	path_length=0
	print "inj_latt"
	thick = 400.0
	#foil = TeapotFoilNode(xmin, xmax, ymin, ymax, thick, "Foil 1")
	foil = TeapotFoilNode(-100, 100, -100, 100, thick, "Foil 1")
	scatterchoice = 0
	foil.setScatterChoice(scatterchoice)
	#addTeapotFoilNode(inj_latt_start,5.3503,foil)	
	#addTeapotFoilNode(inj_latt_start,5.35,foil)	
	nodes = inj_latt_start.getNodes()
	chicane10 = nodes[29]
	chicane11 = nodes[31]
	chicane12 = nodes[34]
	chicane13 = nodes[36]
	chicane10.setParam("kx", strength_chicane10)
	chicane11.setParam("kx", strength_chicane11)
	chicane12.setParam("kx", strength_chicane12)
	chicane13.setParam("kx", strength_chicane13)
	chicane10.setWaveform(chicanewave)
	chicane11.setWaveform(chicanewave)
	chicane12.setWaveform(chicanewave)
	chicane13.setWaveform(chicanewave)
	chicane11.setWaveform(chicanewave)	
	for node in nodes:
		pass
		if node.getName().strip() == "DH_A12":
			print node.getName().strip()
			print node.getnParts()
			#node.setnParts(10)
		if args.printNodes==True:
			path_length=path_length+node.getLength()
			print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%inj_latt_start.getNodePositionsDict()[node], " path_length= ",path_length
			#print "There are ", node.getNumberOfBodyChildren()," child nodes."
			i=i+1
	
			
	#nodes = inj_latt.getNodes()
	#chicane11 = nodes[1]
	#chicane12 = nodes[4]
	
	
	#chicane11.setParam("kx", strength_chicane11)
	#chicane12.setParam("kx", strength_chicane12)
	
	#chicane11.setWaveform(chicanewave)
	#chicane12.setWaveform(chicanewave)
	#for node in kickernode:
	#print "node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
	#print "There are ", node.getNumberOfBodyChildren()," child nodes."
	
	
	
	
	
	#-------------------------------
	#  Lattice is ready
	#-------------------------------
	
	
	#myPrintNode=Print_Node("MyPrintNode",True,args.fileName)
	fileOut=open("%s/print_begQ_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_beg_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_beg_DH11_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_beg_DH10_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_postS_DH11_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/print_end_DH11_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_end_DH10_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/print_beg_b23_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_mid_b23_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_end_b23_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_beg_DH12_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_postS_DH12_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/print_end_DH12_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_beg_DH13_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/print_end_DH13_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/print_end_DB_WASTE_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()		
	myPrintNode_begQ=Print_Node("MyPrintNode_BegQ_%d"%(currentPart),True,"%s/print_begQ_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_beg=Print_Node("MyPrintNode_Beg_%d"%(currentPart),True,"%s/print_beg_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_beg_DH11=Print_Node("MyPrintNode_beg_DH11_%d"%(currentPart),True,"%s/print_beg_DH11_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_beg_DH10=Print_Node("MyPrintNode_beg_DH10_%d"%(currentPart),True,"%s/print_beg_DH10_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_postS_DH11=Print_Node("MyPrintNode_postS_DH11_%d"%(currentPart),True,"%s/print_postS_DH11_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_end_DH11=Print_Node("MyPrintNode_end_DH11_%d"%(currentPart),True,"%s/print_end_DH11_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_end_DH10=Print_Node("MyPrintNode_end_DH10_%d"%(currentPart),True,"%s/print_end_DH10_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_beg_b23=Print_Node("MyPrintNode_beg_b23_%d"%(currentPart),True,"%s/print_beg_b23_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_mid_b23=Print_Node("MyPrintNode_mid_b23_%d"%(currentPart),True,"%s/print_mid_b23_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_end_b23=Print_Node("MyPrintNode_end_b23_%d"%(currentPart),True,"%s/print_end_b23_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_beg_DH12=Print_Node("MyPrintNode_beg_DH12_%d"%(currentPart),True,"%s/print_beg_DH12_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_postS_DH12=Print_Node("MyPrintNode_postS_DH12_%d"%(currentPart),True,"%s/print_postS_DH12_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_end_DH12=Print_Node("MyPrintNode_end_DH12_%d"%(currentPart),True,"%s/print_end_DH12_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_beg_DH13=Print_Node("MyPrintNode_beg_DH13_%d"%(currentPart),True,"%s/print_beg_DH13_%d.txt"%(outputDirectory,currentPart))
	myPrintNode_end_DH13=Print_Node("MyPrintNode_end_DH13_%d"%(currentPart),True,"%s/print_end_DH13_%d.txt"%(outputDirectory,currentPart))	
	myPrintNode_end_DB_WASTE=Print_Node("MyPrintNode_end_DB_WASTE_%d"%(currentPart),True,"%s/print_end_DB_WASTE_%d.txt"%(outputDirectory,currentPart))	

	fileOut=open("%s/emmit_begQ_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_beg_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_beg_DH11_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_beg_DH10_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/emmit_postS_DH11_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/emmit_end_DH11_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_end_DH10_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_beg_b23_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_mid_b23_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_end_b23_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_beg_DH12_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_end_DH12_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_postS_DH12_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_beg_DH13_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()
	fileOut=open("%s/emmit_end_DH13_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	fileOut=open("%s/emmit_end_DB_WASTE_%d.txt"%(outputDirectory,currentPart),'w')
	fileOut.close()	
	myEmitNode_begQ=Calc_Emit("MyEmitNode_BegQ_%d"%(currentPart),True,"%s/emmit_begQ_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_beg=Calc_Emit("MyEmitNode_Beg_%d"%(currentPart),True,"%s/emmit_beg_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_beg_DH11=Calc_Emit("myEmitNode_beg_DH11_%d"%(currentPart),True,"%s/emmit_beg_DH11_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_beg_DH10=Calc_Emit("myEmitNode_beg_DH10_%d"%(currentPart),True,"%s/emmit_beg_DH10_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_postS_DH11=Calc_Emit("myEmitNode_postS_DH11_%d"%(currentPart),True,"%s/emmit_postS_DH11_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_end_DH11=Calc_Emit("myEmitNode_end_DH11_%d"%(currentPart),True,"%s/emmit_end_DH11_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_end_DH10=Calc_Emit("myEmitNode_end_DH10_%d"%(currentPart),True,"%s/emmit_end_DH10_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_beg_b23=Calc_Emit("myEmitNode_beg_b23_%d"%(currentPart),True,"%s/emmit_beg_b23_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_mid_b23=Calc_Emit("myEmitNode_mid_b23_%d"%(currentPart),True,"%s/emmit_mid_b23_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_end_b23=Calc_Emit("myEmitNode_end_b23_%d"%(currentPart),True,"%s/emmit_end_b23_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_beg_DH12=Calc_Emit("myEmitNode_beg_DH12_%d"%(currentPart),True,"%s/emmit_beg_DH12_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_postS_DH12=Calc_Emit("myEmitNode_postS_DH12_%d"%(currentPart),True,"%s/emmit_postS_DH12_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_end_DH12=Calc_Emit("myEmitNode_end_DH12_%d"%(currentPart),True,"%s/emmit_end_DH12_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_beg_DH13=Calc_Emit("myEmitNode_beg_DH13_%d"%(currentPart),True,"%s/emmit_beg_DH13_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_end_DH13=Calc_Emit("myEmitNode_end_DH13_%d"%(currentPart),True,"%s/emmit_end_DH13_%d.txt"%(outputDirectory,currentPart))
	myEmitNode_end_DB_WASTE=Calc_Emit("myEmitNode_end_DB_WASTE_%d"%(currentPart),True,"%s/emmit_end_DB_WASTE_%d.txt"%(outputDirectory,currentPart))
	
	nodes = inj_latt_start.getNodes()
	nodes[0].addChildNode(myPrintNode_begQ,AccNode.ENTRANCE)
	nodes[0].addChildNode(myEmitNode_begQ,AccNode.ENTRANCE)
	#nodes[args.nodeMonitor].addChildNode(myEmitNode,AccNode.EXIT)
	i = 0
	path_length=0
	for node in nodes:
		pass
		if node.getName().strip() == "DB12":
			node.addChildNode(myPrintNode_beg,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_beg,AccNode.ENTRANCE)
		if node.getName().strip() == "DB23":
			#node.setnParts(2)
			node.addChildNode(myPrintNode_beg_b23,AccNode.ENTRANCE)
			node.addChildNode(myPrintNode_mid_b23,AccNode.BODY,1)
			node.addChildNode(myPrintNode_end_b23,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_beg_b23,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_mid_b23,AccNode.BODY,1)
			node.addChildNode(myEmitNode_end_b23,AccNode.EXIT)
			
			if currentPart is nPartsChicane:
				node.addChildNode(myEmitNode_postS_DH11,AccNode.BODY,0)	
				node.addChildNode(myPrintNode_postS_DH11,AccNode.BODY,0)
			
		if node.getName().strip() == "DH_A10":
			node.addChildNode(myPrintNode_beg_DH10,AccNode.ENTRANCE)
			node.addChildNode(myPrintNode_end_DH10,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_beg_DH10,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_end_DH10,AccNode.EXIT)			
		if node.getName().strip() == "DH_A11":
			node.addChildNode(myPrintNode_beg_DH11,AccNode.ENTRANCE)
			node.addChildNode(myPrintNode_end_DH11,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_beg_DH11,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_end_DH11,AccNode.EXIT)
			
			if currentPart is not nPartsChicane:
				node.addChildNode(myEmitNode_postS_DH11,AccNode.BODY,currentPart)
				node.addChildNode(myPrintNode_postS_DH11,AccNode.BODY,currentPart)

		if node.getName().strip() == "DB23":
			#node.setnParts(2)
			node.addChildNode(myPrintNode_beg_b23,AccNode.ENTRANCE)
			node.addChildNode(myPrintNode_mid_b23,AccNode.BODY,1)
			node.addChildNode(myPrintNode_end_b23,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_beg_b23,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_mid_b23,AccNode.BODY,1)
			node.addChildNode(myEmitNode_end_b23,AccNode.EXIT)
			pass
				
		if node.getName().strip() == "DH_A12":
			#print node.getName().strip()
			#print node.getnParts()
			node.addChildNode(myPrintNode_beg_DH12,AccNode.ENTRANCE)
			node.addChildNode(myPrintNode_end_DH12,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_beg_DH12,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_end_DH12,AccNode.EXIT)
			if currentPart is not nPartsChicane:
				node.addChildNode(myEmitNode_postS_DH12,AccNode.BODY,3)	
				node.addChildNode(myPrintNode_postS_DH12,AccNode.BODY,0)
			else:
				node.addChildNode(myEmitNode_postS_DH12,AccNode.BODY,3)	
				node.addChildNode(myPrintNode_postS_DH12,AccNode.BODY,0)				
		if node.getName().strip() == "DH_A13":
			#print node.getName().strip()
			#print node.getnParts()
			node.addChildNode(myPrintNode_beg_DH13,AccNode.ENTRANCE)
			node.addChildNode(myPrintNode_end_DH13,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_beg_DH13,AccNode.ENTRANCE)
			node.addChildNode(myEmitNode_end_DH13,AccNode.EXIT)		
			#node.setnParts(10)
			
		if node.getName().strip() == "DB_Waste":
			node.addChildNode(myPrintNode_end_DB_WASTE,AccNode.EXIT)		
			
			node.addChildNode(myEmitNode_end_DB_WASTE,AccNode.EXIT)				
	
	#================Do some turns===========================================
	
	xOffset=args.xOffset
	pxOffset=args.pxOffset
	yOffset=args.yOffset
	pyOffset=args.pyOffset
	for i in range(bunch_in.getSize()):
		bunch_in.x(i,bunch_in.x(i)+xOffset)
		bunch_in.px(i,bunch_in.px(i)+pxOffset)
		bunch_in.y(i,bunch_in.y(i)+yOffset)
		bunch_in.py(i,bunch_in.py(i)+pyOffset)		
	#track through drift after 3rd chicane and foil at end
	inj_latt_start.trackBunch(bunch_in, paramsDict)

	
	#===========Dump bunch infomration=======================================
	#bunch_pyorbit_to_orbit(inj_latt.getLength(), bunch_in, "mainbunch.dat")
	#bunch_pyorbit_to_orbit(inj_latt.getLength(), lostbunch, "lostbunch.dat")
	print "Stop."



