#############################################################
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

from orbit.teapot import GeneralDipole
from orbit.teapot import YDipole
from orbit.teapot import XDipole
from orbit.teapot import GeneralDipoleStrip
from orbit.teapot import GeneralDipoleStripSeperateField
from orbit.teapot import GeneralDipoleNoStripSeperateField
from KevinPython.printNode import Print_Node
from KevinPython.calculateEmit import Calc_Emit

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

from orbit_utils import Function
from KevinPython.function_stripping import probabilityStripping
from KevinPython.function_strippingIncludeChicaneField import probabilityStrippingWithChicane
from orbit.teapot import addDipoleStripperNode

import argparse

print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--fileName", dest='fileName', default="outputAddMagnet.txt", help="file to print node info into")
parser.add_argument("--fileName2", dest='fileName2', default="outputAddMagnetEmitNoSpaceCharge.txt", help="file to print node info into")
parser.add_argument("--nParts",type=int, dest='nParts', default=10000, help="number of particles")
#parser.add_argument("--turns",type=int, dest='turns', default=100, help="number of complete orbits")
parser.add_argument("--turns",type=int, dest='turns', default=1, help="number of complete orbits")
#parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=35, help="What node to monitor")#35 should be currently used
#parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=37, help="What node to monitor")
parser.add_argument("--printNodes",type=bool, dest='printNodes', default=True, help="print node list")
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="use stripper dipoles")

parser.add_argument("--useChicaneScaleFile",type=bool, dest='useChicaneScaleFile', default=True, help="whether or not to use chicane scales from file")
parser.add_argument("--usePrintNode",type=bool, dest='usePrintNode', default=False, help="whether or not to print bunch info to file")
parser.add_argument("--pencilBeam",type=bool, dest='pencilBeam', default=False, help="Use a single macroparticle beam")
parser.add_argument("--addChicaneFieldToStripper",type=bool, dest='addChicaneFieldToStripper', default=True, help="Include the chicane fields in the stripper if stripper is inside chicane")
parser.add_argument("--bunchFromFile",type=bool, dest='bunchFromFile', default=False, help="Create bunch reading particles from file")
parser.add_argument("--useSecondaryFoil",type=bool, dest='useSecondaryFoil', default=False, help="use secondary foil in lattice")
parser.add_argument("--bunchFromFileName", dest='bunchFromFileName', default="InitialBunches/print_beg_0.txt", help="What File to read bunch from")
parser.add_argument("--outputDirectory", dest='outputDirectory', default="WasteBeamSplitGeneralClosedBeamNewY_Inject3", help="Where to put output")

parser.add_argument("--stripperLength1",type=float, dest='stripperLength1', default=0.06, help="length of first stripper dipole")
parser.add_argument("--stripperStrengthMax1",type=float, dest='stripperStrengthMax1', default=1.3, help="Maximum field Strength of first stripper dipole")
parser.add_argument("--stripperStrengthMin1",type=float, dest='stripperStrengthMin1', default=0.2, help="Minimum field Strength of first stripper dipole")
parser.add_argument("--cutLength1",type=float, dest='cutLength1', default=0.03, help="length the field ramps up linearly from min to max strength")
parser.add_argument("--fieldDirection1",type=float, dest='fieldDirection1', default=math.pi/2., help="The direction of the field of the first stripper dipole (0=positive x, Pi/2=positive y)")

parser.add_argument("--stripperLength2",type=float, dest='stripperLength2', default=0.06, help="length of second stripper dipole")
parser.add_argument("--stripperStrengthMax2",type=float, dest='stripperStrengthMax2', default=1.3, help="Maximum field Strength of second stripper dipole")
parser.add_argument("--stripperStrengthMin2",type=float, dest='stripperStrengthMin2', default=0.2, help="Minimum field Strength of second stripper dipole")
parser.add_argument("--cutLength2",type=float, dest='cutLength2', default=0.03, help="length the field ramps up linearly from min to max strength")
parser.add_argument("--fieldDirection2",type=float, dest='fieldDirection2', default=math.pi/2., help="The direction of the field of the second stripper dipole (0=positive x, Pi/2=positive y)")

parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-1., help="scaleChicane10")
parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-1., help="scaleChicane11")
parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-1., help="scaleChicane12")
parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-1., help="scaleChicane13")

parser.add_argument("--xOffset",type=float, dest='xOffset', default=0., help="x injection offset")
parser.add_argument("--pxOffset",type=float, dest='pxOffset', default=0., help="px injection offset")
parser.add_argument("--yOffset",type=float, dest='yOffset', default=0., help="y injection offset")
parser.add_argument("--pyOffset",type=float, dest='pyOffset', default=0., help="py injection offset")


args = parser.parse_args()
#=====Main bunch parameters============
intensity = 7.8e13
turns = args.turns
#macrosperturn = 260
macrosperturn = args.nParts
macrosize = intensity/turns/macrosperturn

#where to pull chicane scale values from if useChicaneScales is true.

#outputDirectoryChicaneScales="WasteBeamSplitGeneralNewStripperChicaneFieldAddedCleanOpposite"
outputDirectoryChicaneScales="InjectBeam3"
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

#the default chicane kick strength array
chicaneStrengthArray=[-0.041456,0.052434,0.0298523,-0.0398609]

usePrintNode=args.usePrintNode
#this sets how to divide up chicane2/11 in terms of where 1st stripper is placed.
nPartsChicane=1
#nPartsChicane=0
nPartsChicane2=0
outputDirectory=args.outputDirectory
if not os.path.exists(outputDirectory):
	os.mkdir(outputDirectory)
	
#loop over all configurations of where 1st stripper dipole with be placed
#currentPart=0 places it at begginging of chicane2/11
#currentPart=nPartsChicane places it immediately after chicane2/11
#0<currentPart<nPartsChicane places it currentPart/nPartsChicane fractionally into chicane2/11
for currentPart in range(1,nPartsChicane+1):
	for currentPart2 in range(-1,nPartsChicane2+1):
		inj_latt_start = teapot.TEAPOT_Ring()
		print "Read MAD."
		#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
		inj_latt_start.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers_ClosedWaste.LAT","RING")
		print "Lattice=",inj_latt_start.getName()," length [m] =",inj_latt_start.getLength()," nodes=",len(inj_latt_start.getNodes())
	
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
		
		nodes = inj_latt_start.getNodes()
		hkick10 = nodes[10]
		vkick10 = nodes[12]
		hkick11	= nodes[14]
		vkick11 = nodes[16]
		#vkick12 = nodes[49]
		#hkick12 = nodes[51]
		#vkick13 = nodes[53]
		#hkick13 = nodes[55]
		
		vkick10.setParam("ky", strength_vkicker10)
		hkick10.setParam("kx", strength_hkicker10)
		vkick11.setParam("ky", strength_vkicker11)
		hkick11.setParam("kx", strength_hkicker11)
		#vkick12.setParam("ky", strength_vkicker12)
		#hkick12.setParam("kx", strength_hkicker12)
		#vkick13.setParam("ky", strength_vkicker13)
		#hkick13.setParam("kx", strength_hkicker13)
		
		vkick10.setWaveform(kickerwave)
		hkick10.setWaveform(kickerwave)
		vkick11.setWaveform(kickerwave)
		hkick11.setWaveform(kickerwave)
		#vkick12.setWaveform(kickerwave)
		#hkick12.setWaveform(kickerwave)
		#vkick13.setWaveform(kickerwave)
		#hkick13.setWaveform(kickerwave)
			
		inj_latt_start.initialize()
		
		#-----------
		#Initial Distribution Functions
		#------------------------------
		
		e_kin_ini = 1.3 # in [GeV]
		mass =  0.93827231 #0.939294    # in [GeV]
		gamma = (mass + e_kin_ini)/mass
		beta = math.sqrt(gamma*gamma - 1.0)/gamma
		print "relat. gamma=",gamma
		print "relat.  beta=",beta
		
		
		#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
		(alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)
		
		(alphaX,betaX,emittX) = (.224, 10.5, 1.445)
		(alphaY,betaY,emittY) = (.224, 10.5, 1.445)
		
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
		
		bunch_in=Bunch()
		#use pencil beam
		if args.pencilBeam: 
			bunch_in.addParticle(0,0,0,0,0,0)
		#read initial bunch in from file
		elif args.bunchFromFile:
			openedFile=open("%s"%(args.bunchFromFileName),'r')
			lines=openedFile.readlines()
			for line in lines:
				coordToken=line.split("=")[1]
				coordArr=coordToken.split(",")
				bunch_in.addParticle(float(coordArr[0].strip('()\n ')),float(coordArr[1].strip('()\n ')),float(coordArr[2].strip('()\n ')),float(coordArr[3].strip('()\n ')),float(coordArr[4].strip('()\n ')),float(coordArr[5].strip('()\n ')))
		#generate initial bunch
		else:
			bunch_in = bunch_gen.getBunch(nParticles = args.nParts, distributorClass = WaterBagDist3D)
			#bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
			#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)
		
		
		bunch_in.mass(mass) #mass
		bunch_in.macroSize(macrosize)
		energy = e_kin_ini # 1.0 #Gev
		bunch_in.getSyncParticle().kinEnergy(energy)
		bunch_in.charge(1)
	
		#apply offset to initial beam
		xOffset=args.xOffset
		pxOffset=args.pxOffset
		yOffset=args.yOffset
		pyOffset=args.pyOffset	
		#if reading bunch from file the offset should already have been added
		if not args.bunchFromFile:
			for i in range(bunch_in.getSize()):
				bunch_in.x(i,bunch_in.x(i)+xOffset)
				bunch_in.px(i,bunch_in.px(i)+pxOffset)
				bunch_in.y(i,bunch_in.y(i)+yOffset)
				bunch_in.py(i,bunch_in.py(i)+pyOffset)	
				
		paramsDict = {}
		#create failed to strip bunches
	
		lostbunch = Bunch()
		paramsDict["lostbunch"]=lostbunch
		paramsDict["bunch"]= bunch_in
		lostbunch.addPartAttr("LostParticleAttributes") 
		
		#erase files for writing emmittance info immediately prior to stripping dipoles and immediately after stripping
		#erase files for writing emmittance info immediately prior to stripping dipoles and immediately after stripping
		fileOut=open("%s/emmit_DH11_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_DH12_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		if usePrintNode:
			fileOut=open("%s/print_DH11_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()		
			fileOut=open("%s/print_DH12_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
		
		myEmitNode_DH11_3pre=Calc_Emit("MyEmitNode_DH11_3pre_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_DH11_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_DH12_3pre=Calc_Emit("MyEmitNode_DH12_3pre_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_DH12_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		
		if usePrintNode:
			myPrintNode_DH12_3pre=Print_Node("MyPrintNode_DH12_3pre_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_DH12_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_DH11_3pre=Print_Node("MyPrintNode_DH11_3pre_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_DH11_3pre_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			
		fileOut=open("%s/emmit_postS_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()	
		fileOut=open("%s/emmit_postS_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()	
		myEmitNode_postS_DH11=Calc_Emit("myEmitNode_postS_DH11_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_postS_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))	
		myEmitNode_postS_DH12=Calc_Emit("myEmitNode_postS_DH12_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_postS_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		
		
		print "inj_latt"
		#create secondary foil
		thick = 400.0
		foil = TeapotFoilNode(-100, 100, -100, 100, thick, "Foil 2")
		scatterchoice = 0
		foil.setScatterChoice(scatterchoice)
		#place foil at end of first lattice (ie after drift DB34) if secondary foil being used
		if args.useSecondaryFoil:
			addTeapotFoilNode(inj_latt_start,inj_latt_start.getLength(),foil)	
	
		nodes = inj_latt_start.getNodes()
		
		#set the strength of the chicane kicks
		if args.useChicaneScaleFile:
			openedFile=open("%s/ChicaneScales_%d_%d_%d_%d.txt"%(outputDirectoryChicaneScales,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'r')
			line=openedFile.readline()
			print line
			theScales=line.split(",")
			chicaneScale10=-float(theScales[0].strip())
			chicaneScale11=-float(theScales[1].strip())
			chicaneScale12=-float(theScales[2].strip())
			chicaneScale13=-float(theScales[3].strip())
			openedFile.close()
			
		strength_chicane10 = chicaneStrengthArray[0]*chicaneScale10
		strength_chicane11 = chicaneStrengthArray[1]*chicaneScale11
		strength_chicane12 = chicaneStrengthArray[2]*chicaneScale12
		strength_chicane13 = chicaneStrengthArray[3]*chicaneScale13
			
		#find chicanes
		chicane10 = None
		chicane11 = None
		chicane12 = None
		chicane13 = None	
		driftDB23=None
		for index in range(len(inj_latt_start.getNodes())):
			if nodes[index].getName().strip() == "DH_A10":
				chicane10=nodes[index]
			elif nodes[index].getName().strip() == "DH_A11":
				chicane11=nodes[index]
			elif nodes[index].getName().strip() == "DH_A12":
				chicane12=nodes[index]
			elif nodes[index].getName().strip() == "DH_A13":
				chicane13=nodes[index]	
			elif nodes[index].getName().strip() == "DB23":
				driftDB23=nodes[index]			
		chicanewave = flatTopWaveform(1.0)
		print chicane10.getName()
		print strength_chicane10
		chicane10.setParam("kx", strength_chicane10)
		chicane11.setParam("kx", strength_chicane11)
		chicane12.setParam("kx", strength_chicane12)
		chicane13.setParam("kx", strength_chicane13)
		chicane10.setWaveform(chicanewave)
		chicane11.setWaveform(chicanewave)
		chicane12.setWaveform(chicanewave)
		chicane13.setWaveform(chicanewave)
	
		#calculate where to place 1st stripper dipole
		position=-100.
		if currentPart==-1:
			#position =inj_latt_start.getNodePositionsDict()[nodes[0]][1]-theEffLength
			position =inj_latt_start.getNodePositionsDict()[chicane11][0]-theEffLength
		elif currentPart==0:
			position =inj_latt_start.getNodePositionsDict()[chicane11][0]
		elif currentPart is nPartsChicane:
			position =inj_latt_start.getNodePositionsDict()[chicane11][1]
		else :
			position =inj_latt_start.getNodePositionsDict()[chicane11][0]+chicane11.getLength()*currentPart/nPartsChicane
			
		#create stripping dipole magnetic field
		theEffLength=args.stripperLength1
		fieldStrength=args.stripperStrengthMax1
		fieldStrengthMin=args.stripperStrengthMin1
		#fieldStrength=.4
		#fieldStrengthMin=.4
		cutLength=args.cutLength1
		#fieldDirection=math.pi/2.
		fieldDirection=args.fieldDirection1
		
		sp = bunch_in.getSyncParticle()
		beta= sp.beta()
		gamma=sp.gamma()
		
		print "beta= %f"%beta
		print "gamma= %f"%gamma
		c=299792458
		
		rigidity= sp.momentum()/(c/math.pow(10.,9))
		print "rigidity= %f"%rigidity	
		
		if args.doDipoleKickers:
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + theEffLength
			(node_start_ind,node_stop_ind,z,ind) = (-1,-1, 0., 0)
			for nodeCurrent in inj_latt_start.getNodes():
				if(position_start >= z and position_start <= z + nodeCurrent.getLength()):
					node_start_ind = ind
				if(position_stop > z and position_stop <= z + nodeCurrent.getLength()):
					node_stop_ind = ind
				ind += 1
				z += nodeCurrent.getLength()	
				
			if node_start_ind!=node_stop_ind:
				#the stripping dipole spans more than 1 node
				print "something is going to be broken"
				sys.exit(0)
			nodeCurrent=inj_latt_start.getNodes()[node_start_ind]
			xkickerField=0.
			ykickerField=0.
			#nothing to change because no field in drift
			if(isinstance(nodeCurrent,DriftTEAPOT)):
				print "stripper dipole is in drift"
			#in kicker so add kicker field to stripping field
			elif (isinstance(nodeCurrent,KickTEAPOT)):
				print "stripper dipole is in kick node"
				if args.addChicaneFieldToStripper:
					#compute field to create kick
					length=nodeCurrent.getLength()
					kx=nodeCurrent.getParam("kx")
					ykickerField=-kx*rigidity/length
					ky=nodeCurrent.getParam("ky")
					xkickerField=ky*rigidity/length
					
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
			#field starts at fieldStrengthMin, ramps up to fieldStrength linearly from x=0 to x=cutlength. THen from x=cutlength to end of magnet it is constant at fieldStrength
			def pieceWiseField2(x):
				if x<cutLength :
					return (fieldStrength-fieldStrengthMin)/cutLength*x+fieldStrengthMin
				elif x>=cutLength:
					return fieldStrength
				pass
			magneticFieldx= Function()
			magneticFieldy= Function()
			#number of parts to break stripper dipole into
			n=1000
			maxValue=theEffLength
			step=maxValue/n
			
			#actually creates the field functions
			for i in range(n):
				x = step*i;
				#y = constantField(x)
				y = pieceWiseField2(x)
				magneticFieldx.add(x,y*math.cos(fieldDirection)+xkickerField)
				magneticFieldy.add(x,y*math.sin(fieldDirection)+ykickerField)
			
			myDipole_DH_A11=GeneralDipoleNoStripSeperateField(magneticFieldx,magneticFieldy,n,maxValue,gamma,beta,"Dipole_DH_A11")	
			myDipole_DH_A11.addChildNode(myEmitNode_DH11_3pre,AccNode.ENTRANCE)
			myDipole_DH_A11.addChildNode(myEmitNode_postS_DH11,AccNode.EXIT)
			
			addDipoleStripperNode(inj_latt_start,position,myDipole_DH_A11)
						
		#add second stripper dipole without stripping
		if args.doDipoleKickers:
			theEffLength=args.stripperLength2
			fieldStrength=args.stripperStrengthMax2
			fieldStrengthMin=args.stripperStrengthMin2
			#fieldStrength=.4
			#fieldStrengthMin=.4
			cutLength=args.cutLength2
			#fieldDirection=math.pi/2.
			fieldDirection=args.fieldDirection2
			position=-100.
			#place second stripper 5/6 of the way into chicane3/12. temporary position for consistency
			#position =inj_latt_start.getNodePositionsDict()[chicane12][0]+chicane12.getLength()*5./6.	
			if currentPart2==-1:
				#position =inj_latt_start.getNodePositionsDict()[driftDB23][1]-theEffLength
				position =inj_latt_start.getNodePositionsDict()[chicane12][0]-theEffLength
			elif currentPart2==0:
				position =inj_latt_start.getNodePositionsDict()[chicane12][0]
			elif currentPart2 is nPartsChicane2:
				position =inj_latt_start.getNodePositionsDict()[chicane12][1]
			else :
				position =inj_latt_start.getNodePositionsDict()[chicane12][0]+chicane12.getLength()*currentPart2/nPartsChicane2		
			#check if we are in kicker or drift
			position_start = position
			position_stop = position + theEffLength
			(node_start_ind,node_stop_ind,z,ind) = (-1,-1, 0., 0)
			for nodeCurrent in inj_latt_start.getNodes():
				if(position_start >= z and position_start <= z + nodeCurrent.getLength()):
					node_start_ind = ind
				if(position_stop > z and position_stop <= z + nodeCurrent.getLength()):
					node_stop_ind = ind
				ind += 1
				z += nodeCurrent.getLength()	
				
			if node_start_ind!=node_stop_ind:
				print "something is going to be broken2"
				sys.exit(0)
			nodeCurrent=inj_latt_start.getNodes()[node_start_ind]
			xkickerField=0.
			ykickerField=0.
			#nothing to change
			if(isinstance(nodeCurrent,DriftTEAPOT)):
				print "stripper dipole is in drift"
			elif (isinstance(nodeCurrent,KickTEAPOT)):
				print "stripper dipole is in kick node"
				if args.addChicaneFieldToStripper:
					#compute field to create kick
					length=nodeCurrent.getLength()
					kx=nodeCurrent.getParam("kx")		
					ykickerField=-kx*rigidity/length
					ky=nodeCurrent.getParam("ky")
					xkickerField=ky*rigidity/length
					
						
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
			magneticFieldx= Function()
			magneticFieldy= Function()
			n=1000
			maxValue=theEffLength
			step=maxValue/n
			
			for i in range(n):
				x = step*i;
				#y = constantField(x)
				y = pieceWiseField2(x)
				magneticFieldx.add(x,y*math.cos(fieldDirection)+xkickerField)
				magneticFieldy.add(x,y*math.sin(fieldDirection)+ykickerField)
					
			myDipole_DH_A12=GeneralDipoleNoStripSeperateField(magneticFieldx,magneticFieldy,n,maxValue,gamma,beta,"Dipole_DH_A12")
			myDipole_DH_A12.addChildNode(myEmitNode_DH12_3pre,AccNode.ENTRANCE)
			myDipole_DH_A12.addChildNode(myEmitNode_postS_DH12,AccNode.EXIT)
			
			addDipoleStripperNode(inj_latt_start,position,myDipole_DH_A12)	
			
		#print lattice info
		i = 0
		path_length=0
		for node in nodes:
			if args.printNodes==True:
				path_length=path_length+node.getLength()
				print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%inj_latt_start.getNodePositionsDict()[node], " path_length= ",path_length
				i=i+1	
	
	
		#------------------------------
		#  Lattice is ready
		#-------------------------------
		
		#create all the other emmit calculation nodes at various locations along lattice
		if usePrintNode:
			fileOut=open("%s/print_beg_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_begQ_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()				
			fileOut=open("%s/print_beg_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_postS_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()	
			fileOut=open("%s/print_end_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_beg_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_mid_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_end_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_beg_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_postS_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()	
			fileOut=open("%s/print_end_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_beg_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()
			fileOut=open("%s/print_end_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()	
			fileOut=open("%s/print_end_DB_WASTE_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
			fileOut.close()		
			myPrintNode_beg=Print_Node("MyPrintNode_Beg_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_beg_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_begQ=Print_Node("MyPrintNode_BegQ_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_begQ_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_beg_DH11=Print_Node("MyPrintNode_beg_DH11_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_beg_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_postS_DH11=Print_Node("MyPrintNode_postS_DH11_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_postS_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_end_DH11=Print_Node("MyPrintNode_end_DH11_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_end_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_beg_b23=Print_Node("MyPrintNode_beg_b23_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_beg_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_mid_b23=Print_Node("MyPrintNode_mid_b23_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_mid_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_end_b23=Print_Node("MyPrintNode_end_b23_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_end_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_beg_DH12=Print_Node("MyPrintNode_beg_DH12_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_beg_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_postS_DH12=Print_Node("MyPrintNode_postS_DH12_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_postS_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_end_DH12=Print_Node("MyPrintNode_end_DH12_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_end_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_beg_DH13=Print_Node("MyPrintNode_beg_DH13_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_beg_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
			myPrintNode_end_DH13=Print_Node("MyPrintNode_end_DH13_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_end_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))	
			myPrintNode_end_DB_WASTE=Print_Node("MyPrintNode_end_DB_WASTE_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/print_end_DB_WASTE_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))	
		
		fileOut=open("%s/emmit_beg_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_begQ_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()			
		fileOut=open("%s/emmit_beg_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_end_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_beg_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_mid_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_end_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_beg_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_end_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_beg_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()
		fileOut=open("%s/emmit_end_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()	
		fileOut=open("%s/emmit_end_DB_WASTE_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'w')
		fileOut.close()	
		myEmitNode_beg=Calc_Emit("MyEmitNode_Beg_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_beg_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_begQ=Calc_Emit("MyEmitNode_BegQ_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_begQ_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_beg_DH11=Calc_Emit("myEmitNode_beg_DH11_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_beg_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_end_DH11=Calc_Emit("myEmitNode_end_DH11_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_end_DH11_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_beg_b23=Calc_Emit("myEmitNode_beg_b23_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_beg_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_mid_b23=Calc_Emit("myEmitNode_mid_b23_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_mid_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_end_b23=Calc_Emit("myEmitNode_end_b23_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_end_b23_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_beg_DH12=Calc_Emit("myEmitNode_beg_DH12_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_beg_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_end_DH12=Calc_Emit("myEmitNode_end_DH12_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_end_DH12_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_beg_DH13=Calc_Emit("myEmitNode_beg_DH13_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_beg_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_end_DH13=Calc_Emit("myEmitNode_end_DH13_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_end_DH13_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))
		myEmitNode_end_DB_WASTE=Calc_Emit("myEmitNode_end_DB_WASTE_%d_%d_%d_%d"%(currentPart,currentPart2,nPartsChicane,nPartsChicane2),True,"%s/emmit_end_DB_WASTE_%d_%d_%d_%d.txt"%(outputDirectory,currentPart,currentPart2,nPartsChicane,nPartsChicane2))		
		
		#place the emmit and print monitoring child nodes into the appropriate places of the lattice
		nodes = inj_latt_start.getNodes()
		if usePrintNode:
			nodes[0].addChildNode(myPrintNode_begQ,AccNode.ENTRANCE)
		nodes[0].addChildNode(myEmitNode_begQ,AccNode.ENTRANCE)
		i = 0
		path_length=0
		for node in nodes:
			pass
			if node.getName().strip() == "DB12":
				if usePrintNode:
					node.addChildNode(myPrintNode_beg,AccNode.ENTRANCE)
				node.addChildNode(myEmitNode_beg,AccNode.ENTRANCE)	
			if node.getName().strip() == "DB23":
				if usePrintNode:
					node.addChildNode(myPrintNode_beg_b23,AccNode.ENTRANCE)
					node.addChildNode(myPrintNode_end_b23,AccNode.EXIT)		
				
				node.addChildNode(myEmitNode_beg_b23,AccNode.ENTRANCE)
				node.addChildNode(myEmitNode_end_b23,AccNode.EXIT)
				
				
			if node.getName().strip() == "DH_A11":
				if usePrintNode:
					node.addChildNode(myPrintNode_beg_DH11,AccNode.ENTRANCE)
					node.addChildNode(myPrintNode_end_DH11,AccNode.EXIT)		
				
				node.addChildNode(myEmitNode_beg_DH11,AccNode.ENTRANCE)
				node.addChildNode(myEmitNode_end_DH11,AccNode.EXIT)
				
				if currentPart is not nPartsChicane:
					if usePrintNode:
						node.addChildNode(myPrintNode_postS_DH11,AccNode.BODY,currentPart)
			if node.getName().strip() == "DH_A12":
				if usePrintNode:
					node.addChildNode(myPrintNode_beg_DH12,AccNode.ENTRANCE)
					node.addChildNode(myPrintNode_end_DH12,AccNode.EXIT)		
				
				node.addChildNode(myEmitNode_beg_DH12,AccNode.ENTRANCE)
				node.addChildNode(myEmitNode_end_DH12,AccNode.EXIT)	
	
			if node.getName().strip() == "DH_A13":
				if usePrintNode:
					node.addChildNode(myPrintNode_beg_DH13,AccNode.ENTRANCE)
					node.addChildNode(myPrintNode_end_DH13,AccNode.EXIT)		
				
				node.addChildNode(myEmitNode_beg_DH13,AccNode.ENTRANCE)
				node.addChildNode(myEmitNode_end_DH13,AccNode.EXIT)		
				
			if node.getName().strip() == "DB_Waste":
				if usePrintNode:
					node.addChildNode(myPrintNode_end_DB_WASTE,AccNode.EXIT)		
				node.addChildNode(myEmitNode_end_DB_WASTE,AccNode.EXIT)			
		#================Do some turns===========================================
		
	
		
		#track through drift after 3rd chicane and foil at end
		inj_latt_start.trackBunch(bunch_in, paramsDict)
	
		
		#===========Dump bunch infomration=======================================
		#bunch_pyorbit_to_orbit(inj_latt.getLength(), bunch_in, "mainbunch.dat")
		#bunch_pyorbit_to_orbit(inj_latt.getLength(), lostbunch, "lostbunch.dat")
		print "Stop."
	
	

