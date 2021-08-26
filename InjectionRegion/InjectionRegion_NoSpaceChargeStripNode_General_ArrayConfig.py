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
#from KevinPython.function_stripping import probabilityStripping
#from KevinPython.function_strippingIncludeChicaneField import probabilityStrippingWithChicane
from orbit.teapot import addDipoleStripperNode

import argparse

from ConfigureFileClass import ConfigureFileReader
from MagneticFieldClass import MagneticField
from KevinPython.changeCharge import Change_Charge_Child

#finds index of first instance of node with name
#if it isnt found returns -1
def findReferenceNode(theLattice,nameToFind):
	counter=0
	for node in theLattice.getNodes():
		if node.getName().strip() == nameToFind:
			return counter			
		counter+=1	
	return -1
print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--printNodes",type=bool, dest='printNodes', default=True, help="print node list")
parser.add_argument("--addChicaneFieldToStripper",type=bool, dest='addChicaneFieldToStripper', default=True, help="Include the chicane fields in the stripper if stripper is inside chicane")

parser.add_argument("--outputDirectory", dest='outputDirectory', default="Test_InjectBeam3_General", help="Where to put output")
parser.add_argument("--magneticFieldFile", dest='magneticFieldFile', default="MagneticFieldFiles/magneticFieldUpUp.txt", help="infoOnMagneticField")
parser.add_argument("--beamLatticeFile", dest='beamLatticeFile', default="BeamLatticeFiles/InjectBeam.txt", help="infoOnInitalBeamAndLattices")
parser.add_argument("--chicaneScaleDirectory", dest='chicaneScaleDirectory', default="InjectBeam3_ChangeOffset", help="Where to get chicane scales from")

args = parser.parse_args()

beamLatticeDictionary=ConfigureFileReader(args.beamLatticeFile)
beamLatticeDictionary.printDictionary()

startingNumber=8
#the fudgeFactor is the spacing used when inserting the stripping dipole immediately before or after a node
#without it rounding errors can place its end or begginning extending past the edge of reference node. So it is shifted the fudgeFactor from the edge
fudgeFactor=.00000001
#=====Main bunch parameters============
intensity = 7.8e13
turns = 1
#macrosperturn = 260
macrosperturn = 1
if beamLatticeDictionary.hasKey("nParts"):
	macrosperturn=int(beamLatticeDictionary.getValue("nParts"))
macrosize = intensity/turns/macrosperturn

#where to pull chicane scale values from if useChicaneScales is true.
#outputDirectoryChicaneScales="WasteBeamSplitGeneralNewStripperChicaneFieldAddedClean"
outputDirectoryChicaneScales=args.chicaneScaleDirectory

chicaneScale10=-1.
chicaneScale11=-1.
chicaneScale12=-1.
chicaneScale13=-1.

#the default chicane kick strength array
chicaneStrengthArray=[-0.041456,0.052434,0.0298523,-0.0398609]
#change them if doing soemthing different IE PPU
if beamLatticeDictionary.hasKey("chicaneKickStrength"):
	chicaneKickArray=beamLatticeDictionary.getArray("chicaneKickStrength")
	chicaneStrengthArray[0]=-float(chicaneKickArray[0])
	chicaneStrengthArray[1]=-float(chicaneKickArray[1])
	chicaneStrengthArray[2]=-float(chicaneKickArray[2])
	chicaneStrengthArray[3]=-float(chicaneKickArray[3])
usePrintNode=False
if beamLatticeDictionary.hasKey("usePrintNode") and beamLatticeDictionary.getValue("usePrintNode")=="True":
	usePrintNode=True

doDipoleStrippers=False
#print "beamLatticeDictionary.getValue(\"doDipoleStrippers\")",beamLatticeDictionary.getValue("doDipoleStrippers")
#print "hi=%s"%beamLatticeDictionary.getValue("doDipoleStrippers").strip()
#print beamLatticeDictionary.getValue("doDipoleStrippers")
if beamLatticeDictionary.hasKey("doDipoleStrippers") and beamLatticeDictionary.getValue("doDipoleStrippers")=="True":
	#print "lets go"
	doDipoleStrippers=True
numberOfStripperDipoles=-1
if doDipoleStrippers:
	if beamLatticeDictionary.hasKey("numberOfStripperDipoles"):
		numberOfStripperDipoles=int(beamLatticeDictionary.getValue("numberOfStripperDipoles"))
	else:
		print "numberOfStripperDipoles not in beamLattice config file, exiting"
		sys.exit(0)
outputDirectory=args.outputDirectory
if not os.path.exists(outputDirectory):
	os.mkdir(outputDirectory)
	
useSecondaryFoil=False
latticeIndexToAddFoilTo=-1
if beamLatticeDictionary.hasKey("useFoil") and beamLatticeDictionary.getValue("useFoil")=="True":
	useSecondaryFoil=True
	#adds foil to end of this lattice
	latticeIndexToAddFoilTo=int(beamLatticeDictionary.getValue("latticeToAddFoilTo"))
	
useChicaneScaleFile=False
chicaneScaleFile_UseScales=False
chicaneScaleFile_UsePX=False
chicaneScaleFile_UsePY=False
chicaneScaleFile_UsePX_Position=4
chicaneScaleFile_UsePY_Position=5

if beamLatticeDictionary.hasKey("useChicaneScaleFile") and beamLatticeDictionary.getValue("useChicaneScaleFile")=="True":
	useChicaneScaleFile=True
	if beamLatticeDictionary.hasKey("chicaneScaleFile_UseScales") and beamLatticeDictionary.getValue("chicaneScaleFile_UseScales")=="True":
		chicaneScaleFile_UseScales=True
	if beamLatticeDictionary.hasKey("chicaneScaleFile_UsePX") and beamLatticeDictionary.getValue("chicaneScaleFile_UsePX")=="True":
		chicaneScaleFile_UsePX=True
	if beamLatticeDictionary.hasKey("chicaneScaleFile_UsePY") and beamLatticeDictionary.getValue("chicaneScaleFile_UsePY")=="True":
		chicaneScaleFile_UsePY=True	
	if beamLatticeDictionary.hasKey("chicaneScaleFile_UsePX_Position"):
		chicaneScaleFile_UsePX_Position=int(beamLatticeDictionary.getValue("chicaneScaleFile_UsePX_Position"))
	if beamLatticeDictionary.hasKey("chicaneScaleFile_UsePY_Position"):
		chicaneScaleFile_UsePY_Position=int(beamLatticeDictionary.getValue("chicaneScaleFile_UsePY_Position"))			
#------------------------------
#Initial Distribution Functions
#------------------------------

e_kin_ini = float(beamLatticeDictionary.getValue("e_kin_ini"))
mass =  float(beamLatticeDictionary.getValue("mass"))
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
c=299792458
momentum=gamma*beta*mass
rigidity=momentum/(c/math.pow(10.,9))
print "relat. gamma=",gamma
print "relat.  beta=",beta
print "relat.  mom=",momentum
print "rigidity=",rigidity


#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
(alphaZ,betaZ,emittZ) = ( float(beamLatticeDictionary.getValue("alphaZ")), float(beamLatticeDictionary.getValue("alphaZ")), float(beamLatticeDictionary.getValue("alphaZ")))

(alphaX,betaX,emittX) = ( float(beamLatticeDictionary.getValue("alphaX")), float(beamLatticeDictionary.getValue("betaX")), float(beamLatticeDictionary.getValue("emittX")))
(alphaY,betaY,emittY) = ( float(beamLatticeDictionary.getValue("alphaY")), float(beamLatticeDictionary.getValue("betaY")), float(beamLatticeDictionary.getValue("emittY")))

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
#loop over all configurations of where 1st stripper dipole with be placed
#currentPart=0 places it at begginging of chicane2/11
#currentPart=nPartsChicane places it immediately after chicane2/11
#0<currentPart<nPartsChicane places it currentPart/nPartsChicane fractionally into chicane2/11
#for currentPart in range(-1,nPartsChicane+1):

#this is a temporary fix to not use chicane scale files ending in NA_NA_NA_NA
useChicane_NA_NA_NA_NA_File=False
theLatticesList=[]
#for now strippers must be in first lattice
latticeFileNameList=beamLatticeDictionary.getArray("lattice")
for currentLatticeName in latticeFileNameList:
	currentLatticeIndex=len(theLatticesList)
	inj_latt_start = teapot.TEAPOT_Ring()
	print "Read MAD."
	#this lattice has the injection region from the start of the drift prior to chicane2 up to and including the drift after chicane3
	inj_latt_start.readMAD(currentLatticeName,"RING")
	print "Lattice=",inj_latt_start.getName()," length [m] =",inj_latt_start.getLength()," nodes=",len(inj_latt_start.getNodes())
	
	#turn off kickers
	kickerwave = flatTopWaveform(1.0)
	strength_hkicker10 = 0
	strength_hkicker13 = strength_hkicker10
	strength_hkicker11 = 0
	strength_hkicker12 = strength_hkicker11
	strength_vkicker10 = 0
	strength_vkicker13 = strength_vkicker10
	strength_vkicker11 = 0
	strength_vkicker12 = strength_vkicker11		
	nodes = inj_latt_start.getNodes()
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
		print "turning off a kicker"
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
	if useSecondaryFoil and currentLatticeIndex==latticeIndexToAddFoilTo:
		#create secondary foil
		thick = 400.0
		foil = TeapotFoilNode(-100, 100, -100, 100, thick, "Foil 2")
		scatterchoice = 0
		foil.setScatterChoice(scatterchoice)
		addTeapotFoilNode(inj_latt_start,inj_latt_start.getLength(),foil)	
	
	#set the strength of the chicane kicks
	if useChicaneScaleFile:
		if not doDipoleStrippers and useChicane_NA_NA_NA_NA_File:
			openedFile=open("%s/ChicaneScales_%s_%s_%s_%s.txt"%(outputDirectoryChicaneScales,"NA","NA","NA","NA"),'r')
		else:
			openedFile=open("%s/ChicaneScales.txt"%(outputDirectoryChicaneScales),'r')
		line=openedFile.readline()
		print line
		theScales=line.split(",")
		if chicaneScaleFile_UseScales:
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
	nodes = inj_latt_start.getNodes()
	for index in range(len(inj_latt_start.getNodes())):
		if nodes[index].getName().strip() == "DH_A10":
			chicane10=nodes[index]
		elif nodes[index].getName().strip() == "DH_A11":
			chicane11=nodes[index]
		elif nodes[index].getName().strip() == "DH_A12":
			chicane12=nodes[index]
		elif nodes[index].getName().strip() == "DH_A13":
			chicane13=nodes[index]	
		
	chicanewave = flatTopWaveform(1.0)
	if chicane10 is not None:
		chicane10.setParam("kx", strength_chicane10)
		chicane10.setWaveform(chicanewave)
	if chicane11 is not None:
		chicane11.setParam("kx", strength_chicane11)
		chicane11.setWaveform(chicanewave)
	if chicane12 is not None:
		chicane12.setParam("kx", strength_chicane12)
		chicane12.setWaveform(chicanewave)
	if chicane13 is not None:
		chicane13.setParam("kx", strength_chicane13)
		chicane13.setWaveform(chicanewave)				
	magneticFields=[]
	if doDipoleStrippers:
		for index in range(numberOfStripperDipoles):
			scale=1
			offset=0
			if beamLatticeDictionary.hasKey("stripper%d"%(index+1)):
				tempFieldClass=None
				#stripperX=["<magnetic field file>","<Dipole Name>",isStripper<True/False>,<Name of node to place it in reference too>,<position with respect to reference node>]
				#<position with respect to reference node> can be "before"(places it right before reference), "after"(places it right after reference), or float [0-1](places it this fraction into reference)				
				tempFieldArray=beamLatticeDictionary.getArray("stripper%d"%(index+1))
				tempFieldClass=MagneticField(tempFieldArray[0],tempFieldArray[1],tempFieldArray[2],tempFieldArray[3],tempFieldArray[4],tempFieldArray[5])
				magneticFields.append(tempFieldClass)
			else:
				print "stripper%d not in beamLattice config file, exiting"%(index+1)
				sys.exit(0)	
			
			if useChicaneScaleFile and beamLatticeDictionary.hasKey("chicaneScaleFile_useStripperLength%d"%(index+1)) and beamLatticeDictionary.getValue("chicaneScaleFile_useStripperLength%d"%(index+1))=="True":
				openedFile=open("%s/ChicaneScales.txt"%(outputDirectoryChicaneScales),'r')
				line=openedFile.readline()
				#print line
				theScales=line.split(",")
				scale=float(theScales[startingNumber+index*2].strip())			
				magneticFields[index].setStripperLength(magneticFields[index].getStripperLength()*scale)
			if useChicaneScaleFile and beamLatticeDictionary.hasKey("chicaneScaleFile_useStripperAngle%d"%(index+1)) and beamLatticeDictionary.getValue("chicaneScaleFile_useStripperAngle%d"%(index+1))=="True":
				openedFile=open("%s/ChicaneScales.txt"%(outputDirectoryChicaneScales),'r')
				line=openedFile.readline()
				#print line
				theScales=line.split(",")
				offset=float(theScales[startingNumber+1+index*2].strip())			
				magneticFields[index].setFieldDirection(magneticFields[index].getFieldDirection()+offset)	
			if useChicaneScaleFile and beamLatticeDictionary.hasKey("chicaneScaleFile_useStripperOffset%d"%(index+1)) and beamLatticeDictionary.getValue("chicaneScaleFile_useStripperOffset%d"%(index+1))=="True":
				openedFile=open("%s/ChicaneScales.txt"%(outputDirectoryChicaneScales),'r')
				line=openedFile.readline()
				#print line
				theScales=line.split(",")
				offset=float(theScales[startingNumber+numberOfStripperDipoles*2+index].strip())			
				magneticFields[index].setNodePositionOffset(magneticFields[index].getNodePositionOffset()+offset)					
				
		for index in range(numberOfStripperDipoles):
			refNodeIndex=findReferenceNode(inj_latt_start,magneticFields[index].getRefNodeName().strip())
			if refNodeIndex>=0:
				position=-100.
				if magneticFields[index].getNodePosition().lower()=="before":
					position =inj_latt_start.getNodePositionsDict()[inj_latt_start.getNodes()[refNodeIndex]][0]-float(magneticFields[index].getStripperLength())-fudgeFactor
				elif magneticFields[index].getNodePosition().lower()=="after":
					position =inj_latt_start.getNodePositionsDict()[inj_latt_start.getNodes()[refNodeIndex]][1]+fudgeFactor
				else:
					position =inj_latt_start.getNodePositionsDict()[inj_latt_start.getNodes()[refNodeIndex]][0]+inj_latt_start.getNodes()[refNodeIndex].getLength()*float(magneticFields[index].getNodePosition())
					
				position=position+magneticFields[index].getNodePositionOffset()
				#check if we are in kicker or drift
				position_start = position
				position_stop = position + float(magneticFields[index].getStripperLength())
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
					print "something is going to be broken1"
					print position_start
					print position_stop
					sys.exit(0)
				nodeCurrent=inj_latt_start.getNodes()[node_start_ind]
				xkickerField=0.
				ykickerField=0.
				#nothing to change because no field in drift
				if(isinstance(nodeCurrent,DriftTEAPOT)):
					magneticFields[index].setIsInsideChicane(False)
					#print "stripper dipole is in drift"
				#in kicker so add kicker field to stripping field
				elif (isinstance(nodeCurrent,KickTEAPOT)):
					#print "stripper dipole is in kick node"
					if args.addChicaneFieldToStripper:
						magneticFields[index].setIsInsideChicane(True)
						if nodeCurrent.getName().strip() == "DH_A10":
							magneticFields[index].setChicaneItIsInside(0)
						elif nodeCurrent.getName().strip() == "DH_A11":
							magneticFields[index].setChicaneItIsInside(1)
						elif nodeCurrent.getName().strip() == "DH_A12":
							magneticFields[index].setChicaneItIsInside(2)
						elif nodeCurrent.getName().strip() == "DH_A13":
							magneticFields[index].setChicaneItIsInside(3)							
						#compute field to create kick
						length=nodeCurrent.getLength()
						kx=nodeCurrent.getParam("kx")
						ykickerField=-kx*rigidity/length
						ky=nodeCurrent.getParam("ky")
						xkickerField=ky*rigidity/length
				#actually creates the field functions
				nParts=int(magneticFields[index].getNParts())
				lengthStripper=float(magneticFields[index].getStripperLength())
				stepSize=lengthStripper/nParts
				fieldDirection=float(magneticFields[index].getFieldDirection())
				magneticFieldx= Function()
				magneticFieldy= Function()
				for i in range(nParts):
					x = stepSize*i;
					#y = constantField(x)
					y = float(magneticFields[index].getValueOfField(x))
					magneticFieldx.add(x,y*math.cos(fieldDirection)+xkickerField)
					magneticFieldy.add(x,y*math.sin(fieldDirection)+ykickerField)
				
				if magneticFields[index].getIsStripper()=="True":
					myDipole_DH_A11=GeneralDipoleStripSeperateField(magneticFieldx,magneticFieldy,nParts,lengthStripper,gamma,beta,magneticFields[index].getNodeName(),float(magneticFields[index].getFixedStrippingLength()))
				else:
					myDipole_DH_A11=GeneralDipoleNoStripSeperateField(magneticFieldx,magneticFieldy,nParts,lengthStripper,gamma,beta,magneticFields[index].getNodeName())
					if magneticFields[index].getFoilTest():
						chargeChangeNode=Change_Charge_Child("foil1",1)
						myDipole_DH_A11.addChildNode(chargeChangeNode,AccNode.EXIT)						
				#print "xkickerField=",xkickerField
				myDipole_DH_A11.setChicaneFieldx(xkickerField)
				myDipole_DH_A11.setChicaneFieldy(ykickerField)
				addDipoleStripperNode(inj_latt_start,position,myDipole_DH_A11)					

	i = 0
	path_length=0
	for node in nodes:
		if args.printNodes==True:
			path_length_center=path_length+node.getLength()/2.
			path_length=path_length+node.getLength()
			print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%inj_latt_start.getNodePositionsDict()[node], " path_length= ",path_length, " path_length_center= ",path_length_center
			
			i=i+1	
		if usePrintNode:
			myPrintNodeBeg=None
			myPrintNodeEnd=None
			if not doDipoleStrippers and useChicane_NA_NA_NA_NA_File:
				fileOut=open("%s/print_beg_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()),'w')
				fileOut.close()		
				fileOut=open("%s/print_end_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()),'w')
				fileOut.close()		
				myPrintNodeBeg=Print_Node("MyPrintNode_beg_%s_NA_NA_NA_NA"%(node.getName()),True,"%s/print_beg_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()))
				myPrintNodeEnd=Print_Node("MyPrintNode_end_%s_NA_NA_NA_NA"%(node.getName()),True,"%s/print_end_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()))
			else:
				fileOut=open("%s/print_beg_%s.txt"%(outputDirectory,node.getName()),'w')
				fileOut.close()		
				fileOut=open("%s/print_end_%s.txt"%(outputDirectory,node.getName()),'w')
				fileOut.close()	
				myPrintNodeBeg=Print_Node("MyPrintNode_beg_%s"%(node.getName()),True,"%s/print_beg_%s.txt"%(outputDirectory,node.getName()))
				myPrintNodeEnd=Print_Node("MyPrintNode_end_%s"%(node.getName()),True,"%s/print_end_%s.txt"%(outputDirectory,node.getName()))

			node.addChildNode(myPrintNodeBeg,AccNode.ENTRANCE)
			node.addChildNode(myPrintNodeEnd,AccNode.EXIT)
		myEmitNodeBeg=None
		myEmitNodeEnd=None
		if not doDipoleStrippers and useChicane_NA_NA_NA_NA_File:
			fileOut=open("%s/emmit_beg_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()),'w')
			fileOut.close()		
			fileOut=open("%s/emmit_end_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()),'w')
			fileOut.close()		
			myEmitNodeBeg=Calc_Emit("MyEmitNode_beg_%s_NA_NA_NA_NA"%(node.getName()),True,"%s/emmit_beg_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()))
			myEmitNodeEnd=Calc_Emit("MyEmitNode_end_%s_NA_NA_NA_NA"%(node.getName()),True,"%s/emmit_end_%s_NA_NA_NA_NA.txt"%(outputDirectory,node.getName()))					
		else:
			fileOut=open("%s/emmit_beg_%s.txt"%(outputDirectory,node.getName()),'w')
			fileOut.close()		
			fileOut=open("%s/emmit_end_%s.txt"%(outputDirectory,node.getName()),'w')
			fileOut.close()	
			myEmitNodeBeg=Calc_Emit("MyEmitNode_beg_%s"%(node.getName()),True,"%s/emmit_beg_%s.txt"%(outputDirectory,node.getName()))
			myEmitNodeEnd=Calc_Emit("MyEmitNode_end_%s"%(node.getName()),True,"%s/emmit_end_%s.txt"%(outputDirectory,node.getName()))
		node.addChildNode(myEmitNodeBeg,AccNode.ENTRANCE)
		node.addChildNode(myEmitNodeEnd,AccNode.EXIT)					
				
	
	#------------------------------
	#  Lattice is ready
	#-------------------------------
	theLatticesList.append(inj_latt_start)

		
			

		
#create the initial Bunch

bunch_in=Bunch()
bunchFromFile=False
#use pencil beam
if beamLatticeDictionary.hasKey("pencilBeam") and beamLatticeDictionary.getValue("pencilBeam")=="True":
	bunch_in.addParticle(0,0,0,0,0,0)
#read initial bunch in from file
elif beamLatticeDictionary.hasKey("bunchFromFile") and beamLatticeDictionary.getValue("bunchFromFile")=="True":
	bunchFromFile=True
	openedFile=open("%s"%(beamLatticeDictionary.getValue("bunchFromFileName")),'r')
	lines=openedFile.readlines()
	for line in lines:
		coordToken=line.split("=")[1]
		coordArr=coordToken.split(",")
		bunch_in.addParticle(float(coordArr[0].strip('()\n ')),float(coordArr[1].strip('()\n ')),float(coordArr[2].strip('()\n ')),float(coordArr[3].strip('()\n ')),float(coordArr[4].strip('()\n ')),float(coordArr[5].strip('()\n ')))
#generate initial bunch
else:
	bunch_in = bunch_gen.getBunch(nParticles = macrosperturn, distributorClass = WaterBagDist3D)
	#bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
	#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)


bunch_in.mass(mass) #mass
bunch_in.macroSize(macrosize)
energy = e_kin_ini # 1.0 #Gev
bunch_in.getSyncParticle().kinEnergy(energy)
initial_charge=-1
if beamLatticeDictionary.hasKey("initial_charge"):
	initial_charge=int(beamLatticeDictionary.getValue("initial_charge"))
else:
	print "initial_charge not set using initial_charge=",initial_charge
bunch_in.charge(initial_charge)

#apply offset to initial beam
xOffset=0.25671
pxOffset=-.042
yOffset=0.046
pyOffset=0

#use default values if not set in file
if beamLatticeDictionary.hasKey("xOffset"):
	xOffset=float(beamLatticeDictionary.getValue("xOffset"))
if beamLatticeDictionary.hasKey("pxOffset"):
	pxOffset=float(beamLatticeDictionary.getValue("pxOffset"))
if beamLatticeDictionary.hasKey("yOffset"):
	yOffset=float(beamLatticeDictionary.getValue("yOffset"))
if beamLatticeDictionary.hasKey("pyOffset"):
	pyOffset=float(beamLatticeDictionary.getValue("pyOffset"))
	
#overide offsets with values from chicane file if set to use them
if useChicaneScaleFile and (chicaneScaleFile_UsePX or chicaneScaleFile_UsePY):
	if not doDipoleStrippers and useChicane_NA_NA_NA_NA_File:
		openedFile=open("%s/ChicaneScales_NA_NA_NA_NA.txt"%(outputDirectoryChicaneScales),'r')
	else:
		openedFile=open("%s/ChicaneScales_%d_%d_%d_%d.txt"%(outputDirectoryChicaneScales,currentPart,currentPart2,nPartsChicane,nPartsChicane2),'r')
	line=openedFile.readline()
	print line
	theScales=line.split(",")
	if chicaneScaleFile_UsePX:
		pxOffset=float(theScales[chicaneScaleFile_UsePX_Position].strip())
	if chicaneScaleFile_UsePY:
		pyOffset=float(theScales[chicaneScaleFile_UsePY_Position].strip())				
	openedFile.close()
#if reading bunch from file the offset should already have been added
if not bunchFromFile:
	for i in range(bunch_in.getSize()):
		bunch_in.x(i,bunch_in.x(i)+xOffset)
		bunch_in.px(i,bunch_in.px(i)+pxOffset)
		bunch_in.y(i,bunch_in.y(i)+yOffset)
		bunch_in.py(i,bunch_in.py(i)+pyOffset)	
		
paramsDict = {}
#create failed to strip bunches
#these are only used if stripper is actually a stripper
firstChicaneFail = Bunch()
firstChicaneFail.charge(-1)
secondChicaneFail = Bunch()
secondChicaneFail.charge(0)
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch_in
paramsDict["firstChicaneFail"]=firstChicaneFail
paramsDict["secondChicaneFail"]= secondChicaneFail
lostbunch.addPartAttr("LostParticleAttributes") 

latticeIndexToMakeBunchChargePlusOneAtEndOf=-1
if beamLatticeDictionary.hasKey("latticeIndexToMakeBunchChargePlusOneAtEndOf"):
	latticeIndexToMakeBunchChargePlusOneAtEndOf=int(beamLatticeDictionary.getValue("latticeIndexToMakeBunchChargePlusOneAtEndOf"))
#track through all the lattices
for currentLatticeIndex in range(len(theLatticesList)):
	theLatticesList[currentLatticeIndex].trackBunch(bunch_in, paramsDict)
	#set bunch charge to +1 at end of lattice. More desired approach would be to create child node that does this.
	if currentLatticeIndex==latticeIndexToMakeBunchChargePlusOneAtEndOf:
		bunch_in.charge(1)

#===========Dump bunch infomration=======================================
#bunch_pyorbit_to_orbit(inj_latt.getLength(), bunch_in, "mainbunch.dat")
#bunch_pyorbit_to_orbit(inj_latt.getLength(), lostbunch, "lostbunch.dat")
print "Stop."
