#############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################

import random
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
parser.add_argument("--outputDirectory", dest='outputDirectory', default="Test_InjectBeam3_General", help="Where to put output")
parser.add_argument("--inputDirectoryPrefix", dest='inputDirectoryPrefix', default="Test_InjectBeam3_General", help="Where to put output")
parser.add_argument("--inputDirectorySuffix", dest='inputDirectorySuffix', default="Test_InjectBeam3_General", help="Where to put output")
parser.add_argument("--configFile", dest='configFile', default="BeamLatticeFiles/InjectBeam.txt", help="configure file")


args = parser.parse_args()

beamLatticeDictionary=ConfigureFileReader(args.configFile)
beamLatticeDictionary.printDictionary()


nParticles = 1
if beamLatticeDictionary.hasKey("nParts"):
	nParticles=int(beamLatticeDictionary.getValue("nParts"))


nStates=1
if beamLatticeDictionary.hasKey("nStates"):
	nStates=int(beamLatticeDictionary.getValue("nStates"))
else:
	print "nStates not in beamLattice config file, exiting"
	sys.exit(0)
outputDirectory=args.outputDirectory

if not os.path.exists(outputDirectory):
	os.mkdir(outputDirectory)
	
usePrintNode=False
if beamLatticeDictionary.hasKey("usePrintNode") and beamLatticeDictionary.getValue("usePrintNode")=="True":
	usePrintNode=True
	
configFileNameArray=args.configFile.split("/")
configFileName=configFileNameArray[len(configFileNameArray)-1]
print configFileName
outputDirectory=outputDirectory+"/%s"%(configFileName.split(".")[0])
print outputDirectory
if not os.path.exists(outputDirectory):
	os.mkdir(outputDirectory)
#loop over all configurations of where 1st stripper dipole with be placed
#currentPart=0 places it at begginging of chicane2/11
#currentPart=nPartsChicane places it immediately after chicane2/11
#0<currentPart<nPartsChicane places it currentPart/nPartsChicane fractionally into chicane2/11
#for currentPart in range(-1,nPartsChicane+1):

startOfNamesArray=["beg_","end_"]
nodeList=beamLatticeDictionary.getArray("nodeNames")
for currentNode in nodeList:
	for startOfName in startOfNamesArray:
		bunch_in=Bunch()
		
		for index in range(nStates):
			if beamLatticeDictionary.hasKey("state%d"%(index+1)):
				stateName=beamLatticeDictionary.getValue("state%d"%(index+1))
				stateWeight=float(beamLatticeDictionary.getValue("weight%d"%(index+1)))
				
				openedFile=open("%s%s_%s/print_%s%s.txt"%(args.inputDirectoryPrefix,stateName,args.inputDirectorySuffix,startOfName,currentNode),'r')
				lines=openedFile.readlines()
				counter=0
				for line in lines:
					coordToken=line.split("=")[1]
					coordArr=coordToken.split(",")
					if (random.random()<stateWeight):
						bunch_in.addParticle(float(coordArr[0].strip('()\n ')),float(coordArr[1].strip('()\n ')),float(coordArr[2].strip('()\n ')),float(coordArr[3].strip('()\n ')),float(coordArr[4].strip('()\n ')),float(coordArr[5].strip('()\n ')))
					counter=counter+1
					if counter>=nParticles:
						break
			else:
				print "state%d not in beamLattice config file, exiting"%(index+1)
				sys.exit(0)			

		if usePrintNode:
			myPrintNodeBeg=None
			
			fileOut=open("%s/print_%s%s.txt"%(outputDirectory,startOfName,currentNode),'w')
			fileOut.close()		
			myPrintNodeBeg=Print_Node("MyPrintNode_%s%s"%(startOfName,currentNode),True,"%s/print_%s%s.txt"%(outputDirectory,startOfName,currentNode))

		myEmitNodeBeg=None

		fileOut=open("%s/emmit_%s%s.txt"%(outputDirectory,startOfName,currentNode),'w')
		fileOut.close()		
		myEmitNodeBeg=Calc_Emit("MyEmitNode_%s%s"%(startOfName,currentNode),True,"%s/emmit_%s%s.txt"%(outputDirectory,startOfName,currentNode))

		paramsDict = {}
		paramsDict["bunch"]= bunch_in	
		if usePrintNode:
			myPrintNodeBeg.track(paramsDict)
		myEmitNodeBeg.track(paramsDict)
		
print "Stop."
