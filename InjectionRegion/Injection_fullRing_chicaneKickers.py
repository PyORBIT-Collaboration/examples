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

from orbit.teapot import GeneralDipole
from orbit.teapot import YDipole
from orbit.teapot import XDipole
from KevinPython.printNode import Print_Node
from KevinPython.calculateEmit import Calc_Emit
import argparse

print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--outputDirectory", dest='outputDirectory', default="OG_Injection_Kicker", help="Where to put output")
parser.add_argument("--usePrintNode",type=bool, dest='usePrintNode', default=False, help="whether or not to print bunch info to file")
parser.add_argument("--nParts",type=int, dest='nParts', default=10000, help="number of particles")
parser.add_argument("--turns",type=int, dest='turns', default=1, help="number of complete orbits")
parser.add_argument("--printNodes",type=bool, dest='printNodes', default=True, help="print node list")
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=False, help="print node list")
parser.add_argument("--doNormalKickers",type=bool, dest='doNormalKickers', default=False, help="print node list")
parser.add_argument("--chicaneScaleDirectory", dest='chicaneScaleDirectory', default="Method2_strippersNotClosed_ArrayConfig_FoilTest/Method2_Part2_strippersNotClosed_FloatLength_FoilTest", help="Where to get chicane scales from")
args = parser.parse_args()

useChicaneScaleFile=True
outputDirectoryChicaneScales=args.chicaneScaleDirectory

chicaneScale10=-1.
chicaneScale11=-1.
chicaneScale12=-1.
chicaneScale13=-1.

#the default chicane kick strength array
chicaneStrengthArray=[-0.041456,0.052434,0.0298523,-0.0398609]

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
print "Read MAD."
teapot_latt.readMAD("MAD_Injection_Region_Lattice/OG_Injection_InjectionRegionOnly_Chicane_Replaced_With_Kickers_onlyChicane2.LAT","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())
nodes2 = teapot_latt.getNodes()
#numberOfCustomDipoles=2

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


lattlength = teapot_latt.getLength()
sp = b.getSyncParticle()
kickerwave = rootTWaveform(sp, lattlength, duration, startamp, endamp)
chicanewave = flatTopWaveform(1.0)

nodes = teapot_latt.getNodes()
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


#set the strength of the chicane kicks
if useChicaneScaleFile:
	
	openedFile=open("%s/ChicaneScales.txt"%(outputDirectoryChicaneScales),'r')
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
nodes = teapot_latt.getNodes()
for index in range(len(teapot_latt.getNodes())):
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
#chicane10.setWaveform(chicanewave)
#chicane11.setWaveform(chicanewave)
#chicane12.setWaveform(chicanewave)
#chicane13.setWaveform(chicanewave)

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
xcenterpos = 0.0468
#xcenterpos = 0.0
xcentermom = 0.00
ycenterpos = 0.0492
#ycenterpos = 0.0
ycentermom = 0.00

zlim = 120. * lattlength/360.
zmin = -zlim
zmax = zlim
tailfraction = 0
emean = sp.kinEnergy()
efac = 0.784
esigma = 0.0015*efac
etrunc = 1.
emin = sp.kinEnergy() - 0.0025*efac
emax = sp.kinEnergy() + 0.0025*efac
ecmean = 0
ecsigma = 0.0015*efac
ectrunc = 1.
ecmin = -0.0035*efac
ecmax = 0.0035*efac
ecdrifti = 0
ecdriftf = 0
turns = 1000.
tturn = lattlength / (sp.beta() * 2.998e8)
drifttime= 1000.*turns*tturn
ecparams = (ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime)

esnu = 100.
esphase = 0.
esmax = 0
nulltime = 0
esparams = (esnu, esphase, esmax, nulltime) 

sp = b.getSyncParticle()

xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom)
xFunc.setSeed(1)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom)
yFunc.setSeed(1)
lFunc = SNSESpreadDist(lattlength, zmin, zmax, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams,1)

#====Injection and foil aperature============

xmin = xcenterpos - 0.0085
xmax = xcenterpos + 0.0085
ymin = ycenterpos - 0.0080
ymax = ycenterpos + 0.100

#print xmin
#print xmax
#print ymin
#print ymax
#=================Add the injection node and foil node==  ==============

nparts = macrosperturn
injectparams = (xmin, xmax, ymin, ymax)
injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
addTeapotInjectionNode(teapot_latt, 0., injectnode) 

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
outputDirectory=args.outputDirectory
i = 0
path_length=0
for node in nodes:
	if args.printNodes==True:
		path_length=path_length+node.getLength()
		print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node], " path_length= ",path_length
		i=i+1	
	if args.usePrintNode:
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

	fileOut=open("%s/emmit_beg_%s.txt"%(outputDirectory,node.getName()),'w')
	fileOut.close()		
	fileOut=open("%s/emmit_end_%s.txt"%(outputDirectory,node.getName()),'w')
	fileOut.close()	
	myEmitNodeBeg=Calc_Emit("MyEmitNode_beg_%s"%(node.getName()),True,"%s/emmit_beg_%s.txt"%(outputDirectory,node.getName()))
	myEmitNodeEnd=Calc_Emit("MyEmitNode_end_%s"%(node.getName()),True,"%s/emmit_end_%s.txt"%(outputDirectory,node.getName()))
	node.addChildNode(myEmitNodeBeg,AccNode.ENTRANCE)
	node.addChildNode(myEmitNodeEnd,AccNode.EXIT)	

#================Do some turns===========================================
for i in range(args.turns):
	teapot_latt.trackBunch(b, paramsDict)


#===========Dump bunch infomration=======================================
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
print "Stop."



