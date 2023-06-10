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

print "Start."
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--fileName", dest='fileName', default="outputAddMagnetInjectionRegionAllNodesOG.txt", help="file to print node info into")
parser.add_argument("--nParts",type=int, dest='nParts', default=1, help="number of particles")
parser.add_argument("--turns",type=int, dest='turns', default=1, help="number of complete orbits")
#parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=66, help="What node to monitor")
parser.add_argument("--nodeMonitor",type=int, dest='nodeMonitor', default=-1, help="What node to monitor")
parser.add_argument("--printNodes",type=bool, dest='printNodes', default=False, help="print node list")
parser.add_argument("--doDipoleKickers",type=bool, dest='doDipoleKickers', default=True, help="print node list")
parser.add_argument("--x",type=float, dest='x', default=0, help="x position")
parser.add_argument("--px",type=float, dest='px', default=0, help="xp rad")
parser.add_argument("--y",type=float, dest='y', default=0, help="y position")
parser.add_argument("--py",type=float, dest='py', default=0, help="yp rad")
parser.add_argument("--z",type=float, dest='z', default=0, help="z position")
parser.add_argument("--pz",type=float, dest='pz', default=0, help="dE position")
parser.add_argument("--doNormalKickers",type=bool, dest='doNormalKickers', default=False, help="print node list")
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-1.015811, help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-1.000257, help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-0.9997909, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-1.008058, help="scaleChicane13")
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-0.9804517, help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-0.9934917 , help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-1.009009, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-1.042836, help="scaleChicane13")
parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=-1.047686, help="scaleChicane10")
parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=-0.9969892 , help="scaleChicane11")
parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=-1.00911, help="scaleChicane12")
parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=-0.9775873, help="scaleChicane13")
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=1., help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=1., help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=1., help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=1., help="scaleChicane13")
#parser.add_argument("--scaleChicane10",type=float, dest='scaleChicane10', default=1.017866, help="scaleChicane10")
#parser.add_argument("--scaleChicane11",type=float, dest='scaleChicane11', default=1.102997, help="scaleChicane11")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=0.6130163, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=0.8541834, help="scaleChicane13")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=0.8739071, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=0.9329522, help="scaleChicane13")
#parser.add_argument("--scaleChicane12",type=float, dest='scaleChicane12', default=1, help="scaleChicane12")
#parser.add_argument("--scaleChicane13",type=float, dest='scaleChicane13', default=1, help="scaleChicane13")
parser.add_argument("--scaleChicane10b",type=float, dest='scaleChicane10b', default=1, help="scaleChicane10b")
parser.add_argument("--scaleChicane11b",type=float, dest='scaleChicane11b', default=1, help="scaleChicane11b")
#parser.add_argument("--scaleChicane12b",type=float, dest='scaleChicane12b', default=1.155756, help="scaleChicane12b")
#parser.add_argument("--scaleChicane13b",type=float, dest='scaleChicane13b', default=1.109332, help="scaleChicane13b")
parser.add_argument("--scaleChicane12b",type=float, dest='scaleChicane12b', default=1, help="scaleChicane12b")
parser.add_argument("--scaleChicane13b",type=float, dest='scaleChicane13b', default=1, help="scaleChicane13b")
args = parser.parse_args()
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
#chicaneScale10=1.029949
#chicaneScale11=1.089071
#chicaneScale12=0.968613
#chicaneScale13=1.062516
#chicaneScale10=1.0
#chicaneScale11=1.0
#chicaneScale12=1.0
#chicaneScale13=1.0
chicaneScale10=args.scaleChicane10*args.scaleChicane10b
chicaneScale11=args.scaleChicane11*args.scaleChicane11b
chicaneScale12=args.scaleChicane12*args.scaleChicane12b
chicaneScale13=args.scaleChicane13*args.scaleChicane13b

print "chicaneScale10= ",chicaneScale10
print "chicaneScale11= ",chicaneScale11
print "chicaneScale12= ",chicaneScale12
print "chicaneScale13= ",chicaneScale13
teapot_latt = teapot.TEAPOT_Ring()
print "Read MAD."
teapot_latt.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers.LAT","RING")
#teapot_latt.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly.LAT","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())
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
			print "total length= ",node.getLength()
			print "segment length= ",node.getLength(3)
			#myDipole_DH_A11=YDipole("Dipole_DH_A11")
			myDipole_DH_A11=GeneralDipole("Dipole_DH_A11")
			myDipole_DH_A11.setMagneticFieldStrength(1.0)
			myDipole_DH_A11.setFieldDirection(math.pi/2)
			myDipole_DH_A11.setEffLength(.0254)
			node.addChildNode(myDipole_DH_A11,AccNode.BODY,3)
		if (node.getName().strip()=="DH_A12"):
			node.setnParts(numberOfParts_DH_A12)
			print "total length= ",node.getLength()
			print "segment length= ",node.getLength(3)
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

strength_chicane10 = -0.041456*chicaneScale10
strength_chicane11 = 0.052434*chicaneScale11
strength_chicane12 = 0.0298523*chicaneScale12
strength_chicane13 = -0.0398609*chicaneScale13

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

chicane10.setParam("kx", strength_chicane10)
chicane11.setParam("kx", strength_chicane11)
chicane12.setParam("kx", strength_chicane12)
chicane13.setParam("kx", strength_chicane13)

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
print "chicane10= ",chicane10.getParam("kx")


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

print "(x,px,y,py,z,pz)= (%f,%f,%f,%f,%f,%f) " %(b.getSyncParticle().x(),b.getSyncParticle().px(),b.getSyncParticle().y(),b.getSyncParticle().py(),b.getSyncParticle().z(),b.getSyncParticle().pz())
print "gamma= %f"%b.getSyncParticle().gamma()
print "beta= %f"%b.getSyncParticle().beta()
print "momentum= %f"%b.getSyncParticle().momentum()
inject = InjectParts(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
inject.addParticles()


teapot_latt.initialize()



#-------------------------------
#  Lattice is ready
#-------------------------------

fileOut=open(args.fileName,'w')
fileOut.close()
myPrintNode=Print_Node("MyPrintNode",True,args.fileName)

nodes = teapot_latt.getNodes()
if args.nodeMonitor >= 0:
	nodes[args.nodeMonitor].addChildNode(myPrintNode,AccNode.EXIT)
else:
	for node in nodes:
		node.addChildNode(myPrintNode,AccNode.EXIT)	
i = 0
path_length=0
for node in nodes:
	pass
	if node.getName().strip() == "DH_A12":
		print node.getName().strip()
		print node.getnParts()
		#node.setnParts(10)
	if args.printNodes==True:
		path_length=path_length+node.getLength()
		print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node], " path_length= ",path_length
		#print "There are ", node.getNumberOfBodyChildren()," child nodes."
		i=i+1

#================Do some turns===========================================
for i in range(args.turns):
	teapot_latt.trackBunch(b, paramsDict)


#===========Dump bunch infomration=======================================
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
print "Stop."



