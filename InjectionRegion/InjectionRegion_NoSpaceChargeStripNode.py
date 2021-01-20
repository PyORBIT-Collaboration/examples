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
from KevinPython.printNode import Print_Node
from KevinPython.calculateEmit import Calc_Emit

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

from orbit.utils import Function
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
parser.add_argument("--xOffset",type=float, dest='xOffset', default=0, help="x injection offset")
parser.add_argument("--pxOffset",type=float, dest='pxOffset', default=0, help="px injection offset")
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


#=====Make a Teapot style lattice======
chicaneScale10=args.scaleChicane10
chicaneScale11=args.scaleChicane11
chicaneScale12=args.scaleChicane12
chicaneScale13=args.scaleChicane13

inj_latt = teapot.TEAPOT_Ring()
print "Read MAD."
inj_latt.readMAD("MAD_Injection_Region_Lattice/InjectionRegionOnly_Chicane_Replaced_With_Kickers_onlyInnerMostChicanes.LAT","RING")
print "Lattice=",inj_latt.getName()," length [m] =",inj_latt.getLength()," nodes=",len(inj_latt.getNodes())

ring_latt = teapot.TEAPOT_Ring()
print "Read MAD."
ring_latt.readMAD("MAD_Full_Ring/SNSring_pyOrbitBenchmark_Chicanes_Kickers_startEndof3rdChicane.LAT","RING")
print "Lattice=",ring_latt.getName()," length [m] =",ring_latt.getLength()," nodes=",len(ring_latt.getNodes())

theEffLength=0.0254
fieldStrength=1.0


nodes2 = inj_latt.getNodes()
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
			myDipole_DH_A11=GeneralDipoleStrip("Dipole_DH_A11")
			myDipole_DH_A11.setMagneticFieldStrength(fieldStrength)
			myDipole_DH_A11.setFieldDirection(math.pi/2)
			myDipole_DH_A11.setEffLength(theEffLength)
			node.addChildNode(myDipole_DH_A11,AccNode.BODY,3)
		if (node.getName().strip()=="DH_A12"):
			node.setnParts(numberOfParts_DH_A12)
			print "total length= ",node.getLength()
			print "segment length= ",node.getLength(3)
			#myDipole_DH_A12=YDipole("Dipole_DH_A12")
			myDipole_DH_A12=GeneralDipoleStrip("Dipole_DH_A12")
			myDipole_DH_A12.setMagneticFieldStrength(-fieldStrength)
			myDipole_DH_A12.setFieldDirection(math.pi/2)
			myDipole_DH_A12.setEffLength(theEffLength)
			node.addChildNode(myDipole_DH_A12,AccNode.BODY,3)	
		
print "counterDH_A11=",counterDH_A11
print "counterDH_A12=",counterDH_A12

nodes2 = ring_latt.getNodes()
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
			myDipole_DH_A11.setMagneticFieldStrength(fieldStrength)
			myDipole_DH_A11.setFieldDirection(math.pi/2)
			myDipole_DH_A11.setEffLength(theEffLength)
			node.addChildNode(myDipole_DH_A11,AccNode.BODY,3)
		if (node.getName().strip()=="DH_A12"):
			node.setnParts(numberOfParts_DH_A12)
			print "total length= ",node.getLength()
			print "segment length= ",node.getLength(3)
			#myDipole_DH_A12=YDipole("Dipole_DH_A12")
			myDipole_DH_A12=GeneralDipole("Dipole_DH_A12")
			myDipole_DH_A12.setMagneticFieldStrength(-fieldStrength)
			myDipole_DH_A12.setFieldDirection(math.pi/2)
			myDipole_DH_A12.setEffLength(theEffLength)
			node.addChildNode(myDipole_DH_A12,AccNode.BODY,3)	
		
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

strength_chicane10 = -0.041456*chicaneScale10
strength_chicane11 = 0.052434*chicaneScale11
strength_chicane12 = 0.0298523*chicaneScale12
strength_chicane13 = -0.0398609*chicaneScale13

#strength_chicane12 = -0.041456*chicaneScale10
#strength_chicane13 = 0.052434*chicaneScale11
#strength_chicane10 = 0.0298523*chicaneScale12
#strength_chicane11 = -0.0398609*chicaneScale13

lattlength = ring_latt.getLength()
sp = b.getSyncParticle()
kickerwave = rootTWaveform(sp, lattlength, duration, startamp, endamp)
chicanewave = flatTopWaveform(1.0)

i = 0
path_length=0
print "inj_latt"
nodes = inj_latt.getNodes()
for node in nodes:
	pass
	if node.getName().strip() == "DH_A12":
		print node.getName().strip()
		print node.getnParts()
		#node.setnParts(10)
	if args.printNodes==True:
		path_length=path_length+node.getLength()
		print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%inj_latt.getNodePositionsDict()[node], " path_length= ",path_length
		#print "There are ", node.getNumberOfBodyChildren()," child nodes."
		i=i+1

i = 0
path_length=0
print "ring_latt"
nodes = ring_latt.getNodes()
for node in nodes:
	pass
	if node.getName().strip() == "DH_A12":
		print node.getName().strip()
		print node.getnParts()
		#node.setnParts(10)
	if args.printNodes==True:
		path_length=path_length+node.getLength()
		print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%ring_latt.getNodePositionsDict()[node], " path_length= ",path_length
		#print "There are ", node.getNumberOfBodyChildren()," child nodes."
		i=i+1
nodes = ring_latt.getNodes()
hkick10 = nodes[669]
vkick10 = nodes[671]
hkick11	= nodes[673]
vkick11 = nodes[675]
vkick12 = nodes[15]
hkick12 = nodes[17]
vkick13 = nodes[19]
hkick13	= nodes[21]

chicane10 = nodes[688]
chicane11 = nodes[690]
chicane12 = nodes[693]
chicane13 = nodes[2]

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


nodes = inj_latt.getNodes()
chicane11 = nodes[1]
chicane12 = nodes[4]


chicane11.setParam("kx", strength_chicane11)
chicane12.setParam("kx", strength_chicane12)

chicane11.setWaveform(chicanewave)
chicane12.setWaveform(chicanewave)
#for node in kickernode:
#print "node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
#print "There are ", node.getNumberOfBodyChildren()," child nodes."

#------------------------------
#Initial Distribution Functions
#------------------------------


sp = b.getSyncParticle()
beta= sp.getBeta()
gamma=sp.getGamma()

print "beta= %f"%beta
print "gamma= %f"%gamma

A1=2.47e-6
A2=4.49e9

def constantField(x);
	return fieldStrength

magneticField= Function()
n=1000
maxValue=theEffLength
step=maxValue/n

for i in range(n):
	x = step*i;
	y = constantField(x)
	magneticField.add(x,y)
	
theStrippingFunctions=probabilityStripping(magneticField,n,maxValue,gamma,beta)
notNormalizedFunction=theStrippingFunctions.getNotNormalizedFunction()
NormalizedFunction=theStrippingFunctions.getNormalizedFunction()
InverseFunction=theStrippingFunctions.getInverseFunction()
#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta 
(alphaX,betaX,emittX) = (-1.9620, 0.1831, 0.21)
(alphaY,betaY,emittY) = ( 1.7681, 0.1620, 0.21)
(alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)

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


bunch_in.mass(0.93827231)
bunch_in.macroSize(macrosize)
energy = 1.0 #Gev
bunch_in.getSyncParticle().kinEnergy(energy)
bunch_in.setCharge(-1)

paramsDict = {}
firstChicaneFail = Bunch()
firstChicaneFail.setCharge(-1)
secondChicaneFail = Bunch()
secondChicaneFail.setCharge(0)
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch_in
paramsDict["firstChicaneFail"]=firstChicaneFail
paramsDict["secondChicaneFail"]= secondChicaneFail
lostbunch.addPartAttr("LostParticleAttributes") 

#====Injection and foil aperature============

xmin = xcenterpos - 0.0085
xmax = xcenterpos + 0.0085
ymin = ycenterpos - 0.0080
ymax = ycenterpos + 0.100


xFunc=notRandom(args.x,args.px)
yFunc=notRandom(args.y,args.py)
lFunc=notRandom(args.z,args.pz)

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
#injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
#addTeapotInjectionNode(teapot_latt, 0., injectnode) 

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
ZtoPhi = 2.0 * math.pi / lattlength;
dESync = 0.0
RF1HNum = 1.0
RF1Voltage = 0.000016
RF1Phase = 0.0
RF2HNum = 2.0
RF2Voltage = -0.000003
RF2Phase = 0.0
length = 0.0

rf1_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RF1HNum, RF1Voltage, RF1Phase, length, "RF1")
rf2_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RF2HNum, RF2Voltage, RF2Phase, length, "RF2")
position1 = 196.0
position2 = 196.5
RFLatticeModifications.addRFNode(teapot_latt, position1, rf1_node)
RFLatticeModifications.addRFNode(teapot_latt, position2, rf2_node)

#----------------------------------------------
#make 2.5D space charge calculator
#----------------------------------------------
#set boundary


#-----------------------------------------------
# Add longitudinal space charge node with Imped
#-----------------------------------------------



#-------------------------------
#  Lattice is ready
#-------------------------------

fileOut=open(args.fileName,'w')
fileOut.close()
fileOut=open(args.fileName2,'w')
fileOut.close()
myPrintNode=Print_Node("MyPrintNode",True,args.fileName)

myEmitNode=Calc_Emit("MyEmitNode",True,args.fileName2)

nodes = teapot_latt.getNodes()
nodes[args.nodeMonitor].addChildNode(myPrintNode,AccNode.EXIT)
nodes[args.nodeMonitor].addChildNode(myEmitNode,AccNode.EXIT)
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



