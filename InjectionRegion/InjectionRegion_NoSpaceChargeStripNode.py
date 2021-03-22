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
parser.add_argument("--printNodes",type=bool, dest='printNodes', default=False, help="print node list")
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

#------------------------------
#Initial Distribution Functions
#------------------------------

e_kin_ini = 1.0 # in [GeV]
mass =  0.93827231 #0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta


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


bunch_in.mass(mass) #mass
bunch_in.macroSize(macrosize)
energy = e_kin_ini # 1.0 #Gev
bunch_in.getSyncParticle().kinEnergy(energy)
bunch_in.charge(-1)

paramsDict = {}
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
InverseFunction=theStrippingFunctions.getInverseFunction()

print "accumlatedSum->getY(maxValue)"
print accumlatedSum.getY(maxValue)
print "CDF->getY(maxValue)"
print CDF.getY(maxValue)
print "InverseFunction.getY(maxValue)"
print InverseFunction.getY(.954)
print deltaxp_rigidity.getY(maxValue)
print deltax_rigidity.getY(maxValue)

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
			myDipole_DH_A11.setFunctionCDF(CDF)
			myDipole_DH_A11.setFunctionInverse(InverseFunction)
			myDipole_DH_A11.setFunctionXPRigidity(deltaxp_rigidity)
			myDipole_DH_A11.setFunctionXRigidity(deltax_rigidity)
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
			myDipole_DH_A12.setFunctionCDF(CDF)
			myDipole_DH_A12.setFunctionInverse(InverseFunction)
			myDipole_DH_A12.setFunctionXPRigidity(deltaxp_rigidity)
			myDipole_DH_A12.setFunctionXRigidity(deltax_rigidity)
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
sp = bunch_in.getSyncParticle()
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





#-------------------------------
#  Lattice is ready
#-------------------------------


#myPrintNode=Print_Node("MyPrintNode",True,args.fileName)
fileOut=open("print_beg.txt",'w')
fileOut.close()
fileOut=open("print_beg_b23.txt",'w')
fileOut.close()
fileOut=open("print_mid_b23.txt",'w')
fileOut.close()
fileOut=open("print_end_b23.txt",'w')
fileOut.close()
fileOut=open("print_beg_DH12.txt",'w')
fileOut.close()
fileOut=open("print_end_DH12.txt",'w')
fileOut.close()
myPrintNode_beg=Print_Node("MyPrintNode_Beg",True,"print_beg.txt")
myPrintNode_beg_b23=Print_Node("MyPrintNode_beg_b23",True,"print_beg_b23.txt")
myPrintNode_mid_b23=Print_Node("MyPrintNode_mid_b23",True,"print_mid_b23.txt")
myPrintNode_end_b23=Print_Node("MyPrintNode_end_b23",True,"print_end_b23.txt")
myPrintNode_beg_DH12=Print_Node("MyPrintNode_beg_DH12",True,"print_beg_DH12.txt")
myPrintNode_end_DH12=Print_Node("MyPrintNode_end_DH12",True,"print_end_DH12.txt")



fileOut=open("emmit_beg.txt",'w')
fileOut.close()
fileOut=open("emmit_beg_b23.txt",'w')
fileOut.close()
fileOut=open("emmit_mid_b23.txt",'w')
fileOut.close()
fileOut=open("emmit_end_b23.txt",'w')
fileOut.close()
fileOut=open("emmit_beg_DH12.txt",'w')
fileOut.close()
fileOut=open("emmit_end_DH12.txt",'w')
fileOut.close()
myEmitNode_beg=Calc_Emit("MyEmitNode_Beg",True,"emmit_beg.txt")
myEmitNode_beg_b23=Calc_Emit("myEmitNode_beg_b23",True,"emmit_beg_b23.txt")
myEmitNode_mid_b23=Calc_Emit("myEmitNode_mid_b23",True,"emmit_mid_b23.txt")
myEmitNode_end_b23=Calc_Emit("myEmitNode_end_b23",True,"emmit_end_b23.txt")
myEmitNode_beg_DH12=Calc_Emit("myEmitNode_beg_DH12",True,"emmit_beg_DH12.txt")
myEmitNode_end_DH12=Calc_Emit("myEmitNode_end_DH12",True,"emmit_end_DH12.txt")

nodes = inj_latt.getNodes()
nodes[0].addChildNode(myPrintNode_beg,AccNode.ENTRANCE)
nodes[0].addChildNode(myEmitNode_beg,AccNode.ENTRANCE)
#nodes[args.nodeMonitor].addChildNode(myEmitNode,AccNode.EXIT)
i = 0
path_length=0
for node in nodes:
	pass
	if node.getName().strip() == "DB23":
		node.setnParts(2)
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
		#node.setnParts(10)
	if args.printNodes==True:
		path_length=path_length+node.getLength()
		print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%ring_latt.getNodePositionsDict()[node], " path_length= ",path_length
		#print "There are ", node.getNumberOfBodyChildren()," child nodes."
		i=i+1

#================Do some turns===========================================
inj_latt.trackBunch(bunch_in, paramsDict)
circulating_bunch=Bunch()
for i in range(args.turns):
	bunch_in.copyBunchTo(circulating_bunch)
	#ring_latt.trackBunch(circulating_bunch, paramsDict)
	pass


#===========Dump bunch infomration=======================================
bunch_pyorbit_to_orbit(inj_latt.getLength(), bunch_in, "mainbunch.dat")
bunch_pyorbit_to_orbit(inj_latt.getLength(), lostbunch, "lostbunch.dat")
print "Stop."



