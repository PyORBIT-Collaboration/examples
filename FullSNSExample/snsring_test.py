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
from injection import InjectParts
from injection import JohoTransverse, JohoLongitudinal, SNSESpreadDist
from kickernodes import XKicker, YKicker
from kickernodes import rootTWaveform, flatTopWaveform
from kickernodes import TeapotXKickerNode, TeapotYKickerNode,addTeapotKickerNode
from orbit.foils import TeapotFoilNode, addTeapotFoilNode
from foil import Foil
from orbit.collimation import TeapotCollimatorNode, addTeapotColimatorNode
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D, Boundary2D
from spacecharge import LSpaceChargeCalc
from orbit.space_charge.sc1d import addLongitudinalSpaceChargeNode, SC1D_AccNode

print "Start."

#=====Main bunch parameters============
intensity = 7.8e13
turns = 1000.0
macrosperturn = 260
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
teapot_latt.readMAD("MAD_Lattice/RealInjection/SNSring_pyOrbitBenchmark.LAT","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

#====Add the injection kickers======

duration = 0.001
startamp = 1.0
endamp = 0.58
deltapos = 0.001
 
strength_hkicker10 = 7.211536E-03
strength_hkicker13 = strength_hkicker10
strength_hkicker11 = -2.278306E-03
strength_hkicker12 = strength_hkicker11
strength_vkicker10 = 4.188402E-03
strength_vkicker13 = strength_vkicker10
strength_vkicker11 = -2.118213E-03
strength_vkicker12 = strength_vkicker11

lattlength = teapot_latt.getLength()
sp = b.getSyncParticle()
kickerwave = rootTWaveform(sp, lattlength, duration, startamp, endamp)

nodes = teapot_latt.getNodes()
hkick10 = nodes[671]
vkick10 = nodes[673]
hkick11	= nodes[675]
vkick11 = nodes[677]
vkick12 = nodes[17]
hkick12 = nodes[19]
vkick13 = nodes[21]
hkick13	= nodes[23]

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
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom)
lFunc = SNSESpreadDist(lattlength, zmin, zmax, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams)

#====Injection and foil aperature============

xmin = xcenterpos - 0.0085
xmax = xcenterpos + 0.0085
ymin = ycenterpos - 0.0080
ymax = ycenterpos + 0.100

#=================Add the injection node and foil node==  ==============

nparts = macrosperturn
injectparams = (xmin, xmax, ymin, ymax)
injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
addTeapotInjectionNode(teapot_latt, 0., injectnode) 

thick = 400.0
foil = TeapotFoilNode(xmin, xmax, ymin, ymax, thick, "Foil 1")
scatterchoice = 0
foil.setScatterChoice(scatterchoice)
addTeapotFoilNode(teapot_latt,0.000001,foil)

#----------------------------------------------
# Add one black absorber collimator to act like
# an aperture
#----------------------------------------------
colllength = 0.00001
ma = 9
density_fac = 1.0
shape = 1
radius = 0.110

collimator = TeapotCollimatorNode(colllength, ma, density_fac, shape, radius, 0., 0., 0., 0., "Collimator 1")
addTeapotColimatorNode(teapot_latt, 0.5, collimator)

#----------------------------------------------
#make 2.5D space charge calculator
#----------------------------------------------

sizeX = 64   #number of grid points in horizontal direction
sizeY = 64  #number of grid points in vertical direction
sizeZ = 1     #number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)
sc_path_length_min = 0.00000001
#scLatticeModifications.setSC2p5DAccNodes(teapot_latt, sc_path_length_min,calc2p5d)

#-----------------------------------------------
# Add longitudinal space charge node with Imped
#-----------------------------------------------

b_a = 10.0/3.0
length = lattlength
nMacrosMin = 1000
useSpaceCharge = 1
nBins= 128     #number of longitudinal slices in the 1D space charge solver
position = 64.0

#SNS Longitudinal Impedance tables. EKicker impedance from private communication
# with J.G. Wang. Seems to be for 7 of the 14 kickers (not sure why).
# Impedance in Ohms/n. Kicker and RF impedances are inductive with real part positive and imaginary is negative by Chao definition.

ZL_EKicker = [complex(42., -182),
			  complex(35, -101.5),
			  complex(30.3333, -74.6667),
			  complex(31.5, -66.5),
			  complex(32.2,-57.4),
			  complex(31.5, -51.333),
			  complex(31, -49),
			  complex(31.5, -46.375),
			  complex(31.8889, -43.556),
			  complex(32.9, -40.6),
			  complex(32.7273, -38.18),
			  complex(32.25, -35.58),
			  complex(34.46, -32.846),
			  complex(35, -30.5),
			  complex(35.4667, -28.),
			  complex(36.75, -25.81),
			  complex(36.647, -23.88),
			  complex(36.944, -21.1667),
			  complex(36.474, -20.263),
			  complex(36.4, -18.55),
			  complex(35.333, -17),
			  complex(35, -14.95),
			  complex(33.478, -13.69),
			  complex(32.375, -11.67),
			  complex(30.8, -10.08),
			  complex(29.615, -8.077),
			  complex(28.519, -6.74),
			  complex(27.5, -5),
			  complex(26.552, -4.103),
			  complex(25.433, -3.266),
			  complex(24.3871, -2.7),
			  complex(23.40625, -2.18)]

ZL_RF = [complex(0.0, 0.0),
		 complex(0.750, 0.0),
		 complex(0.333,0.0),
		 complex(0.250, 0.0),
		 complex(0.200, 0.0),
		 complex(0.167, 0.0),
		 complex(3.214, 0.0),
		 complex(0.188, 0.0),
		 complex(0.167, 0.0),
		 complex(0.150, 0.0),
		 complex(1.000, 0.0),
		 complex(0.125, 0.0),
		 complex(0.115, 0.0),
		 complex(0.143, 0.0),
		 complex(0.333, 0.0),
		 complex(0.313, 0.0),
		 complex(0.294, 0.0),
		 complex(0.278, 0.0),
		 complex(0.263, 0.0),
		 complex(0.250, 0.0),
		 complex(0.714, 0.0),
		 complex(0.682, 0.0),
		 complex(0.652, 0.0),
		 complex(0.625, 0.0),
		 complex(0.600, 0.0),
		 complex(0.577, 0.0),
		 complex(0.536, 0.0),
		 complex(0.536, 0.0),
		 complex(0.517, 0.0),
		 complex(0.500, 0.0),
		 complex(0.484, 0.0),
		 complex(0.469, 0.0)]


Z = []
for i in range(0,32):
	zk = ZL_EKicker[i]
	zrf = ZL_RF[i]
	zreal = (zk.real/1.75 + zrf.real)
	zimag = (zk.imag/1.75 + zrf.imag)
	Z.append(complex(zreal, zimag))

sc1Dnode = SC1D_AccNode(b_a,length, nMacrosMin,useSpaceCharge,nBins)
sc1Dnode.assignImpedance(Z);
addLongitudinalSpaceChargeNode(teapot_latt, position, sc1Dnode)

#-------------------------------
#  Lattice is ready
#-------------------------------

nodes = teapot_latt.getNodes()
i = 0
for node in nodes:
	print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
	print "There are ", node.getNumberOfBodyChildren()," child nodes."
	i=i+1

#================Do some turns===========================================

teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_1.dat")

for i in xrange(99):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_100.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_200.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_300.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_400.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_500.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_600.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_700.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_800.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_900.dat")

for i in xrange(100):
	teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch_1000.dat")

#===========Dump bunch infomration=======================================
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
print "Stop."



