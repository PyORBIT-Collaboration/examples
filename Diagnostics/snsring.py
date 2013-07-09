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
from orbit.collimation import TeapotCollimatorNode, addTeapotCollimatorNode
from orbit.diagnostics import StatLats
from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotStatLatsNodeSet, addTeapotMomentsNodeSet

print "Start."

#=====Main bunch parameters============
b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b
lostbunch.addPartAttr("LostParticleAttributes") 

#=====Make a Teapot style lattice======

teapot_latt = teapot.TEAPOT_Lattice()
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
#xcenterpos = 0.0468
xcenterpos = 0.0
xcentermom = 0.00
#ycenterpos = 0.0492
ycenterpos = 0.0
ycentermom = 0.00
tailfrac = 0.1  # 10% in tails
tailfac = 1.5   # tail emittance is 50% greater than core

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

xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom, tailfrac, tailfac)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom, tailfrac, tailfac)
lFunc = SNSESpreadDist(lattlength, zmin, zmax, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams)

#====Injection and foil aperature============

xmin = xcenterpos - 0.0085
xmax = xcenterpos + 0.0085
ymin = ycenterpos - 0.0080
ymax = ycenterpos + 0.100

#=================Add the injection node and foil node==  ==============

nparts = 20.
injectparams = (xmin, xmax, ymin, ymax)
injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)
addTeapotInjectionNode(teapot_latt, 0., injectnode) 

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
addTeapotCollimatorNode(teapot_latt, 0.5, collimator)

#-------------------------------
#  Lattice is ready
#-------------------------------

nodes = teapot_latt.getNodes()
i = 0
for node in nodes:
	print i, " node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
	print "There are ", node.getNumberOfBodyChildren()," child nodes."
	i=i+1

tunes = TeapotTuneAnalysisNode("tune_analysis")
addTeapotDiagnosticsNode(teapot_latt, 36.066, tunes)


#================Do some turns===========================================

teapot_latt.trackBunch(b, paramsDict)
teapot_latt.trackBunch(b, paramsDict)

for i in xrange(39):
	teapot_latt.trackBunch(b, paramsDict)
addTeapotStatLatsNodeSet(teapot_latt, "statlats")
addTeapotMomentsNodeSet(teapot_latt, "moments", 3)
print 'last turn'
teapot_latt.trackBunch(b, paramsDict)
b.dumpBunch("bunch_36.dat")


#===========Dump bunch infomration=======================================
print "Stop."



