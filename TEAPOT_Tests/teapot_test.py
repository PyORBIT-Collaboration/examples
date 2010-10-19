from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch

print "Start."

b = Bunch()
b.addParticle(1.0e-3,0.0,0.0,0.0,0.0,0.0)
b.addParticle(0.0,1.0e-3,0.0,0.0,0.0,0.0)
b.addParticle(0.0,0.0,1.0e-3,0.0,0.0,0.0)
b.addParticle(0.0,0.0,0.0,1.0e-3,0.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,1.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,0.0,1.0e-3)
b.compress()

syncPart = b.getSyncParticle()
energy = 1.0                          #energy in GeV
#p = syncPart.energyToMomentum(energy)
#syncPart.pz(p)
syncPart.kinEnergy(energy)

lattice = AccLattice("test_lattice")

elem1 = teapot.DriftTEAPOT("drift1")
elem2 = teapot.QuadTEAPOT("quad1")
elem3 = teapot.QuadTEAPOT("quad2")
elem4 = teapot.BendTEAPOT("bend1")
elem5 = teapot.BendTEAPOT("bend2")
elem6 = teapot.MultipoleTEAPOT("sextupole")

lattice.addNode(elem1)
lattice.addNode(elem2)
lattice.addNode(elem3)
lattice.addNode(elem4)
lattice.addNode(elem5)
lattice.addNode(elem6)

#-----------------------------
# Set TEAPOT nodes parameters
#-----------------------------
elem1.setLength(0.2)
elem2.setLength(0.3)
elem3.setLength(0.4)
elem4.setLength(0.5)
elem5.setLength(0.6)
elem6.setLength(0.7)

elem2.setnParts(5)
elem2.addParam("kq",-0.7)

elem3.setnParts(5)
elem3.addParam("kq",+0.7)

elem4.setnParts(11)
elem4.addParam("theta",+0.1)
elem4.addParam("ea1",+0.01)
elem4.addParam("ea2",+0.02)

elem5.setnParts(11)
elem5.addParam("theta",+0.2)
elem5.addParam("ea1",+0.01)
elem5.addParam("ea2",+0.02)

elem6.setnParts(5)
elem6.getParam("poles").append(2)
elem6.getParam("kls").append(0.28)
elem6.getParam("skews").append(0)

lattice.initialize()

print "==============BEFORE============================"
#b.dumpBunch()
print "=========================================="

#=====track action ============
def bodyAction(paramsDict):
	node = paramsDict["node"]
	node.track(paramsDict)


accContainer = AccActionsContainer()
accContainer.addAction(bodyAction,AccActionsContainer.BODY)

paramsDict = {}
paramsDict["bunch"] = b

lattice.trackActions(accContainer,paramsDict)
print "=============AFTER============================="
b.dumpBunch()
print "=========================================="

print "lattice length=",lattice.getLength()
print "beta=",b.getSyncParticle().beta()
print "TEAPOT time[sec]=",b.getSyncParticle().time()
print "SIMPLE time[sec]=",lattice.getLength()/(b.getSyncParticle().beta()*2.99792458e+8)
print "Stop."

#==============BEFORE============================
#==========================================
#=============AFTER=============================
#% PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES
#% BUNCH_ATTRIBUTE_DOUBLE charge   1
#% BUNCH_ATTRIBUTE_DOUBLE classical_radius   1.5347e-18
#% BUNCH_ATTRIBUTE_DOUBLE macro_size   0
#% BUNCH_ATTRIBUTE_DOUBLE mass   0.938272
#%  SYNC_PART_COORDS 0 0 0  x, y, z positions in [m]
#%  SYNC_PART_MOMENTUM 0 0 1.696037912  px, py, pz momentum component in GeV/c
#%  info only: energy of the synchronous particle [GeV] = 1
#%  info only: momentum of the synchronous particle [GeV/c] = 1.696037912
#%  info only: beta=v/c of the synchronous particle = 0.8750256155
#%  info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle = 2.065788684
#%  SYNC_PART_TIME 1.029253363e-08  time in [sec]
#% x[m] px[rad] y[m] py[rad] z[m]  (pz or dE [GeV])
#0.00082842613 -0.00015325645 0 0 0.00029901071 -1.647656e-16 
#0.0023631085 0.00076945561 0 0 0.00043014195 -3.2953121e-16 
#2.3859816e-08 1.3257047e-07 0.0010236879 4.1172715e-05 4.4922855e-08 0 
#-4.7883174e-07 6.8646524e-07 0.0028841539 0.0010928603 1.5869592e-06 -3.2953121e-16 
#0 0 0 0 1 0 
#0.00023639406 0.00019942888 0 0 -0.00041525426 0.001 
#==========================================
#lattice length= 2.7
#beta= 0.875025615499
#TEAPOT time[sec]= 1.02925336251e-08
#SIMPLE time[sec]= 1.02925336251e-08
#Stop.


