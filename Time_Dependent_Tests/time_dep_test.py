import sys
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement, MAD_LattLine
from orbit.time_dep import time_dep, waveform
from bunch import Bunch
"""
#test Parser 'Dict'
parser = MAD_Parser()
parser.parse("rcs1124.dat")

elems = parser.getMAD_Elements()
lines = parser.getMAD_Lines()
print "number of elems",len(elems)
print "number of lines",len(lines)

elemsDict = parser.getMAD_ElementsDict()
linesDict = parser.getMAD_LinesDict()
print "number of elemsDict",len(elemsDict)
print "number of linesDict",len(linesDict)

ring = linesDict["R11"]
elem = elemsDict["LQE"]
#ring = lines[0]
#elem = elems[2]
num = 0
for ring_elem in ring.getElements():
	if ring_elem.getName() == elem.getName():
		num+=1
print "The line",ring.getName()," includes", elem.getName(), " =",num

#test TIME_DEP_Lattice
b = Bunch()
b.addParticle(1.0e-3,0.0,0.0,0.0,0.0,0.0)
b.addParticle(0.0,1.0e-3,0.0,0.0,0.0,0.0)
b.addParticle(0.0,0.0,1.0e-3,0.0,0.0,0.0)
b.addParticle(0.0,0.0,0.0,1.0e-3,0.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,1.0e-3,0.0)
b.addParticle(0.0,0.0,0.0,0.0,0.0,1.0e-3)
b.compress()
syncPart = b.getSyncParticle()
energy = 1.0                          
syncPart.kinEnergy(energy)
"""
time_dep_latt = time_dep.TIME_DEP_Lattice()
time_dep_latt.readMAD("rcs1124.dat","R1")
time_dep_latt.setLatticeOrder()
for node in time_dep_latt.getNodes():
	print "node", node.getName()," = ", node.getType(), "param=", node.getParamsDict()

WaveForm01 = waveform.MagnetWaveform()
WaveForm01.setConstant(2)
time_dep_latt.setTimeDepNode("R1SF02",1,WaveForm01)
node = time_dep_latt.getTimeDepNode("R1SF02",1)
print  "node", node.getName()," = ", node.getType(), "param=", node.getParamsDict()

"""
time_dep_latt.trackBunch(b)
print "Lattice=",time_dep_latt.getName()," length [m] =",time_dep_latt.getLength()
print "Stop."
b.dumpBunch()
"""
