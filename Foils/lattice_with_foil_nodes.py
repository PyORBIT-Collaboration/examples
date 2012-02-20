##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# several collimation nodes
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.foils import TeapotFoilNode
from orbit.foils import addTeapotFoilNode

print "Start."

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("MAD_Lattice/LATTICE","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

n_drifts = 0
drift_length = 0.
for node in teapot_latt.getNodes():
	if(isinstance(node,DriftTEAPOT)):
		n_drifts += 1
		drift_length += node.getLength()
		
print "number of drifts =",n_drifts
print "total drift length =",drift_length
print "============= nodes inside the region ==========="
pos_start = 16.5
pos_stop  = 20.5
# print all nodes around the specified position
for node in teapot_latt.getNodes():
	(node_pos_start,node_pos_stop) = teapot_latt.getNodePositionsDict()[node]
	if(node_pos_start > pos_start and node_pos_stop < pos_stop):
		print "node=",node.getName()," type=",node.getType(),"  pos=",node_pos_start," L=",node.getLength()


xmin = -0.100
xmax = 0.100
ymin = -0.100
ymax = 0.100
# Below is 1000 times the width of normal foil but will do only one turn.
thick = 400000

foil = TeapotFoilNode(xmin, xmax, ymin, ymax, thick, "Foil 1")

addTeapotFoilNode(teapot_latt, 18.5,foil) 
print "===========Lattice modified ======================================="
print "New Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

n_drifts = 0
drift_length = 0.
for node in teapot_latt.getNodes():
	if(isinstance(node,DriftTEAPOT)):
		n_drifts += 1
		drift_length += node.getLength()
		
print "number of drifts =",n_drifts
print "total drift length =",drift_length

print "============= nodes inside the region ==========="
# print all nodes around the specified position
for node in teapot_latt.getNodes():
	(node_pos_start,node_pos_stop) = teapot_latt.getNodePositionsDict()[node]
	if(node_pos_start > pos_start and node_pos_stop < pos_stop):
		print "node=",node.getName()," type=",node.getType(),"  pos=",node_pos_start," L=",node.getLength()


#------------------------------
# Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "Benchmark_Foil"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
# get initial bunch from ORBIT_MPI input
bunch_orbit_to_pyorbit(teapot_latt.getLength(), energy, "Bm_KV_Uniform_10000",b)
#b.readBunch("parts_in.dat")

b.getSyncParticle().kinEnergy(energy)

#=====track bunch through Collimator Node============
paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b

lostbunch.addPartAttr("LostParticleAttributes") 

scatterchoice = 1
foil.setScatterChoice(scatterchoice)
#collimator.trackBunch(b,paramsDict)
foil.track(paramsDict)

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
lostbunch.dumpBunch()
print "Stop."



