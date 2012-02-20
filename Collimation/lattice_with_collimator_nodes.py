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

from orbit.collimation import TeapotCollimatorNode
from orbit.collimation import addTeapotColimatorNode

print "Start."

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("./MAD_Lattice/LATTICE","RING")
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


length = 0.5
ma = 9
density_fac = 1.0
shape = 1
a = 0
b = 0
c = 0
d = 0
angle = 0

collimator = TeapotCollimatorNode(length, ma, density_fac, shape, a, b, c, d, angle, "Collimator 1")

addTeapotColimatorNode(teapot_latt, 18.5,collimator) 
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

print "Stop."


