#!/usr/bin/env python

#--------------------------------------------------------
# Test of the Overlapping Quads Model
#--------------------------------------------------------

import math
import sys
import os

from orbit.py_linac.overlapping_fields import EngeFunction, SNS_MEBT_OverlappingQuadsSubst
from orbit.py_linac.overlapping_fields import OverlappingQuadsController
from orbit.py_linac.overlapping_fields import getGlobalField

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

#-------------------------------------------
# START of Script
#-------------------------------------------

names = ["MEBT",]

py_orbit_sns_home = "../"

#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- the XML file name with the structure
xml_file_name = py_orbit_sns_home+"sns_linac_xml/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print "Linac lattice is ready. L=",accLattice.getLength()

length = accLattice.getLength()
step = 0.001
n_points = int(length/step) + 1
res_z_G0_G_arr = []
for ind in range(n_points):
	z = step*ind
	G0 = getGlobalField(accLattice,z)
	res_z_G0_G_arr.append([z,G0,0.])

#--------------------------------------------------------
# Replace overlapping quads with the special nodes
#--------------------------------------------------------

(ovrlp_quads1,ovrlp_quads2) = SNS_MEBT_OverlappingQuadsSubst(accLattice)

print "Linac lattice is ready. L=",accLattice.getLength()

node_pos_dict = accLattice.getNodePositionsDict()
nodes = accLattice.getNodes()
for node in nodes:
	(posBefore, posAfter) = node_pos_dict[node]
	print "%35s      start stop L = %10.4f %10.4f    %10.4f "%(node.getName(),posBefore,posAfter,(posAfter-posBefore))

for ind in range(n_points):
	z = step*ind
	G = getGlobalField(accLattice,z)
	res_z_G0_G_arr[ind][2] = G

fl_out = open("overlapping_field_test.dat","w")
for [z,G0,G] in res_z_G0_G_arr:
	s = " %12.5f  %12.5f  %12.5f "%(z,G0,G)
	fl_out.write(s+"\n")
fl_out.close()