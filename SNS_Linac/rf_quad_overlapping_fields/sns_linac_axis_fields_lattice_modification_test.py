#! /usr/bin/env python

"""
This script will check the lattice modification for the AxisFieldRF_Gap nodes
"""

import sys
import math
import random
import time

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.py_linac.lattice import Drift, AxisFieldRF_Gap

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.overlapping_fields import GetGlobalQuadGradient
from orbit.py_linac.overlapping_fields import GetGlobalRF_AxisField
from orbit.py_linac.overlapping_fields import SNS_EngeFunctionFactory

from bunch import Bunch


names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
names = ["MEBT","DTL1"]
#names = ["MEBT",]

#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print "Linac initial lattice is ready.  L=",accLattice.getLength()


#---- magn_field_arr[[z,g0,g1],...]
magn_field_arr = []
step = 0.001
n_points = int(accLattice.getLength()/step)
step = accLattice.getLength()/(n_points-1)
for ip in range(n_points):
	z = step*ip
	g0 = GetGlobalQuadGradient(accLattice,z)
	g1 = 0.
	Ez = 0.
	magn_field_arr.append([z,g0,g1,Ez])

#---- longitudinal step along the distributed fields lattice
z_step = 0.005 


#---- axis fields files location 
dir_location = "../sns_rf_fields/"

#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["MEBT","CCL1","CCL2","CCL3","CCL4","SCLMed"])
#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["SCLHigh",])
#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["SCLMed",])
#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["CCL1",])

#Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6"])
#Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1"])
Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1"],[],SNS_EngeFunctionFactory)
#Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT",])
#Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["CCL1","CCL2","CCL3","CCL4"])
#Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["DTL2",])

#Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["MEBT","DTL1",],[],SNS_EngeFunctionFactory)

print "Linac modified lattice is ready. L=",accLattice.getLength()

for ip in range(n_points):
	[z,g0,g1,Ez] = magn_field_arr[ip]
	g1 = GetGlobalQuadGradient(accLattice,z)
	Ez = GetGlobalRF_AxisField(accLattice,z)
	magn_field_arr[ip] = [z,g0,g1,Ez]

fl_out = open("g_fields.dat","w")
for ip in range(n_points):
	[z,g0,g1,Ez] = magn_field_arr[ip]
	fl_out.write(" %12.6f %12.6f  %12.6f    %12.5g "%(z,g0,g1,Ez)+"\n")
fl_out.close()


print "Stop."