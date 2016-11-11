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

from bunch import Bunch


names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]


#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- the XML file name with the structure
xml_file_name = "../../sns_linac_xml/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print "Linac initial lattice is ready.  L=",accLattice.getLength()

#---- axis fields files location 
dir_location = "../../sns_rf_fields/"

Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,dir_location,["MEBT","CCL1","CCL2","CCL3","CCL4","SCLMed"])
#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,dir_location,["SCLHigh",])
#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,dir_location,["SCLMed",])
#Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,dir_location,["CCL1",])

print "Linac modified lattice is ready. L=",accLattice.getLength()

print "Stop."