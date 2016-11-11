#! /usr/bin/env python

"""
This script will check the axis field lengths of RF gap data in the SNS linac lattice,
and print out RF gaps that cover the neighboring nodes which are not Drifts.
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
from orbit.py_linac.lattice import Drift

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator



class RF_AxisFieldsStore:
	
	"""
	The dictionary with the axis field Functions 
	with the BaseR_Gap node names as keys.
	The keys in this dictionary are the file names from BaseRfGap instances.
	"""
	static_axis_field_dict = {}
	
	def __init__(self):
		pass
	
	@classmethod
	def addField(cls,rf_gap,dir_location = ""):
		if(cls.static_axis_field_dict.has_key(rf_gap.getName())): return
		fl_name = dir_location + rf_gap.getParam("EzFile")
		comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
		data_type = mpi_datatype.MPI_DOUBLE
		rank = orbit_mpi.MPI_Comm_rank(comm)
		main_rank = 0
		x_arr = []
		y_arr = []
		if(rank == 0):
			fl_in = open(fl_name,"r")
			lns = fl_in.readlines()
			fl_in.close()
			for ln in lns:
				res_arr = ln.split()
				if(len(res_arr) == 2):
					x = float(res_arr[0])
					y = float(res_arr[1])
					x_arr.append(x)		
					y_arr.append(y)	
		x_arr = orbit_mpi.MPI_Bcast(x_arr,data_type,main_rank,comm)
		y_arr = orbit_mpi.MPI_Bcast(y_arr,data_type,main_rank,comm)
		function = Function()
		for ind in range(len(x_arr)):
			function.add(x_arr[ind],y_arr[ind])
		cls.static_axis_field_dict[rf_gap.getName()] = function
		
	@classmethod
	def getAxisFieldFunction(cls,gap_name):
		if(cls.static_axis_field_dict.has_key(gap_name)):
			return cls.static_axis_field_dict[gap_name]
		else:
			return None

	@classmethod
	def getSize(cls):
		return len(cls.static_axis_field_dict.keys())
		
names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]

#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- the XML file name with the structure
xml_file_name = "../../sns_linac_xml/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print "Linac lattice is ready. L=",accLattice.getLength()

integral = GaussLegendreIntegrator(500)
	
#---- get all rf gaps in the lattice
dir_location = "../../sns_rf_fields/"
rf_gaps = accLattice.getRF_Gaps()
print "debug n rf gaps=",len(rf_gaps)
axis_field_dict = {}
axis_field_func_arr = []
for rf_gap in rf_gaps:
	RF_AxisFieldsStore.addField(rf_gap,dir_location)
	function = RF_AxisFieldsStore.getAxisFieldFunction(rf_gap.getName())
	bad_nodes = []
	axis_field_func_arr.append([rf_gap,function,bad_nodes])
print "debug nGaps=",len(axis_field_func_arr)," store size=",RF_AxisFieldsStore.getSize()

node_pos_dict = accLattice.getNodePositionsDict()

"""
#---- print out distances between RF gaps 
pos_old = 0.
for [rf_gap,function,bad_nodes] in axis_field_func_arr:
	(posBefore, posAfter) = node_pos_dict[rf_gap]
	print "gap=",rf_gap.getName()," dist[mm]= %8.1f"%((posBefore-pos_old)*1000.)
	pos_old = posBefore
sys.exit(1)
"""

delta_max = 0.
(rf_gap_max_0,rf_gap_max_1) = (None,None)

#----------------------------------------------------------------
# This part will print rf gap names that have overlapping fields
# and the value of overlapping.
#----------------------------------------------------------------
for [rf_gap,function,bad_nodes] in axis_field_func_arr:
	x_min0 = function.getMinX()
	x_max0 = function.getMaxX()
	(posBefore0, posAfter0) = node_pos_dict[rf_gap]	
	for [rf_gap1,function1,bad_nodes1] in axis_field_func_arr:
		if(rf_gap1 == rf_gap): continue
		x_min1 = function1.getMinX()
		x_max1 = function1.getMaxX()
		(posBefore1, posAfter1) = node_pos_dict[rf_gap1]
		if(posBefore1 > posBefore0):
			if(posBefore1 + x_min1 < posAfter0 + x_max0):
				delta = posBefore1 + x_min1 - (posAfter0 + x_max0)
				if(math.fabs(delta_max) < math.fabs(delta)):
					delta_max = delta
					(rf_gap_max_0,rf_gap_max_1) = (rf_gap,rf_gap1)
				print "debug overlap gap0=", rf_gap.getName()," gap1=",rf_gap1.getName()," dist[mm]= %12.5g"%(delta*1000.)

if(math.fabs(delta_max) > 0.):
	print "debug max diff=",delta_max," for gaps= ",(rf_gap_max_0.getName(),rf_gap_max_1.getName())

print "=================================================================="

#----------------------------------------------------------------
# This part will print for each rf gap the names of the nodes that
# are covered (even only in part) by the gaps' RF field
#----------------------------------------------------------------
for [rf_gap,function,bad_nodes] in axis_field_func_arr:
	x_min = function.getMinX()
	x_max = function.getMaxX()
	rf_gap_ind = accLattice.getNodeIndex(rf_gap)
	(posBefore0, posAfter0) = node_pos_dict[rf_gap]
	bad_nodes = []
	max_dist_ovrlp = 0.
	for node in accLattice.getNodes():
		(posBefore, posAfter) = node_pos_dict[node]
		if(posAfter < posBefore0):
			dist = posAfter - ( posBefore0 + x_min)
			if(dist > 0):
				if(not isinstance(node,Drift)):
					bad_nodes.append([node,dist])
		if(posBefore > posAfter0):
			dist = posBefore - (posAfter0 + x_max)
			if(dist < 0):
				if(not isinstance(node,Drift)):
					bad_nodes.append([node,dist])
	if(len(bad_nodes) > 0):
		print "----------- RF gap=",rf_gap.getName()
		for [node,dist] in bad_nodes:
			print "   debug node = ",node.getName()," dist[mm]= %8.1f "%(dist*1000.)," L= %8.1f "%((x_max-x_min)*1000.)," (x_min,x_max) =  (%8.1f,%8.1f) "%(x_min*1000.,x_max*1000.)

print "=================================================================="
print "Stop."
	
	
	
