#!/usr/bin/env python

#--------------------------------------------------------
# Test of the Linac Aperture for quads
# The example includes quads with overlapping fields and the lost bunch.
# The user can comment out the parts with these quads or the lost bunch 
# in the paramsDict dictionary.
#--------------------------------------------------------

import math
import sys
import os

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory
from orbit.py_linac.overlapping_fields import SNS_MEBT_OverlappingQuadsSubst

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import GetLostDistributionArr
from orbit.py_linac.lattice_modifications import AddScrapersAperturesToLattice
from orbit.py_linac.lattice_modifications import AddMEBTChopperPlatesAperturesToSNS_Lattice

from orbit.lattice import AccLattice, AccActionsContainer

from orbit.py_linac.lattice import BaseLinacNode

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D
from orbit.bunch_generators import TwissAnalysis

from bunch import Bunch
from bunch import BunchTwissAnalysis


class LostBunchDumpNode(BaseLinacNode):
	"""
	The class LostBunchDumpNode writes the information from the lostbunch
	to the file.
	"""
	def __init__(self, name = "LostBunchDump"):
		BaseLinacNode.__init__(self,name)
		self.file_name = "lost_bunch.dat"
	
	def setFileName(self,file_name):
		self.file_name = file_name
		
	def getFileName(self):
		return self.file_name
	
	def track(self, paramsDict):
		if(paramsDict.has_key("lostbunch")):
			lostbunch = paramsDict["lostbunch"]
			lostbunch.dumpBunch(self.file_name)

	def trackDesign(self, paramsDict):
		"""
		This method does nothing for this class.
		"""
		pass

#-------------------------------------------
# START of Script
#-------------------------------------------

comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
size = orbit_mpi.MPI_Comm_size(comm)		
data_type = mpi_datatype.MPI_DOUBLE		
main_rank = 0		


names = ["MEBT",]

py_orbit_sns_home = "../"

#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- the XML file name with the structure
xml_file_name = py_orbit_sns_home+"sns_linac_xml/sns_linac.xml"
xml_file_name = py_orbit_sns_home+"sns_linac_xml/sns_linac_with_aprt.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print "1. Linac lattice is ready. L=",accLattice.getLength()

#--------------------------------------------------------
# Replace overlapping quads with the special nodes
#--------------------------------------------------------

(ovrlp_quads1,ovrlp_quads2) = SNS_MEBT_OverlappingQuadsSubst(accLattice)

print "2. Linac lattice is ready. L=",accLattice.getLength()

node_pos_dict = accLattice.getNodePositionsDict()
nodes = accLattice.getNodes()
for node in nodes:
	(posBefore, posAfter) = node_pos_dict[node]
	print "%35s      start stop L = %10.4f %10.4f    %10.4f "%(node.getName(),posBefore,posAfter,(posAfter-posBefore))

#----------------------------------------------------------
#      1. add Aperture nodes to the quads in the linac lattice
#      2. add Aperture nodes to the MEBT chopper entrance/exit plates 
#         in the linac lattice
#      3. add Aperture nodes to the MEBT scrapers (H and V)
#      4. add a LostBunchDumpNode to one of the chopper plate
#----------------------------------------------------------

print "===== Aperture Nodes ======="
aprtNodes = Add_quad_apertures_to_lattice(accLattice)
aprtNodes = AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice,aprtNodes)
x_size = 0.042
y_size = 0.042
aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:H_SCRP",x_size,y_size,aprtNodes)

x_size = 0.042
y_size = 0.042
aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:V_SCRP",x_size,y_size,aprtNodes)

for node in aprtNodes:
	print "aprt=",node.getName()," pos =",node.getPosition()

#---- set the Lost Bunch Dump Node to see the one of the apertures effect 
aprt_node_1 = None
for node in aprtNodes:	
	if(node.getName().find("MEBT:ChpPlt:Entr") >= 0):
		aprt_node_1 = node
		break
if(aprt_node_1 != None):
	dumpNode = LostBunchDumpNode()
	dumpNode.setFileName("lost_bunch_chpplt.dat")
	dumpNode.setSequence(aprt_node_1.getSequence())
	aprt_node_1.addChildNode(dumpNode,node.EXIT)

#-----------------------------------------------------------
#    Bunch Generation
#-----------------------------------------------------------
bunch = Bunch()
#set H- mass
bunch.mass(0.9382723 + 2*0.000511)
bunch.charge(-1.0)
bunch.getSyncParticle().kinEnergy(0.0025)

# electron charge in SI
si_e_charge = 1.6021773e-19
frequency = 402.5e+6
beam_current = 0.038

N_particles = 10000
macrosize = (beam_current/frequency)
macrosize /= (math.fabs(bunch.charge())*si_e_charge)
macrosize /= N_particles

(alphaX,betaX,emittX) = (-1.39, 0.126, 3.67*1.0e-6)
(alphaY,betaY,emittY) = ( 2.92, 0.281, 3.74*1.0e-6)
(alphaZ,betaZ,emittZ) = ( 0.0 , 117.0, 0.0166*1.0e-6)

#---- we increase the emittances to see apertures effects
emittX *= 5.
emittY *= 10.

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

distributor = WaterBagDist3D(twissX,twissY,twissZ)

for ind in range(N_particles):
	(x,xp,y,yp,z,dE) = distributor.getCoordinates()
	(x,xp,y,yp,z,dE) = orbit_mpi.MPI_Bcast((x,xp,y,yp,z,dE),data_type,main_rank,comm)
	if(ind%size == rank):
		bunch.addParticle(x,xp,y,yp,z,dE)
		
nParticlesGlobal = bunch.getSizeGlobal()
if(rank == 0):
	print "total number of particles =",nParticlesGlobal
bunch.macroSize(macrosize)

#set up design
accLattice.trackDesignBunch(bunch)

paramsDict = {}
lost_parts_bunch = Bunch()
paramsDict["lostbunch"] = lost_parts_bunch

actionContainer = AccActionsContainer("Test Design Bunch Tracking")

twiss_analysis = BunchTwissAnalysis()

def action_entrance(paramsDict):
	if(isinstance(paramsDict["parentNode"],AccLattice)):
		node = paramsDict["node"]
		pos = paramsDict["path_length"]
		bunch = paramsDict["bunch"]
		twiss_analysis.analyzeBunch(bunch)
		x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
		y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
		z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
		s = " %35s  %4.5f    %5.3f  %5.3f   %5.3f     %10.6f   %8d "%(node.getName(),pos,x_rms,y_rms,z_rms,eKin,bunch.getSizeGlobal())
		print s	

actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
accLattice.trackBunch(bunch, paramsDict = paramsDict, actionContainer = actionContainer)

lost_parts_bunch.dumpBunch("lostbunch.dat")

aprtNodes_loss_arr = GetLostDistributionArr(aprtNodes,lost_parts_bunch)
total_loss = 0.
for [aprtNode,loss] in aprtNodes_loss_arr:
	print "aprt. node= %30s "%aprtNode.getName()," pos= %9.3f "%aprtNode.getPosition()," loss= %6.0f "%loss
	total_loss += loss
print "Total loss=",total_loss
