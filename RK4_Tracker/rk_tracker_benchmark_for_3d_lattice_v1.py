#-----------------------------------------------------
# Track ORBIT two particles through the TEAPOT lattice
# and through the 3D lattice with Field Source class
# instances that should be equivalent to the hard edge
# quads and kickers.
# This variant does not include any kickers fields.
# Tracking will be exactly along z-axis of the TEAPOT 
# lattice.
# User can play with number of steps and the accuracy eps parameter
# in the Runge-Kutta tracker of 4th order.
#
# Author: Andrei Shishlo
#-----------------------------------------------------

import sys
import math
import random 

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice import ThickKick
from orbit.py_linac.lattice import Drift, Quad, Bend, MarkerLinacNode

from bunch import Bunch
from bunch import BunchTwissAnalysis

from trackerrk4 import RungeKuttaTracker

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import KVDist2D, KVDist3D
from orbit.bunch_generators import GaussDist2D, GaussDist3D
from orbit.bunch_generators import WaterBagDist2D, WaterBagDist3D

from orbit_utils import Function

#--------------------------------------------------------
# Classes for 3D tracking
#--------------------------------------------------------
from orbit_utils import Matrix, PhaseVector
from orbit_utils import FieldSourceContainer
from orbit_utils import field_sources
from field_sources import QuadFieldSource
from field_sources import DipoleFieldSource
from field_sources import MagnetFieldSourceGrid3D

from trackerrk4 import RungeKuttaTracker
from trackerrk4 import PyExternalEffects

random.seed(100)

seq_names = ["TEST_SEQ"]
#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- the XML file name with the structure
xml_file_name = "test_3d_lattice.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(seq_names,xml_file_name)

print "Linac lattice is ready. L=",accLattice.getLength()
print "========================================================"

kickers = accLattice.getNodesOfClass(ThickKick)
n_kicker_parts = 100
print "debug n_kicker_parts = ",n_kicker_parts
for kicker in kickers:
	kicker.setFieldBx(0.0)
	kicker.setFieldBy(0.0)
	kicker.setnParts(n_kicker_parts)
	print "kicker=",kicker.getName()," field (X,Y) [T] = (%6.5f,%6.5f) "%(kicker.getFieldBx(),kicker.getFieldBy())
	
print "========================================================"
quads = accLattice.getNodesOfClass(Quad)

#--- The last quad will have zero field.
#--- we will use it in the other examples 
#--- to provide tilted exit plane.
#--- It is better to have exit from the 3D region
#--- in the free space.
quads[len(quads) - 1].setParam("dB/dr",0.)

for quad in quads:
	quad.setnParts(10)
	print "quad  =",quad.getName()," gradient[T/m] = %+6.5f "%quad.getParam("dB/dr")
print "========================================================"


bunch_ini = Bunch()
eKin = 1.3
bunch_ini.getSyncParticle().kinEnergy(eKin)

#---- There are 2 macro-particles added to the bunch
#---- User can change coordinates to see the difference 
#---- between TEAPOT and 3D RK4 tracking.
(x,xp,y,yp,z,dE) = (0.,0.,0.,0.,0.,0.)
bunch_ini.addParticle(x,xp,y,yp,z,dE)

(x,xp,y,yp,z,dE) = (0.05,0.0,0.05,0.0,0.,0.001)
bunch_ini.addParticle(x,xp,y,yp,z,dE)

bunch_tmp = Bunch()
bunch_ini.copyBunchTo(bunch_tmp)

start_ind = accLattice.getNodeIndex(quads[0])
stop_ind  = accLattice.getNodeIndex(quads[len(quads)-1])

accLattice.trackDesignBunch(bunch_tmp,None,None, index_start = start_ind, index_stop = stop_ind)

traj_x_function  = Function()
traj_xp_function = Function()

#----- if you want more beam size points along the trajectory - change pos_step here
paramsDict = {"old_pos":-1.,"count":0,"pos_step":0.1,"path_length":0.}
actionContainer = AccActionsContainer("TEAPOT Bunch Tracking")

#----- start of the 1st quad
pos_start = - quads[0].getLength()/2

def action_account(paramsDict):
	node = paramsDict["node"]
	name = node.getName()		
	bunchI = paramsDict["bunch"]
	pos = paramsDict["path_length"]
	if(paramsDict["old_pos"] == pos): return
	if(paramsDict["old_pos"] + paramsDict["pos_step"] > pos): return
	paramsDict["old_pos"] = pos
	paramsDict["count"] += 1
	traj_x_function.add(pos+pos_start,bunchI.x(0)*1000.)
	traj_xp_function.add(pos+pos_start,bunchI.xp(0)*1000.)
	st = "i= %3d "%paramsDict["count"] + " pos= %8.3f "%(pos+pos_start)
	(x,xp,y,yp,z,dE) = (bunchI.x(0)*1000.,bunchI.xp(0)*1000.,bunchI.y(0)*1000.,bunchI.yp(0)*1000.,bunchI.z(0)*1000.,bunchI.dE(0)*1000.)
	st += " (x,xp,y,yp,z,dE) = ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f )"%(x,xp,y,yp,z,dE)
	(x,xp,y,yp,z,dE) = (bunchI.x(1)*1000.,bunchI.xp(1)*1000.,bunchI.y(1)*1000.,bunchI.yp(1)*1000.,bunchI.z(1)*1000.,bunchI.dE(1)*1000.)
	st += " ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f )"%(x,xp,y,yp,z,dE)
	print st

actionContainer.addAction(action_account, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_account, AccActionsContainer.BODY)
actionContainer.addAction(action_account, AccActionsContainer.EXIT)

print "============================================================================================"
accLattice.trackBunch(bunch_tmp,paramsDict = paramsDict, actionContainer = actionContainer, index_start = start_ind, index_stop = stop_ind)
print "============================================================================================"

print "============================ final TEAPOT bunch ==========================="
st = ""
(x,xp,y,yp,z,dE) = (bunch_tmp.x(0)*1000.,bunch_tmp.xp(0)*1000.,bunch_tmp.y(0)*1000.,bunch_tmp.yp(0)*1000.,bunch_tmp.z(0)*1000.,bunch_tmp.dE(0)*1000.)
st += " (x,xp,y,yp,z,dE) = ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f)"%(x,xp,y,yp,z,dE)
(x,xp,y,yp,z,dE) = (bunch_tmp.x(1)*1000.,bunch_tmp.xp(1)*1000.,bunch_tmp.y(1)*1000.,bunch_tmp.yp(1)*1000.,bunch_tmp.z(1)*1000.,bunch_tmp.dE(1)*1000.)
st += " ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f)"%(x,xp,y,yp,z,dE)
print st
print "==========================================================================="

#----- 3D quads
qv01_3d = QuadFieldSource()
qh02_3d = QuadFieldSource()
qv03_3d = QuadFieldSource()
qh04_3d = QuadFieldSource()

qv01_3d.length(quads[0].getLength())
qh02_3d.length(quads[1].getLength())
qv03_3d.length(quads[2].getLength())
qh04_3d.length(quads[3].getLength())

qv01_3d.gradient(accLattice.getNodeForName("QV01").getParam("dB/dr"))
qh02_3d.gradient(accLattice.getNodeForName("QH02").getParam("dB/dr"))
qv03_3d.gradient(accLattice.getNodeForName("QV03").getParam("dB/dr"))
qh04_3d.gradient(accLattice.getNodeForName("QH04").getParam("dB/dr"))

qv01_pos_end   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QV01")][1]
qh02_pos_end   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QH02")][1]
qv03_pos_end   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QV03")][1]
qh04_pos_end   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QH04")][1]


qv01_pos_start   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QV01")][0]
qh02_pos_start   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QH02")][0]
qv03_pos_start   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QV03")][0]
qh04_pos_start   = accLattice.getNodePositionsDict()[accLattice.getNodeForName("QH04")][0]

qv01_pos_center   = (qv01_pos_start + qv01_pos_end)/2.0

#---- qv01_3d    coordinate transformation from the lattice system
transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()
dist_from_qv01_center = 0.
transfCoordsMatrix.set(2,3,-dist_from_qv01_center)
qv01_3d.transormfMatrix(transfCoordsMatrix)

#---- qh02_3d    coordinate transformation from the lattice system
transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()
dist_from_qv01_center = qh02_pos_end - qv01_pos_end
transfCoordsMatrix.set(2,3,-dist_from_qv01_center)
qh02_3d.transormfMatrix(transfCoordsMatrix)

#---- qv03_3d    coordinate transformation from the lattice system
transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()
dist_from_qv01_center = qv03_pos_end - qv01_pos_end
transfCoordsMatrix.set(2,3,-dist_from_qv01_center)
qv03_3d.transormfMatrix(transfCoordsMatrix)

#---- qh04_3d    coordinate transformation from the lattice system
transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()
dist_from_qv01_center = qh04_pos_end - qv01_pos_end
transfCoordsMatrix.set(2,3,-dist_from_qv01_center)
qh04_3d.transormfMatrix(transfCoordsMatrix)

kickers = accLattice.getNodesOfClass(ThickKick)
kickers_3d = []

sizeX = 2.0
sizeY = 2.0
for kicker in kickers:
	kicker_3d = DipoleFieldSource()
	kicker_pos = (accLattice.getNodePositionsDict()[kicker][0] + accLattice.getNodePositionsDict()[kicker][1])/2
	kicker_length = kicker.getLength()
	kicker_3d.sizesXYZ(sizeX,sizeY,kicker_length)
	kicker_3d.fieldsXYZ(0.,-kicker.getFieldBy(),0.)
	transfCoordsMatrix = Matrix(4,4)
	transfCoordsMatrix.unit()
	transfCoordsMatrix.set(2,3,(qv01_pos_center - kicker_pos))
	#print "debug (qv01_pos_center - kicker_pos)=",(qv01_pos_center - kicker_pos)," kicker_length=",kicker_length
	kicker_3d.transormfMatrix(transfCoordsMatrix)
	print "kicker=",kicker.getName()," field (X,Y,Z) = ",kicker_3d.fieldsXYZ()
	kickers_3d.append(kicker_3d)


#----- field source with all Field Sources from magnets
field_container = FieldSourceContainer()

field_container.addFieldSource(qv01_3d)
field_container.addFieldSource(qh02_3d)
field_container.addFieldSource(qv03_3d)
field_container.addFieldSource(qh04_3d)

for kicker_3d in kickers_3d:
	field_container.addFieldSource(kicker_3d)
#-----------------------------------------------------

#---- Now let's define RK4 tracker and entrance and exit planes
#---- as (r_vector - r0_vector)*n_vector = 0
#---- where n_vector is a normal vector to the plane
#---- For all planes these normal vectors should be directed to the 
#---- outer space.
length = qh04_pos_end - qv01_pos_start
tracker = RungeKuttaTracker(length)

#------------------------------------------------------------------------
#----- This is accuracy parameter in meters for the synchronous particle
#----- tracking using RK4. In normal circumstances it defines the
#----- size of time steps for integration, but here it is useless
#----- because the field on the axis is equal to zero.
tracker.spatialEps(0.001)

#----- palnes setup
tracker.entrancePlane(0.,0.,-1.,(qv01_pos_start - qv01_pos_center))
tracker.exitPlane(0.,0.,1.,-(qh04_pos_end - qv01_pos_center))
print "Tracker Entrance plane (a,b,c,d)=",tracker.entrancePlane()
print "Tracker Exit     plane (a,b,c,d)=",tracker.exitPlane()

# we have to specify an initial guess for number of time steps
tracker.stepsNumber(10000)

bunch_tmp = Bunch()
bunch_ini.copyBunchTo(bunch_tmp)

n_space_parts = 1
(a_entr,b_entr,c_entr,d_entr) = tracker.entrancePlane() 
(a_exit,b_exit,c_exit,d_exit) = tracker.exitPlane()
d_start = - d_entr
d_stop  = d_exit
space_step = (d_stop + (-d_start))/n_space_parts

print "====================== Initial bunch before 3D tracking=============================="
st = " pos = %8.5f "%d_entr
(x,xp,y,yp,z,dE) = (bunch_tmp.x(0)*1000.,bunch_tmp.xp(0)*1000.,bunch_tmp.y(0)*1000.,bunch_tmp.yp(0)*1000.,bunch_tmp.z(0)*1000.,bunch_tmp.dE(0)*1000.)
st += " (x,xp,y,yp,z,dE) = ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f)"%(x,xp,y,yp,z,dE)
(x,xp,y,yp,z,dE) = (bunch_tmp.x(1)*1000.,bunch_tmp.xp(1)*1000.,bunch_tmp.y(1)*1000.,bunch_tmp.yp(1)*1000.,bunch_tmp.z(1)*1000.,bunch_tmp.dE(1)*1000.)
st += " ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f)"%(x,xp,y,yp,z,dE)
print st
print "======================================================================================"

print "debug d_start=",d_start," d_stop=",d_stop," space_step=",space_step," n_space_parts=",n_space_parts
for space_ind in range(n_space_parts):
	d_plane_start = -(d_start + space_step*space_ind)
	d_plane_stop  = -d_plane_start + space_step
	#print "debug d_plane_start=",d_plane_start," d_plane_stop=",d_plane_stop
	tracker.entrancePlane(a_entr,b_entr,c_entr,d_plane_start)
	tracker.exitPlane(a_exit,b_exit,c_exit,d_plane_stop)
	bunch_synch_part = Bunch()
	bunch_tmp.copyEmptyBunchTo(bunch_synch_part)	
	tracker.trackBunch(bunch_synch_part,field_container)
	synch_part_pvector = bunch_synch_part.getSyncParticle().pVector()
	synch_part_rvector = bunch_synch_part.getSyncParticle().rVector()
	momentum = math.sqrt(synch_part_pvector[0]**2 + synch_part_pvector[1]**2 + synch_part_pvector[2]**2)
	norm_final_v = [synch_part_pvector[0]/momentum,synch_part_pvector[1]/momentum,synch_part_pvector[2]/momentum]
	d_parameter = norm_final_v[0]*synch_part_rvector[0] + norm_final_v[1]*synch_part_rvector[1] + norm_final_v[2]*synch_part_rvector[2]	
	tracker.exitPlane(norm_final_v[0],norm_final_v[1],norm_final_v[2],-d_parameter)
	tracker.trackBunch(bunch_tmp,field_container)
	#----------------------------------------------------------------------------
	st = "i = %3d "%space_ind + " pos = %8.5f "%d_parameter
	(x,xp,y,yp,z,dE) = (bunch_tmp.x(0)*1000.,bunch_tmp.xp(0)*1000.,bunch_tmp.y(0)*1000.,bunch_tmp.yp(0)*1000.,bunch_tmp.z(0)*1000.,bunch_tmp.dE(0)*1000.)
	st += " (x,xp,y,yp,z,dE) = ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f)"%(x,xp,y,yp,z,dE)
	(x,xp,y,yp,z,dE) = (bunch_tmp.x(1)*1000.,bunch_tmp.xp(1)*1000.,bunch_tmp.y(1)*1000.,bunch_tmp.yp(1)*1000.,bunch_tmp.z(1)*1000.,bunch_tmp.dE(1)*1000.)
	st += " ( %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f)"%(x,xp,y,yp,z,dE)
	print st


print "======================================================================================"
print "Stop."

