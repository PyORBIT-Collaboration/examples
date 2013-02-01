
import math
import sys

from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.space_charge.directforce2p5d import directforceAccNodes, directforceLatticeModifications
from spacecharge import SpaceChargeForceCalc2p5D

print "Start."

#------------------------------------------
#
#          Now let's make a lattice
#------------------------------------------

def getLattice(lattice_length, n_parts):
	elem = teapot.DriftTEAPOT("a drift")
	elem.setLength(lattice_length)
	elem.setnParts(n_parts)	
	teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
	teapot_lattice.addNode(elem)
	teapot_lattice.initialize()
	return teapot_lattice

lattice_length = 248.0    # the length of the drift
n_parts = 1  # number of parts on what the drift will be chopped, or the number of SC nodes
lattice = getLattice(lattice_length,n_parts)

#------------------------------
#Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "Benchmark_SpaceCharge"

total_macroSize=1.0e+14
b.mass(0.93827231)
energy = 1.0 #Gev
bunch_orbit_to_pyorbit(lattice.getLength(), energy, "Bm_KV_Uniform_10",b)
b.dumpBunch("bunch_init.dat")
b.getSyncParticle().kinEnergy(energy)
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize/nParticlesGlobal)
print total_macroSize/nParticlesGlobal
#-----------------------------------
# Add Direct Force Space Charge node
#-----------------------------------

nMacrosMin = 1
sc_path_length_min = 0.00000001
sizeX = 9   #number of grid points in horizontal direction
sizeY = 9   #number of grid points in vertical direction
sizeZ = 1   #number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeForceCalc2p5D(sizeX,sizeY,sizeZ)
calc2p5d.trackBunch(b, lattice_length)
#directforceLatticeModifications.setDirectForce2p5DAccNodes(teapot_latt, sc_path_length_min, calc2p5d)



#Now let's calculate the exact force grid and compare with the force grid from the code
rhogrid = calc2p5d.getRhoGrid()
file_out = open("exactforce.dat","w")
for i in range(rhogrid.getSizeX()):
	for j in range(rhogrid.getSizeY()):
			x0=rhogrid.getGridX(i)
			y0=rhogrid.getGridY(j)
			forcex = 0
			forcey = 0
			for k in range(rhogrid.getSizeY()):
				for l in range(rhogrid.getSizeX()):
					xlocal = rhogrid.getGridX(k)
					ylocal = rhogrid.getGridY(l)
					x = xlocal - x0
					y = ylocal - y0
					r = math.sqrt((x0 - xlocal)**2 + (y0 - ylocal)**2)
					theta=math.atan2(y,x)
					if(r != 0.):
						forcex += -rhogrid.getValueOnGrid(k,l) * 1/(r) * math.cos(theta)
						forcey += -rhogrid.getValueOnGrid(k,l) * 1/(r) * math.sin(theta)
			file_out.write(str(i) + " " + str(j) + " " + str(forcex) + " " + str(forcey) + "\n")
file_out.close()

#Get the force grid from the space charge calculator
file_out = open("pyforce.dat","w")
forcegridx = calc2p5d.getForceGridX()
forcegridy = calc2p5d.getForceGridY()
for i in range(forcegridx.getSizeX()):
	for j in range(forcegridx.getSizeY()):
		file_out.write(str(i) + " " + str(j) + " " + str(forcegridx.getValueOnGrid(i,j)) + " " + str(forcegridy.getValueOnGrid(i,j)) + "\n")		
file_out.close()

print "===========Lattice modified ======================================="
print "New Lattice=",lattice.getName()," length [m] =",lattice.getLength()," nodes=",len(lattice.getNodes())

#------------------------------------

paramsDict = {}
paramsDict["bunch"]= b
print "tracking done"
b.dumpBunch("pybunch.dat")
bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch.dat")
print "Stop."



