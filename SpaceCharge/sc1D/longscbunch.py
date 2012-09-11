
import math
import sys

from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from spacecharge import LSpaceChargeCalc
from orbit.space_charge.sc1d import addLongitudinalSpaceChargeNode	
from orbit.space_charge.sc1d import SC1D_AccNode
	

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
runName = "Benchmark_Collimator"

total_macroSize=1.0e+10
b.mass(0.93827231)
energy = 1.0 #Gev
bunch_orbit_to_pyorbit(lattice.getLength(), energy, "Bm_KV_Uniform_10000",b)
b.getSyncParticle().kinEnergy(energy)
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize/nParticlesGlobal)

#=====track bunch through Foil============

b_a = 10.0/3.0
length = 248.0
nMacrosMin = 1
useSpaceCharge = 1
nBins= 128     #number of longitudinal slices in the 1D space charge solver
position = 0.0

#SNS Longitudinal Impedance tables. EKicker impedance from private communication
# with J.G. Wang. Seems to be for 7 of the 14 kickers (not sure why).  
# Impedance in Ohms/n. Kicker and RF impedances are inductive with real part positive and imaginary is negative by Chao definition. 

ZL_EKicker = [complex(42., -182),
	 complex(35, -101.5), 
	 complex(30.3333, -74.6667), 
	 complex(31.5, -66.5),
	 complex(32.2,-57.4), 
	 complex(31.5, -51.333), 
	 complex(31, -49), 
	 complex(31.5, -46.375),
	 complex(31.8889, -43.556), 
	 complex(32.9, -40.6), 
	 complex(32.7273, -38.18),
	 complex(32.25, -35.58),
	 complex(34.46, -32.846),
	 complex(35, -30.5),
	 complex(35.4667, -28.),
	 complex(36.75, -25.81),
	 complex(36.647, -23.88),
	 complex(36.944, -21.1667),
	 complex(36.474, -20.263),
	 complex(36.4, -18.55),
	 complex(35.333, -17),
	 complex(35, -14.95),
	 complex(33.478, -13.69),
	 complex(32.375, -11.67),
	 complex(30.8, -10.08),
	 complex(29.615, -8.077),
	 complex(28.519, -6.74),
	 complex(27.5, -5),
	 complex(26.552, -4.103),
	 complex(25.433, -3.266),
	 complex(24.3871, -2.7),
	 complex(23.40625, -2.18)]

ZL_RF = [complex(0.0, 0.0),
		 complex(0.750, 0.0), 
		 complex(0.333,0.0),
		 complex(0.250, 0.0),
		 complex(0.200, 0.0),
		 complex(0.167, 0.0),
		 complex(3.214, 0.0),
		 complex(0.188, 0.0),
		 complex(0.167, 0.0),
		 complex(0.150, 0.0),
		 complex(1.000, 0.0),
		 complex(0.125, 0.0),
		 complex(0.115, 0.0),
		 complex(0.143, 0.0),
		 complex(0.333, 0.0),
		 complex(0.313, 0.0),
		 complex(0.294, 0.0),
		 complex(0.278, 0.0),
		 complex(0.263, 0.0),
		 complex(0.250, 0.0),
		 complex(0.714, 0.0),
		 complex(0.682, 0.0),
		 complex(0.652, 0.0),
		 complex(0.625, 0.0),
		 complex(0.600, 0.0),
		 complex(0.577, 0.0),
		 complex(0.536, 0.0),
		 complex(0.536, 0.0),
		 complex(0.517, 0.0),
		 complex(0.500, 0.0),
		 complex(0.484, 0.0),
		 complex(0.469, 0.0)]
			 

Z = []
for i in range(0,32):
	zk = ZL_EKicker[i]
	zrf = ZL_RF[i]
	#Multiply by 10 to make the effect bigger for benchmark purpose.
	zreal = 10.0*(zk.real/1.75 + zrf.real)
	zimag = 10.0*(zk.imag/1.75 + zrf.imag)
	#zreal = 10.0*(zk.real/1.75)
	#zimag = 10.0*(zk.imag/1.75)
	Z.append(complex(zreal, zimag))

sc1Dnode = SC1D_AccNode(b_a,length,nMacrosMin,useSpaceCharge,nBins)
sc1Dnode.assignImpedance(Z);


addLongitudinalSpaceChargeNode(lattice, position, sc1Dnode)

print "===========Lattice modified ======================================="
print "New Lattice=",lattice.getName()," length [m] =",lattice.getLength()," nodes=",len(lattice.getNodes())

paramsDict = {}
paramsDict["bunch"]= b
b.dumpBunch("bunch_init.dat")
print "here about to track"
sc1Dnode.trackBunch(b)
print "tracking done"
#for i in range (0:10):
#	lattice.trackBunch(b, paramsDict)
b.dumpBunch("bunch_final.dat")
bunch_pyorbit_to_orbit(lattice.getLength(), b, "pybunch_final.dat")
print "Stop."



