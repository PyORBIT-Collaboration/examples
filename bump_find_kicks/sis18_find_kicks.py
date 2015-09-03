##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice and add injection 
##############################################################
import sys
import pickle
import math
import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize
from pylab import *


from bunch import Bunch
from bunch import BunchTwissAnalysis

# lattice, teapot class
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch, BunchTwissAnalysis
from orbit.errors import AddErrorNode


from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotStatLatsNodeSet, addTeapotMomentsNodeSet
from orbit.diagnostics import BPMSignal
from orbit.diagnostics import TeapotBPMSignalNode

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import  KVDist2D

from scipy.optimize import leastsq

from orbit.bumps import simpleBump, closed_orbit_bumps
#======================plot===============
def plot_info(lattice):
	xAvg = []
	yAvg = []
	s = []
	#============Set BPMs into lattice===============
	length = round(lattice.getLength()/12,4)
	for i in range(12):
		name = TeapotBPMSignalNode("BPM")
		addTeapotDiagnosticsNode(lattice, i*length, name)

	xAvg = []
	yAvg = []
	s = []

	lattice.trackBunch(bunch)
	lattice.trackBunch(bunch)

	#============Get the BPM signal===============
	nodes = lattice.getNodes()
	for i in range(1):
		lattice.trackBunch(bunch)
		for node in nodes:
				if node.getType()=="BPMSignal":
					xAvg.append(node.getSignal()[0])
					yAvg.append(node.getSignal()[1])
					s.append(lattice.getNodePositionsDict()[node][0]+lattice.getLength()*(i))
					#print node.getSignal()[0]*1e3, node.getSignal()[1]*1e3
					
	fig, (ax1, ax2) = plt.subplots(2, sharex='col')

	ax1.plot(array(s),array(xAvg)*1e3, '-rx')
	ax2.plot(array(s),array(yAvg)*1e3, '-bx')

	ax1.set_ylim(-100,100)
	ax2.set_ylim(-100,100)

	ax1.set_ylabel('x [mm]')
	ax2.set_ylabel('xp [mm]')
	ax2.set_xlabel('s [m]')
#======================plot===============

#============add particles to beam===============
def set_bunch(off,lattice,bunch):
	bunch.deleteAllParticles()
	x_offset = off[0]
	xp_offset = off[1]
	y_offset = off[2]
	yp_offset = off[3]

	matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)
	(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
	(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
	(arrDispX,arrDispPrimeX) = matrix_lattice.getRingDispersionDataX()
	(arrDispY,arrDispPrimeY) = matrix_lattice.getRingDispersionDataY()

	alpax= arrPosAlphaX[0][1]
	betax = arrPosBetaX[0][1]

	alpay= arrPosAlphaY[0][1]
	betay = arrPosBetaY[0][1]

	# ini. 2D match distribution
	twissX = TwissContainer(alpha = alpax, beta = betax, emittance = emittance_x/2)
	twissY = TwissContainer(alpha = alpay, beta = betay, emittance = emittance_y/2)
	dist = KVDist2D(twissX,twissY)


	for i in range(NPIC):
		(x,xp,y,yp) = dist.getCoordinates()
		bunch.addParticle(x+x_offset,xp+xp_offset,y+y_offset,yp+yp_offset,0.0,0.0)
#============add particles to beam===============

#=====set up bunch stuff============
bunch = Bunch()
print "Read Bunch."
NPIC = 500 

A = 1
Z = 1
energy = 11.4e-3*A # Gev
intensity = 1.4e10
emittance_x = 50e-6
emittance_y = 50e-6
bunch.mass(0.93827231)
bunch.getSyncParticle().kinEnergy(energy)
#=====set up bunch stuff============

	
print "Generate Lattice."
lattice = TEAPOT_Lattice("no_sc_lattice")

case = "inj"
#case = "mid"

if case == "inj":
	lattice.readMAD(".sis18_inj.lat","SIS18") # lattice start at injection point
if case == "mid":
	lattice.readMAD(".sis18_inj.lat","SIS18_MID") # lattice start at injection point plus half circumference 


# get the lattice function for the kick and bpms, jet only horizontal 
find = close_orbit_bumps()

xc0 = 70e-3
xcs0 = 0.0e-3 # or 7.0e-3

# set variable, list [[node.name], [variable, start value], ...]
kicker1 = [["S11MB1"], ["kx",-0.00111]]
kicker2 = [["S12MB2"], ["kx",0.00]]
kicker3 = [["S01MB3"], ["kx",0.0333]]
kicker4 = [["S03MB4"], ["kx", -0.00111]]
variables = [kicker1, kicker2, kicker3, kicker4]

# set constraints, list [[node.name], [coordinate (x or xp), value], ...]
con1 = [["S12DX5H"], ["x", xc0]]
con2 = [["S11MB1"], ["x", 0.0]]
con3 = [["S03MB4"], ["x", 0.0]]
con4 = [["S12DX5H"], ["xp", xcs0]]
constraints = [con1, con2, con3, con4]


find.lattice_function(lattice, bunch, variables, constraints)

kick_nelder = find.find_kicks(method="Simplex")
kick_leastsq = find.find_kicks(method="LeastSq")

# analytic solution for four bumps like the injection for the GSI SIS18,
if xcs0 == 0.0:
	kick_ana = [-0.002353486972582688, 0.011463597538482626, 0.00421523569549529, -0.006945791341904927] # xc0 = 70e-3, xcs0 = 0.0e-3
if xcs0 == 7.0e-3:
	kick_ana = [-0.008555530928460921, 0.01670459776629772, 0.002272837942467913, -0.015418194930452534] # xc0 = 70e-3, xcs0 = 7.0e-3

print "solution:", kick_ana
print "python nelder", kick_nelder
print "python leastsq", kick_leastsq


kick = kick_ana
nodes = lattice.getNodes()
for node in nodes:
	if node.getName() == "S11MB1":
		node.setParam("kx",kick[0])
	if node.getName() == "S12MB2":
		node.setParam("kx",kick[1])
	if node.getName() == "S01MB3":
		node.setParam("kx",kick[2])
	if node.getName() == "S03MB4":
		node.setParam("kx",kick[3])

if case == "inj":
	pos = [xc0,xcs0,0.0,0.0]  # for SIS18 lattice
if case =="mid":
	pos = [0.0,0.0,0.0,0.0]  # for SIS18_MID lattice

set_bunch(pos,lattice,bunch)
lattice.trackBunch(bunch)
plot_info(lattice)
show()
exit()

