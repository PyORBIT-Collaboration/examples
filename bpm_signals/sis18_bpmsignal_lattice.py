##############################################################
# This script reads the input MAD file with lattice information including monitors,
# creates the TEAPOT lattice plus get the information 
##############################################################
import sys
import math
import numpy as np
from pylab import *


from bunch import Bunch
from bunch import BunchTwissAnalysis

# lattice, teapot class
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice


from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.diagnostics import BPMSignal
from orbit.diagnostics import TeapotBPMSignalNode

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import KVDist2D


#=====set up bunch stuff============

bunch = Bunch()
print "Read Bunch."
NPIC = 1000 

A = 1
Z = 1
energy = 11.4e-3*A # Gev
intensity = 1.4e10
		
		
emittance_x = 50e-6
emittance_y = 50e-6
		
bunch.mass(0.93827231)
bunch.getSyncParticle().kinEnergy(energy)
	
print "Generate Lattice."
lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD("./lattice/sis18_inj.lat","SIS18_BPM_MID")
offset = 0.#70e-3



matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
(arrDispX,arrDispPrimeX) = matrix_lattice.getRingDispersionDataX()
(arrDispY,arrDispPrimeY) = matrix_lattice.getRingDispersionDataY()

alpax= arrPosAlphaX[0][1]
betax = arrPosBetaX[0][1]

alpay= arrPosAlphaY[0][1]
betay = arrPosBetaY[0][1]

#============add particles to beam===============
# ini. 3D match distribution
twissX = TwissContainer(alpha = alpax, beta = betax, emittance = emittance_x/2)
twissY = TwissContainer(alpha = alpay, beta = betay, emittance = emittance_y/2)
dist = KVDist2D(twissX,twissY)


for i in range(NPIC):
	(x,xp,y,yp) = dist.getCoordinates()
	bunch.addParticle(x+offset,xp,y,yp,0.0,0.0)
#============add particles to beam===============



xAvg = []
yAvg = []
s = []


#============Get BPMs signal===============
nodes = lattice.getNodes()
for i in range(1):
	lattice.trackBunch(bunch)
	for node in nodes:
			if node.getType()=="monitor":
				xAvg.append(node.getParam("xAvg"))
				yAvg.append(node.getParam("yAvg"))
				s.append(lattice.getNodePositionsDict()[node][0]+lattice.getLength()*(i))

#======================plot===============
fig, (ax1, ax2) = plt.subplots(2, sharex='col')

ax1.plot(array(s),array(xAvg)*1e3, '-rx')
ax2.plot(array(s),array(yAvg)*1e3, '-bx')

ax1.set_ylim(-100,100)
ax2.set_ylim(-100,100)

ax1.set_ylabel('x [mm]')
ax2.set_ylabel('y [mm]')
ax2.set_xlabel('s [m]')


show()

print "done"	


