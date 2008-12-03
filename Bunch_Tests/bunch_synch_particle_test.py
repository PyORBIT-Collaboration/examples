import sys
import math
from bunch import Bunch
from bunch import SyncParticle

import orbit_mpi

#-----------------------------------------------------
#Memory leak test - it gets the synch. particle 
#-----------------------------------------------------

print "Start.1"

b = Bunch()

print "Start.2"

syncPart = b.getSyncParticle()
print "Start.3"

m = syncPart.mass()
print "m=",m
syncPart.x(0.1)
syncPart.y(0.2)
syncPart.z(0.3)
print "x=",syncPart.x()
print "y=",syncPart.y()
print "z=",syncPart.z()
syncPart.px(1.1)
syncPart.py(1.2)
syncPart.pz(1.3)
print "px=",syncPart.px()
print "py=",syncPart.py()
print "pz=",syncPart.pz()

print "pVect=",syncPart.pVector()
print "rVect=",syncPart.rVector()
print "nxVect=",syncPart.nxVector()
(px,py,pz) = syncPart.pVector()
(nx,ny,nz) = syncPart.nxVector()
print "p*n=",(px*nx+py*ny+pz*nz)
print "abs(n)=",(nx*nx+ny*ny+nz*nz)
print "=============================="
syncPart.pVector((2,3,4))
syncPart.rVector((5,6,7))
syncPart.nxVector((0,1,0))
print "pVect=",syncPart.pVector()
print "rVect=",syncPart.rVector()
print "nxVect=",syncPart.nxVector()
print "=============================="
(px,py,pz) = syncPart.pVector()
(nx,ny,nz) = syncPart.nxVector()
print "p*n=",(px*nx+py*ny+pz*nz)
print "abs(n)=",(nx*nx+ny*ny+nz*nz)
print "=============================="

p = syncPart.momentum()
print "p=",p

orbit_mpi.finalize("The test is done!")

#memory test
#the amount of memory should be the same
nIter = 100000000
for i in xrange(nIter):
	syncPart = b.getSyncParticle()
	m = syncPart.mass()
	p = syncPart.momentum()
	syncPart.px(syncPart.px())
	syncPart.py(syncPart.py())
	syncPart.pz(syncPart.pz())
	syncPart.x(syncPart.x())
	syncPart.y(syncPart.y())
	syncPart.z(syncPart.z())
	beta = syncPart.beta()
	g = syncPart.gamma()
	e = syncPart.kinEnergy()
	p = syncPart.energyToMomentum(e)
	e = syncPart.momentumToEnergy(p)
	if i % 1000 == 0:
		print "i=",i

print "Stop."
