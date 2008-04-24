#-----------------------------------------------------
#Track bunch through the external field
#-----------------------------------------------------
import sys
from bunch import Bunch

from tracker3dfield import RungeKuttaTracker
from orbit_utils import PyBaseFieldSource

class FieldSource(PyBaseFieldSource):
	def __init__(self):
		PyBaseFieldSource.__init__(self)
	
	def getElectricField(self,x,y,z,t):
		#print "Efield x,y,z, t=",(x,y,z,t)
		fx = 0.
		fy = 0.
		fz = 0.
		return (fx,fy,fz)
		
	def getMagneticField(self,x,y,z,t):
		#print "Bfield x,y,z, t=",(x,y,z,t)
		fx = 0.
		fy = 1.
		fz = 0.
		return (fx,fy,fz)		

print "Start."

b = Bunch()
b.addParticle(0.,0.,0.,0.,0.,0.)
b.compress()


syncPart = b.getSyncParticle()
syncPart.kinEnergy(1.0)
print "synch. part. m=",syncPart.mass()
print "synch. part. p=",syncPart.momentum()


fS = FieldSource()

tracker = RungeKuttaTracker(2.0)

print "Start tracking."
tracker.track(b,fS)
print "Stop tracking."

print "Stop."

sys.exit(1)
