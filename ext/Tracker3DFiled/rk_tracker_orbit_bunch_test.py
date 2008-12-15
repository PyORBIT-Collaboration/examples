#-----------------------------------------------------
#Track traditional ORBIT bunch through the external field
#-----------------------------------------------------
import sys
from bunch import Bunch

from trackerrk4 import RungeKuttaTracker
from trackerrk4 import PyExternalEffects
from orbit_utils import PyBaseFieldSource

#the implementation of the field source
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

class ExternalEffects(PyExternalEffects):
	def __init__(self):
		PyExternalEffects.__init__(self)
		self.count = 0
		
	def setupEffects(self,bunch):
		print "Set up external effects, bunch momentum[GeV/c]=",bunch.getSyncParticle().momentum()
		
	def finalizeEffects(self,bunch):
		print "Finilize external effects, bunch momentum[GeV/c]=",bunch.getSyncParticle().momentum()
		print "Finilize external effects, bunch size=",bunch.getSize()
		print "Finilize external effects, (x,xp,y,yp,z,dE)=",(bunch.x(0),bunch.xp(0),bunch.y(0),bunch.yp(0),bunch.z(0),bunch.dE(0))
		
		
	def applyEffects(self,bunch, index, inVct, outVct, t, t_step, field_source, tracker):
		self.count = self.count + 1
		
print "Start."

b = Bunch()
b.addParticle(0.,0.,0.,0.,0.,0.)
b.compress()


syncPart = b.getSyncParticle()
syncPart.kinEnergy(1.0)
print "synch. part. m=",syncPart.mass()
print "synch. part. p=",syncPart.momentum()


fS = FieldSource()
extEff = ExternalEffects()
print "ExternalEffects name=",extEff.name()

tracker = RungeKuttaTracker(2.0)
print "Entrance plane (a,b,c,d)=",tracker.entrancePlane()
print "Exit     plane (a,b,c,d)=",tracker.exitPlane()
print "Spatial Eps[m]=",tracker.spatialEps(0.0001)
print "Length[m]=",tracker.length()
print "steps number=",tracker.stepsNumber()
print "tracker.isOutside(0,0,2)=",tracker.isOutside(0,0,2)
print "tracker.isAfterEntrance(0,0,2)=",tracker.isAfterEntrance(0,0,2)
print "tracker.isBeforeExit(0,0,2)=",tracker.isBeforeExit(0,0,2)

print "Start tracking."
tracker.trackBunch(b,fS,extEff)
print "Stop tracking."

print "time step=",tracker.timeStep()
print "steps number=",tracker.stepsNumber()
print "Stop."

sys.exit(1)
