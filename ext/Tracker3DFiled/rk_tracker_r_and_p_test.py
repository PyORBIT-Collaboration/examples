#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------
import sys
import math
from bunch import Bunch

from tracker3dfield import RungeKuttaTracker
from tracker3dfield import PyExternalEffects
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
		pass
		
	def finalizeEffects(self,bunch):
		pass		
		
	def applyEffects(self,bunch, index, inVct, outVct, t, t_step, field_source, tracker):
		print " %5d "%self.count," %8.6f  %8.6f  %8.6f "%(bunch.x(0),bunch.y(0),bunch.z(0))
		self.count = self.count + 1
		pass
		
print "Start."


b = Bunch()
print "Part. m=",b.mass()
print "Part. q=",b.charge()

TK = 1.0
E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())
c_light = 2.99792458e+8

print "TK[GeV] = ",TK
print "P[GeV/c] = ",P

b.addParticle(0.,0.,0.,0.,0.,P)
b.compress()

# radius estimation
R = P*1.0e9/(c_light*b.charge()*1.0)
print "R[m] = ",R
n_step = 1000
time_step = (2*3.1415926*R/(c_light*P/E))/n_step

fS = FieldSource()
extEff = ExternalEffects()
print "ExternalEffects name=",extEff.name()

tracker = RungeKuttaTracker(10000.0)
print "Entrance plane (a,b,c,d)=",tracker.entrancePlane()
print "Exit     plane (a,b,c,d)=",tracker.exitPlane()
print "Length[m]=",tracker.length()

print "Start tracking."
print "================================================"
print "Step_Index    x     y   z "
tracker.track(b,time_step*n_step, time_step,fS,extEff)
print "================================================"
print "Stop tracking."

print "time step=",tracker.timeStep()
print "Stop."

sys.exit(1)
