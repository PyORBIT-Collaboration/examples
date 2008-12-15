#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------
import sys
import math

from bunch import *
from tracker3dfield import *
from orbit_utils import *
import mygra
import os
import orbit_mpi


#print dir()

print "Start."


b = Bunch()
b.addPartAttr("Amplitudes")

b.charge(0)
print "Part. m=",b.mass()
print "Part. q=",b.charge()

TK = 1.0
E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())
c_light = 2.99792458e+8

print "TK[GeV] = ",TK
print "P[GeV/c] = ",P

b.addParticle(0.,0.,0.,0.,0.,P*0.)
b.compress()

b.partAttrValue("Amplitudes",0,1,1.0)
print "AttrValue=", b.partAttrValue("Amplitudes",0,1)
print "AttrSize=",b.getPartAttrSize("Amplitudes")

# radius estimation
R = P*1.0e9/(c_light*(b.charge()+1.0)*1.0)
print "R[m] = ",R
n_step = 10000
time_step = (2*3.1415926*R/(c_light*P/E))/n_step/1;
time_step=(2*3.1415926/1e+12)/n_step;

fS=CppBaseFieldSource()


First = LasStripExternalEffects(0.0005,1,102.5e-9,"/home/tg4/transitions/",3)
 


#First = LasStripExternalEffects(0.0005,1,102.5e-9)
First.name("first_effect")


print "ExternalEffects name=",First.name()



tracker = RungeKuttaTracker(1000.0)
print "Entrance plane (a,b,c,d)=",tracker.entrancePlane()
print "Exit     plane (a,b,c,d)=",tracker.exitPlane()
print "Length[m]=",tracker.length()

print "Start tracking."
print "==========================================================================================="
print "Step_Index    x                                y                              z "
tracker.track(b,0,10000*time_step*n_step, time_step,fS,First)
print "==========================================================================================="
print "Stop tracking.",time_step*n_step

print "time step=",tracker.timeStep()
print "Stop."
print "time_fl=",time_step*n_step
#print b.partAttrValue("Amplitudes",0,2)*b.partAttrValue("Amplitudes",0,2)+ b.partAttrValue("Amplitudes",0,3)*b.partAttrValue("Amplitudes",0,3)
print "AttrValue=", 1-b.partAttrValue("Amplitudes",0,1)


mygra.PlotPopl()

os.system('gthumb image.png')

os.remove('data_ampl.txt')


#orbit_mpi.finalize()
sys.exit(1)
