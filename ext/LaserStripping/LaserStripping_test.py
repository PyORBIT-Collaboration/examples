#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------
addr="/home/tg4/workspace/PyOrbit/ext/laserstripping/"
import sys
sys.path.append(addr+"/PyModules/")
import math

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *

import mygra
import os
import orbit_mpi


#print sys.path
#print dir(math)



#Parameters of laser beam and particle
z0=1.0e-2
TK = 1.0*(1+0.0000)
fx=fy=-1
alpha=39.39
alpha=39.5
wx=59.0e-6
wy=850.0e-6
lambd=355.0e-9
power=1000000.
n_step = 1000000





print "Start."

b = Bunch()
b.addPartAttr("Amplitudes")

b.charge(0)
print "Part. m=",b.mass()
print "Part. q=",b.charge()


E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())


print "TK[GeV] = ",TK
print "P[GeV/c] = ",P



b.addParticle(0.,0.,0.,0.,-z0,P)
b.partAttrValue("Amplitudes",0,1,1.0)
print "AttrValue=", b.partAttrValue("Amplitudes",0,1)
print "AttrSize=",b.getPartAttrSize("Amplitudes")




fS=LSFieldSource()

LFS=HermiteGaussianLFmode(math.sqrt(power),0,0,wx,wy,fx,fy,lambd) # (sqrt(P),n,m,wx,wy,fx,fy,lambda)


kz=-1/math.tan(math.radians(alpha))


LFS.setLaserFieldOrientation(0.,0.,0.,#{x0,y0,z0}
                             -1.,0.,kz,#{kx,ky,kz}     
                             1.,0.,1/kz,#{mx,my,mz}    Please be shure that kx*mx+ky*my+kz*mz==0
                             0.,1.,0.)#{Ex,Ey,Ez}    Please be shure that kx*Ex+ky*Ey+kz*Ez==0




First = LasStripExternalEffects(LFS,addr+"/transitions/",3,1000.) 


#First = LasStripExternalEffects(0.0005,1,102.5e-9)
#First.name("first_effect")
#print "ExternalEffects name=",First.name()



tracker = RungeKuttaTracker(1000.0)


print "Start tracking."
print "==========================================================================================="
print "Step_Index    x                                y                              z "


#n_step = 100000
#time_step = (2*3.1415926*R/(c_light*P/E))/n_step/1;
#time_step=(2*3.1415926/1e+12)/n_step;
#par=10;
#tracker.track(b,-par*time_step*n_step,2*par*time_step*n_step, time_step,fS,First)


vz=299792458*P/math.sqrt(b.mass()*b.mass()+P*P)


time_step = (2*z0/vz)/n_step;

tracker.track(b,0,time_step*n_step, time_step,fS,First)



print "==========================================================================================="
print "Stop tracking.",time_step*n_step


#print b.partAttrValue("Amplitudes",0,2)*b.partAttrValue("Amplitudes",0,2)+ b.partAttrValue("Amplitudes",0,3)*b.partAttrValue("Amplitudes",0,3)
print "AttrValue=", 1-b.partAttrValue("Amplitudes",0,1)


graph = mygra.PlotPopl()

os.system('gthumb /home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/image.png')
os.remove(addr+"/working_dir/data_ampl.txt")


#orbit_mpi.finalize()
sys.exit(1)
