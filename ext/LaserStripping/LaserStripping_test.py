#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------
import os
orbit_path = os.environ["ORBIT_ROOT"]
addr=orbit_path+"/ext/laserstripping/"
import sys
sys.path.append(addr+"/PyModules/")
import math

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *


import os
import orbit_mpi


#print sys.path
#print dir(math)
#print dir()


#Parameters of laser beam and particle

#la0=102.5175436158152e-9
z0=0.5e-3
Bx=1./10000.
n_states=3
TK = 1.0*(1+0.0001)
fx=fy=-0.00
la=355.0e-9
wx=94.4e-06 
wy=850.0e-6
power=1000000.
n_step = 100000





print "Start."







b = Bunch()



print "Part. m=",b.mass()
print "Part. q=",b.charge()


E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())


print "TK[GeV] = ",TK
print "P[GeV/c] = ",P

#alpha=360*math.acos((b.mass()*(la/la0-1)-TK)/P)/2/math.pi



#print "AttrValue=", b.partAttrValue("Amplitudes",0,1)
#print "AttrSize=",b.getPartAttrSize("Amplitudes")






Stark=HydrogenStarkParam(addr+"/transitions/",n_states)

#print "alpha=",alpha
vz=299792458*P/math.sqrt(b.mass()*b.mass()+P*P)
time_step = (2*z0/vz)/n_step;
    
  
    


for i in range(1):
    
    
#    Bx=0.001*i
#    wx=59.0e-7*i

#    Bx=0.001
    
    b.addParticle(0.,0.,0.,0.,-z0,P)
    b.charge(0)
    b.addPartAttr("Amplitudes")
    b.partAttrValue("Amplitudes",0,1,1.0)
    
    fS=LSFieldSource(0.,0.,0.,Bx,0.,0.)
    
    la0= 2*math.pi*5.291772108e-11/7.297352570e-3/(Stark.getStarkEnergy(b.mass(),
                           2,0,0,#{n1,n2,m}
                           0.,0.,0.,#{Ex,Ey,Ez}
                           Bx,0.,0.,#{Bx,By,Bz}
                           0.,0.,P)#{px,py,pz}
    -Stark.getStarkEnergy(b.mass(),
                           0,0,0,#{n1,n2,m}
                           0.,0.,0.,#{Ex,Ey,Ez}
                           Bx,0.,0.,#{Bx,By,Bz}
                           0.,0.,P)#{px,py,pz}
    )
#print la0
    alpha=360*math.acos((b.mass()*(la/la0-1)-TK)/P)/2/math.pi
    alpha=39.3850143245


    kz=-1/math.tan(math.radians(alpha))

    

    LFS=HermiteGaussianLFmode(math.sqrt(power),0,0,wx,wy,fx,fy,la) # (sqrt(P),n,m,wx,wy,fx,fy,lambda)
    LFS.setLaserFieldOrientation(0.,0.,0.,#{x0,y0,z0}
                             -1.,0.,kz,#{kx,ky,kz}     
                             1.,0.,1/kz,#{mx,my,mz}    Please be shure that kx*mx+ky*my+kz*mz==0
                             0.,1.,0.)#{Ex,Ey,Ez}    Please be shure that kx*Ex+ky*Ey+kz*Ez==0

    First = LasStripExternalEffects(LFS,Stark,1000.) 




#First = LasStripExternalEffects(0.0005,1,102.5e-9)
#First.name("first_effect")
#print "ExternalEffects name=",First.name()



    tracker = RungeKuttaTracker(1000.0)


#print "Start tracking."
#print "==========================================================================================="
#print "Step_Index    x                                y                              z "





    tracker.track(b,0,time_step*n_step, time_step,fS,First)



#print "==========================================================================================="
#print "Stop tracking.",time_step*n_step


#print b.partAttrValue("Amplitudes",0,2)*b.partAttrValue("Amplitudes",0,2)+ b.partAttrValue("Amplitudes",0,3)*b.partAttrValue("Amplitudes",0,3)
    population=1-b.partAttrValue("Amplitudes",0,1)
#    file=open("020_tuning.txt",'a')
#    file.write('%f'%Bx+"\t")
#    file.write('%f'%population+"\n")
#    file.close() 
    print   wx,"    ",population
    b.deleteParticle(0)

   
print "AttrValue=", population





#orbit_mpi.finalize()
sys.exit(1)
