
#import laserstripping.DM_noLaserField
#reload(laserstripping)
#print dir()




from ext.las_str.SimplexMod import Simplex
import sys, math, os, orbit_mpi
from trackerrk4 import *
from laserstripping import *
from bunch import *
import orbit_mpi
import time,os
from ext.las_str.part_generator import *
from ext.las_str.ls_math import *
from ext.las_str.emittance import *
from ext.las_str.print_mod import printf
from ext.las_str.part_generator import BunchGen


    

bgen = BunchGen()

bgen.N_part = 1000

bgen.TK = 1.0                                    # [GeV]

bgen.mass = 0.938256 + 0.000511                 # [GeV]
bgen.charge = -1                                  # [e]


bgen.alphaX = 0.0                                 # [rad]
bgen.betaX = 10.0                                # [m]
bgen.emtX = 0.2e-6                             # [m*rad] non-norm

bgen.alphaY = -0.2                              # [rad]
bgen.betaY = 0.5                                 # [m]
bgen.emtY = 0.2e-6                             # [m*rad] non-norm


E = bgen.mass + bgen.TK
P = math.sqrt(E*E - bgen.mass*bgen.mass)
vz = 299792458*P/E

sigmaZ_beam = 30e-12*vz/(2*math.sqrt(2*math.log(2)))    #[m]
sigmaP_div_P = 0.5e-3   

bgen.alphaZ = 0.0
bgen.betaZ = sigmaZ_beam/sigmaP_div_P
bgen.emtZ = sigmaZ_beam*sigmaP_div_P


bgen.dispD = 0.0                                  # [m]
bgen.dispDP = -2.6                                # [rad]


b = bgen.getBunch(0)

for i in range(b.getSize()):    #consider all particle as stripped
    b.z(i,b.z(i) - 0.3) 



n_step = 10000


n_states = 3
levels = n_states*(1+n_states)*(1+2*n_states)/6

orbit_path = os.environ["ORBIT_ROOT"]

mag = SNSstrippingMagnet(orbit_path + "/ext/laserstripping/SNSstrippingMagnet/field_data")
St = Stark(orbit_path + "/ext/laserstripping/Hydrogen_data/StarkEG/", n_states)



eff = HminusStripping(1)    #par = 0 calculate probability and par = 1 randomly delete particle and track as neutral just after ionization



tracker = RungeKuttaTracker(1000)


time_step = 0.3/vz/n_step

tracker.track(b,0,time_step*n_step, time_step,mag, eff)









b.removeAllPartAttr()
b.addPartAttr("Populations",{"size":levels+1})  


polar = 1
excitation_eff = 0.9


####--------------------------distribution of populations for parallel and perpendicular polarization of laser field--------------------#### 
n = n_states
for i in range(b.getSize()):
    if (polar==0):#for the second level excitation of two levels with m=0
        for n1 in range(n):
            n2 = n-n1-1
            b.partAttrValue("Populations",i,1+n1+(n**3-n)/3, 3.*(n1-n2)*(n1-n2)/(n*(n*n-1)))

    if (polar==1):
        b.partAttrValue("Populations",i,1, 1.0 - excitation_eff)
        for n1 in range(n-1):
            n2 = n-n1-2
            b.partAttrValue("Populations",i,1+n+n1+(n**3-n)/3, excitation_eff*3.*(n1+1)*(n2+1)/(n*(n*n-1)))
            b.partAttrValue("Populations",i,2-n+n1+(n**3-n)/3, excitation_eff*3.*(n1+1)*(n2+1)/(n*(n*n-1)))





eff = DM_noLaserField(St, 1)


tracker.track(b,0,time_step*n_step, time_step,mag, eff)



#--------------bunch comments------------------
#At this point the bunch b is an output bunch after the stripping magnet

#the z axes location if the beam is about 0.3 m:
count = 0
z_aver = 0
for i in range(b.getSize()):
    z_aver += b.z(i)
z_aver = z_aver/b.getSize()
print "z_aver = ", z_aver

# The bunch particles have parameter of b.flag(i) = 1 that correspornd to the protons and b.flag(i) = 0 corresponds to the unstripped H0 particles
# That number depends mostly on excitation_eff. The strippinf efficiensy count_p / (count_p + count_H)) is about 0.98 for excitation_eff = 1.0 because of spontaneous (H0* -> H0 + photon) radiation losses 
count_H0 = 0
count_p = 0
for i in range(b.getSize()):
    if(b.flag(i) == 1):
        count_p += 1
    else:
        count_H0 += 1

print "number of protons = ", count_p
print "number of H0", count_H0

# the proton beam orbit has the following y-axes angle
count = 0
ang_y_aver = 0
for i in range(b.getSize()):
    if(b.flag(i) == 1):
        ang_y_aver += b.py(i)/b.pz(i)
ang_y_aver /= count_p
print "proton beam angle of orbit = ", ang_y_aver 



