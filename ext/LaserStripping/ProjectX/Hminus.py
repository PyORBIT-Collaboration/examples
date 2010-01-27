
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




orbit_path = os.environ["ORBIT_ROOT"]
addr = orbit_path+"/ext/laserstripping/working_dir/"
trans = orbit_path+"/ext/laserstripping/transitions/"
    



pf = printf("Hminus.dat","By [T]", "rmsAng. [mrad]", "emitt. [mm mrad]","em_norm. [mm mrad]", "efficiency")

#for n_states, polar, By  in [(n_states, polar, By)
#                        for n_states in [2, 3, 4, 5]
#                        for polar in [0, 1]
#                        for By in [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 
#                                   1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1]]:

for Bx in [0.45]:
    
    n_step = 10000000
    Nevol=10000
    
    h = 0.050
    
    z1 = -2*h
    z2 = 10*h
    
    
    
    
    b = Bunch()
    b.charge(0)
    b.mass(0.938256 + 0.000511)  
    
    TK = 8.0
    E = b.mass() + TK
    P = math.sqrt(E*E - b.mass()*b.mass())
    b.addParticle(0,0,0,0,z1,P)
    
    c_light = 2.99792458e+8
    beta = P/E
    gamma = 1/math.sqrt(1-beta*beta)
    vz = beta*c_light
    
 
    

    
    
    mag = FringeField(h,Bx,1)



    evo = RecordEvolution("Populations",0,Nevol)

    eff = HminusStripping(0)

    

    cont_eff = ExtEffectsContainer()
    cont_eff.AddEffect(evo)
    cont_eff.AddEffect(eff)

    
    
    tracker = RungeKuttaTracker(1000)
    
    
    time_step = (z2 - z1)/vz/n_step
    tracker.track(b,0,time_step*n_step, time_step,mag, cont_eff)
        
#    print b.px(0)/b.pz(0)


    ang = [0]*(Nevol+1)
    fd = [0]*(Nevol+1)
    dz = (z2 - z1)/Nevol
    prob = b.partAttrValue("Evolution",0,Nevol)
    


    ang_aver = 0
    for i in range(Nevol-1,-1,-1):
        fd[i] = abs((b.partAttrValue("Evolution",0,i+1) - b.partAttrValue("Evolution",0,i))/dz/prob)
        ang[i] = ang[i+1] + mag.getField(0,0,z1 + i*dz,0)*dz*c_light/(P*1.0e9)
        ang_aver += ang[i]*fd[i]*dz       
      
    ang_aver2 = 0 
    for i in range(Nevol):
        ang_aver2 += (ang[i] - ang_aver)*(ang[i] - ang_aver)*fd[i]*dz
   
    print prob
    
    if(prob>0.01):
#        print Bx,1e3*math.sqrt(ang_aver2),ang_aver2*10*1e6, beta*gamma*ang_aver2*10*1e6,prob
        pf.fdata(Bx,1e3*math.sqrt(ang_aver2),ang_aver2*10*1e6, beta*gamma*ang_aver2*10*1e6,prob)  

    
   
    
    f = open('evol_out.txt','w')
    for i in range(Nevol):
        print >>f, abs(b.partAttrValue("Evolution",0,i))
    f.close()
    
    
"""    
    gen = BunchGen()
    bunch, bunch_unstr = gen.getAutoionizationBunch(100000,b,0)
    
        
    tracker.track(bunch,0,time_step*n_step, time_step*100,mag)
    
    if (bunch.getSize() > 0):
        ang_aver = 0
        for i in range(bunch.getSize()):
            ang_aver += bunch.px(i)/bunch.pz(i)
        ang_aver = ang_aver/bunch.getSize()
        
        ang_aver2 = 0
        for i in range(bunch.getSize()):
            ang_aver2 += (bunch.px(i)/bunch.pz(i) - ang_aver)*(bunch.px(i)/bunch.pz(i) - ang_aver)
        ang_aver2 = ang_aver2/bunch.getSize()
        
        pf.fdata(bunch.getSize(),n_states,polar,By,1e3*math.sqrt(ang_aver2),beta*gamma*ang_aver2*20*1e6,bunch.getSize()/100000.)
        
""" 
 

print "end"

