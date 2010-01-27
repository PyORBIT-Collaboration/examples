
#import laserstripping.DM_noLaserField
#reload(laserstripping)
#print dir()




from ext.las_str.SimplexMod import Simplex
import sys, math, os, orbit_mpi
from trackerrk4 import *
from laserstripping import *
from bunch import *
import time,os
from ext.las_str.part_generator import *
from ext.las_str.emittance import *
from orbit_mpi import mpi_comm,mpi_datatype,mpi_op



if not("n_states" in dir ()):


    n_states = 5
    levels = n_states*(1+n_states)*(1+2*n_states)/6
    
    orbit_path = os.environ["ORBIT_ROOT"]
    addr = orbit_path+"/ext/laserstripping/working_dir/"
    trans = orbit_path+"/ext/laserstripping/transitions/"
    
    St = Stark(trans, n_states)
#    design case
#    magnetic_field2 = RegularGridFS(addr+"D2_2126A_Data.table","m",0.01,0.0001*(46.2/53.1))
#    magnetic_field3 = RegularGridFS(addr+"D3_1449A_Data.table","m",0.01,0.0001*(42.0/28.3))
#    magnetic_field4 = RegularGridFS(addr+"D4_1737A_Data.table","m",0.01,0.0001*(46.2/39.4))

#    real case    
    magnetic_field2 = RegularGridFS(addr+"D2_2126A_Data.table","m",0.01,0.0001)
    magnetic_field3 = RegularGridFS(addr+"D3_1449A_Data.table","m",0.01,0.0001)
    magnetic_field4 = RegularGridFS(addr+"D4_1737A_Data.table","m",0.01,0.0001)

#magnetic_field2 = ConstEMfield(0.,0.,0.,0.,0.25,0.)
    
    magnetic_field2.setFieldOrientation(0,0.023,0,0,0,1,1,0,0)
    magnetic_field3.setFieldOrientation(-0.00772,0.023,1.814,0,0,1,1,0,0)
    magnetic_field4.setFieldOrientation(-0.09155,0.023,3.843,0,0,1,1,0,0)
    
field_cont = FieldSourceContainer()
field_cont.AddFieldSource(magnetic_field2)
field_cont.AddFieldSource(magnetic_field3)
field_cont.AddFieldSource(magnetic_field4)


#execfile("/home/tg4/workspace/PyOrbit/examples/ext/LaserStripping/SNS/tests.py")


print "start"
range_delta = []
for i in range(1500):
    range_delta.append(i/1000.)

#range_delta = [0]

for delta_z in range_delta:

#    z_foil1 = 0.476 + delta_z    #design case
    z_foil1 = 0.30699 + delta_z    #real case
#    z_foil2 = 2.936809
    z_foil2 = 4.430187+0.1  #flange
    
    
#    z_foil1 = -2.
#    z_foil2 = 2. - delta_z
    z_flange = 4.430187
    x_flange = -0.09679
    x_inj = 0.040
    y_inj = 0.023
    
    n_step = 10000
    TK = 1.0
    
    b = Bunch()
    b.charge(1)
    b.mass(0.938256)
    
    E = b.mass() + TK
    pz = math.sqrt(E*E - b.mass()*b.mass())
    beta = pz/E
    vz = beta*299792458
    time_step = (z_foil2 - z_foil1)/vz/n_step

    b.addParticle(x_inj,0.,y_inj,0.,z_foil1,pz)
    

    
    tracker = RungeKuttaTracker(100)
    walls = Walls()
    

    tracker.track(b, 0, time_step*n_step, time_step, field_cont,walls)
    
    print z_foil1, "\t",b.x(0), "\t",b.z(0),b.flag(0)


    
    
    
    
f = open(addr+'field.txt','w')
for i in range(100000):
    print >>f,z_foil1+(z_foil2 - z_foil1)*i/100000,"\t",field_cont.getFields(x_inj,y_inj,z_foil1+(z_foil2 - z_foil1)*i/100000,0)[4]
f.close()

