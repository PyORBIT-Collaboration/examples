
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





print "start"

#execfile("/home/tg4/workspace/PyOrbit/examples/ext/LaserStripping/SNS/RingInjection.py")
    

#z_foil1 = 0.476    #design case
z_foil1 = 0.30699    #real case
z_foil2 = 2.936809



z_flange = 4.430187
x_flange = -0.09679
x_inj = 0.040
y_inj = 0.023
    
    
    
n_step = 1000000
Nevol = 100000
coeff_bunch = 100000
    
    #print field_cont.getFields(0.04321,0,z_foil1,0)
    
    

   
    ######-------------------Bunch generation----------------------------------
    
a = BunchGen()
    
a.N_part = 1
    
a.TK = 1.0                                    # [GeV]
    
a.mass = 0.938256 + 0.000511                  # [GeV]
a.charge = 0                                  # [e]
    
a.alphaX = 0.0586419                                 # [rad]
a.betaX = 10.9436                                # [m]
a.emtX = 0.3e-6                                 # [m*rad]
a.cutOffX = math.sqrt(a.emtX*a.betaX)*3.0     # [m]
    
a.alphaY = 0.0543763                                  # [rad]
a.betaY = 12.7927                                 # [m]
a.emtY = 0.3e-6                                 # [m*rad]
a.cutOffY = math.sqrt(a.emtY*a.betaY)*3.0     # [m]
    
a.relativeSpread = 0.5e-3*0
a.sigma_beam = 0
a.cutOffZ = 3.
    
a.dispD = 0.                                  # [m]
a.dispDP = 0.                                 # [rad]
    

b = a.getBunch(0)

    
for i in range(b.getSize()):
    b.x(i,b.x(i) + x_inj)
    b.y(i,b.y(i) + y_inj)
    b.z(i,b.z(i) + z_foil1)
        
    for i in range(b.getSize()):
        b.px(i,0)
        b.py(i,0) 

    
######-------------------Bunch generation----------------------------------
    


  
    
####--------------------------uniform distribution--------------------####    

b.addPartAttr("Populations",{"size":levels+1})      
const = 0
for i in range(1,n_states+1):
        const += math.pow(i,-2.8)
    
for n in range(1,n_states+1):
    for m in range(-(n-1),(n-1)+1):
        for n1 in range(n-abs(m)-1+1):
            for np in range(b.getSize()):
                b.partAttrValue("Populations",np,1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2,math.pow(n,-2.8)/(const*n*n))
                   
####--------------------------uniform distribution--------------------####    

    



E = b.mass() + a.TK
beta = math.sqrt(E*E - b.mass()*b.mass())/E
vz = beta*299792458
time_step = (z_foil2  - z_foil1)/vz/n_step
    

evo = RecordEvolution("Populations",0,Nevol)
eff = DM_noLaserField(St) 
    
cont_eff = ExtEffectsContainer()
cont_eff.AddEffect(evo)
cont_eff.AddEffect(eff)
    
    
tracker = RungeKuttaTracker(100)
    

tracker.track(b, 0, time_step*n_step, time_step, field_cont, cont_eff)

    
f = open('evol_out.txt','w')
for i in range(Nevol+1):
    print >>f, b.partAttrValue("Evolution",0,i)
f.close()
    
sys.exit(0)
print "begin generating"



n_step = 10000

walls = Walls() 

bunch, bunch_unstr = a.getAutoionizationBunch(coeff_bunch,b,0)
###    bunch = b
###    bunch.charge(+1)
  
print "begin tracking"
time_step = (z_flange + 0.01 - z_foil1)/vz/n_step
tracker.track(bunch, 0, time_step*n_step, time_step, field_cont, walls)

#    f = open('x_xp.txt','w')
#    for i in range(bunch.getSize()):
#        print >>f,bunch.x(i)-x_flange,"\t",bunch.px(i)/bunch.pz(i)
#    f.close()

#    f = open('y_yp.txt','w')
#    for i in range(bunch.getSize()):
#        print >>f,bunch.y(i),"\t",bunch.py(i)/bunch.pz(i)
#    f.close()   


    
rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
    
n1 = b.getSizeGlobal()*coeff_bunch
n2 =  bunch.getSizeGlobal()
    
    
bunch.compress()
n3 = bunch.getSizeGlobal()
    
if (rank==0):
    print z_foil1,"    n_H0 = ",n1,"    n_p = ",n2,"    n_p(loss) = ",n3
    


#print z_foil1,"    ",1.*(n0 - nf)/n0

"""
f = open('flags.txt','w')
for i in range(b.getSize()):
    print >>f,b.flag(i)
f.close()
"""
"""    
    f = open('x_xp.txt','w')
    for i in range(b.getSize()):
        print >>f,bunch.x(i)-x_flange,"\t",bunch.px(i)/bunch.pz(i)
    f.close()

    f = open('y_yp.txt','w')
    for i in range(bunch.getSize()):
        print >>f,bunch.y(i),"\t",bunch.py(i)/bunch.pz(i)
    f.close()
"""


"""
for i in range(b.getSize()):
    print b.x(i)

"""
#print b.partAttrValue("Evolution",0,100)



#b.dumpBunch("bunch_loss_not.dat")


#print b.x(0)
#print b.y(0)
#print b.z(0)
#print b.getSize()
#b.compress()
#print b.getSize()


"""
for i in range(b.getSize()):
    print i," " ,z_foil1+(z_foil2 - z_foil1)*i/1000," ",b.flag(i)
"""

"""
f = open('field.txt','w')
for i in range(100000):
    print >>f,field_cont.getFields(-0.04321,0,z_foil1+(z_foil2 - z_foil1)*i/100000,0)[4]
f.close()
"""            
 
"""
f = open('evol_out.txt','w')
for i in range(Nevol+1):
    print >>f, b.partAttrValue("Evolution",0,i)
f.close()
"""
    

#execfile("/home/tg4/workspace/PyOrbit/examples/ext/LaserStripping/RingInjection.py")
#execfile("/home/tg4/workspace/laser_stripping/src/parameters.py")


