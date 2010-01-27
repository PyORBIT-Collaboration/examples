#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------

import sys, math, os, orbit_mpi

from laserstripping import *
from bunch import *
from orbit_mpi import *
from orbit_utils import *
from ext.las_str.TwoLevelFuncMod import TwoLevelFunc
from ext.las_str.SchredingerFuncMod import SchredingerFunc
from ext.las_str.SimplexMod import Simplex
from ext.las_str.part_generator import BunchGen
from ext.las_str.print_mod import printf



time_start = orbit_mpi.MPI_Wtime()
rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)
orbit_path = os.environ["ORBIT_ROOT"]






#----------------------Beginning of the beam parameters----------------------#

b = BunchGen()

b.N_part = 25

b.TK= 1.0                                    # [GeV]

b.mass = 0.938256 + 0.000511                  # [GeV]
b.charge = 0                                  # [e]

b.alphaX = 0.                                 # [rad]
b.betaX = 25.0                                # [m]
b.emtX = 0.2e-6                             # [m*rad]
b.cutOffX = math.sqrt(b.emtX*b.betaX)*3.0     # [m]


b.alphaY = 2.0                              # [rad]
b.betaY =0.75                                 # [m]
b.emtY = 0.2e-6                             # [m*rad]
b.cutOffY = math.sqrt(b.emtY*b.betaY)*3.0     # [m]



b.relativeSpread = 0.5e-3

b.dispD = 0.                                  # [m]
b.dispDP = 2.6                                # [rad]

b.sigma_beam = 30e-12/(2*math.sqrt(2*math.log(2)))                          #[m]
b.cutOffZ = 3*b.sigma_beam



#----------------------End of the beam parameters----------------------#

#----------------------Beginning of the Laser field and excitation parameters of H0 atom----------------------#


H13 = SchredingerFunc(orbit_path+"/ext/laserstripping/transitions/",3,300)

H13.method = 2







H13.Bx = 0.005e-3    
H13.fS = ConstEMfield(0.,0.,0.,H13.Bx,0.,0.)                                   # [T]
    
#grad_B = 1.0                                    # [T/m]
#H13.fS = QuadEMfield()
#H13.fS.cxBz(grad_B)
#H13.fS.czBx(grad_B)



H13.bunch = b.getBunch(0)


#f = open('bunch_parameters.txt','w')
#for i in range(H13.bunch.getSize()):
#    print >>f, i,"\t",(math.sqrt(H13.bunch.pz(i)*H13.bunch.pz(i) + H13.bunch.mass()*H13.bunch.mass()) - H13.bunch.mass() - 1),"\t",H13.bunch.x(i),"\t",H13.bunch.px(i)/H13.bunch.pz(i),"\t",H13.bunch.y(i),"\t",H13.bunch.py(i)/H13.bunch.pz(i)  
#f.close()


H13.TK = b.TK
H13.sigma_beam = b.sigma_beam

H13.n_sigma = 3
H13.n_step = 1000000

H13.la = 355.0e-9                               # [m]
H13.power = 2.2e6                               # [W]

H13.fx = -4.5                                      # [m]
H13.fy = -4.5                                      # [m]

H13.wx = 370.0e-6                               # [m]
H13.wy = 700.0e-6                               # [m]

    
H13.env_sigma = 55e-12/(2*math.sqrt(2*math.log(2)))  #[m]



H13.rx = 1.69e-3                               # [m]
H13.ry = 0.77e-3                               # [m]

H13.ax = 0.36e-3                               # [rad]
H13.ay = 0.24e-3                               # [rad]

#----------------------End of the Laser field and excitation parameters of H0 atom----------------------#

#H13.bunch.dumpBunch("bunch_init1.dat")




#-------------------definition of the optimization function----------------------------#
#pf = printf("optimiz.dat","N_part","cpu_time", "W[MW]", "wx[um]", "wy[um]", "fx[cm]", "fy[cm]", "Population", "+- Err")
pf = printf("optimiz_sh1.0_test.dat","N_step","cpu_time", "W[MW]","rx[mm]", "ry[mm]", "ax[mrad]", "ay[mrad]", "Population", "+- pop", "P_ioniz", "+- P_ioniz")

#name_args, guess, increments = ['wx','wy'],[300.0e-6, 300.0e-6],[10e-6, 100e-6]
#name_args, guess, increments  = ['wx','wy','fx'], [333.0e-6, 333.0e-6,-2.800], [10e-6, 100e-6,1.00]
name_args, guess, increments  = ['rx','ry','ax','ay'], [0.7e-3,0.7e-3,0.3e-3,0.3e-3], [1.0e-4,1.0e-4,1.0e-5,1.0e-5]
#name_args, guess, increments  = ['rx','ry','ax','ay'], [68.45e-3,0.377e-3,0.2e-3,0.2e-3], [1.0e-4,1.0e-4,1.0e-4,1.0e-4]
#name_args, guess, increments  = ['wx','fx'], [370.0e-6, -4.500], [10e-6,1.00]
#name_args, guess, increments  = ['wx','wy','fx','fy'], [300.6e-6, 300.6e-6,-2.0,-2.0], [10e-6, 10e-6,1.00,1.00]



#powers = [0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.6e6, 0.7e6, 0.8e6, 0.9e6, 1.0e6]
#powers = [0.5e6]
#powers = [0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.6e6, 0.7e6, 0.8e6, 0.9e6, 1.0e6, 1.1e6, 1.2e6, 1.3e6, 1.4e6, 1.5e6, 1.6e6, 1.7e6, 1.8e6, 1.9e6, 2.0e6]

def opt_func(args):

    for i in range(len(args)):
        H13.__dict__[name_args[i]] = args[i]

    H13.ay = 0
        
    
    sum = 0
    for N_p in range(len(powers)):
        H13.power = powers[N_p]
        pop, sigma_pop = H13.population()
        sum += pop
#        pf.fdata("","",H13.power/1.0e+6,"","","","",pop,sigma_pop)
    H13.count += 1
#    pf.fdata(H13.count,orbit_mpi.MPI_Wtime() - time_start,"aver_0.1-1.0",1.0e+6*H13.wx,1.0e+6*H13.wy,H13.fx*100,H13.fy*100,sum/len(powers))
    pf.fdata(H13.count,orbit_mpi.MPI_Wtime() - time_start,"aver_0.1-2.0",1.0e+3*H13.rx,1.0e+3*H13.ry,1.0e+3*H13.ax,1.0e+3*H13.ay,sum/len(powers))

    return -sum/len(powers)

#-------------------definition of the optimization function----------------------------#


def opt_func1(args):

    for i in range(len(args)):  
        H13.__dict__[name_args[i]] = args[i]

        
    
    pop, sigma_pop, p_ioniz, sigma_p_ioniz  = H13.population()
    pf.fdata(H13.count,orbit_mpi.MPI_Wtime() - time_start,H13.power/1.0e+6,1.0e+3*H13.rx,1.0e+3*H13.ry,1.0e+3*H13.ax,1.0e+3*H13.ay,pop,sigma_pop,p_ioniz, sigma_p_ioniz)
    H13.count += 1

    return -pop

#-------------------definition of the optimization function----------------------------#




#-----------------------Beginning of optimization-----------------------------------------

#s = Simplex(opt_func1, guess, increments)
#(values, err, iter) = s.minimize(1e-20, 1000,0)

#-----------------------End of Optimization-----------------------------------------------









#-----------------------Single point calculation-----------------------------------------------------
#print "start"
#print  H13.population() 

"""
for i in range(0,50):
    H13.B_x = i*0.001+1.0e-10       
    Bx = i*0.001+1.0e-10
    H13.fS = ConstEMfield(0.,0.,0.,Bx,0.,0.)
    pop, sigma_pop = H13.population()
    if(rank==0):    print i*0.001,"  ",pop,"  ",sigma_pop
"""
#-----------------------End of single point calculation-----------------------------------------------





#-----------------------Single point calculation with optimization function-----------------------------------------------------

pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
print pop

#-----------------------End of Single point calculation with optimization function-----------------------------------------------








#-----------------------Slide show-----------------------------------------
"""
res_dir = os.environ["ORBIT_ROOT"]+"/ext/laserstripping/working_dir/"
data_name = "data"
pic_name = "pic"


     
for i in range(b.N_part):
     popi = 1 - H13.bunch_target.partAttrValue("Populations",i,0) - H13.bunch_target.partAttrValue("Populations",i,1)
     P = math.sqrt((H13.bunch.mass() + b.TK)*(H13.bunch.mass() + b.TK) - H13.bunch.mass()*H13.bunch.mass())
     PlotPopl([3,1],["eff= %1.6f"%popi,"px= %1.6f"%H13.bunch.px(i),"py= %1.6f"%H13.bunch.py(i),"x= %1.6f"%H13.bunch.x(i),"y= %1.6f"%H13.bunch.y(i),"pz-p0= %1.6f"%(H13.bunch.pz(i) - P),"-px/(pz-p0)= %1.6f"%(-H13.bunch.px(i)/(H13.bunch.pz(i) - P)),"d= %1.6f"%(abs(H13.bunch.px(i)+2.6*(H13.bunch.pz(i) - P))/math.sqrt(1+2.6*2.6))],0.15,res_dir+data_name+"%i.dat"%(i*size+rank),res_dir+pic_name+"%i.png"%(i*size+rank))
     os.remove(res_dir+data_name+"%i.dat"%(i*size+rank))
"""
#-----------------------Slide show------------------------------------------












#pf = printf("magnetic_field.dat","dispDP", "relativeSpread", "population")
#pf = printf("magnetic_field.dat","power [MW]", "B [T]", "population", "+- Err")
"""
for spread, disp  in [(spread,disp)
                        for spread in [0.1e-3,0.2e-3,0.3e-3,0.4e-3,0.5e-3,0.6e-3,0.7e-3,0.8e-3,0.9e-3,1.0e-3,1.1e-3,1.2e-3,1.3e-3,1.4e-3,1.5e-3,1.6e-3,1.7e-3,1.8e-3,1.9e-3,2.0e-3]
                        for disp in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0]]:
"""

"""

for B, power  in [(B,power)
                        for B in [0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.008,0.01]
                        for power in [1.0e6,2.0e6,3.0e6,4.0e6,5.0e6,6.0e6,7.0e6,8.0e6,9.0e6,10.0e6]]:



    By = B                    
    H13.fS = ConstEMfield(0.,0.,0.,0.,By,0.)
    
    H13.power = power
#    H13.bunch = b.getBunch(0)
    pop, sigma_pop = H13.population()
    
    pf.fdata(H13.power/1.0e6,By,pop,sigma_pop)


"""
