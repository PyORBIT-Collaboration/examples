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
from ext.las_str.plot_mod import *


time_start = orbit_mpi.MPI_Wtime()
rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)
orbit_path = os.environ["ORBIT_ROOT"]






#----------------------Beginning of the beam parameters----------------------#

b = BunchGen()

b.N_part = 250

b.TK= 8.0                                     # [GeV]

b.mass = 0.938256 + 0.000511                  # [GeV]
b.charge = 0                                  # [e]

b.alphaX = 0.0                                 # [rad]
b.betaX = 40.0                                # [m]
b.emtX = 0.5e-6                             # [m*rad]
b.cutOffX = math.sqrt(b.emtX*b.betaX)*3.0     # [m]


b.alphaY = 0.0                              # [rad]
b.betaY = 10.0                                 # [m]
b.emtY = 0.5e-6                             # [m*rad]
b.cutOffY = math.sqrt(b.emtY*b.betaY)*3.0     # [m]



b.relativeSpread = 0.25e-3

b.dispD = 0.0                                  # [m]
b.dispDP = 0.0                                # [rad]

b.sigma_beam = 65e-12                           #[m]
b.cutOffZ = 3*b.sigma_beam



#----------------------End of the beam parameters----------------------#

#----------------------Beginning of the Laser field and excitation parameters of H0 atom----------------------#


H13 = SchredingerFunc(orbit_path+"/ext/laserstripping/transitions/",2,1000)

H13.method = 1




H13.Bx = 0.001
H13.fS = ConstEMfield(0.,0.,0.,H13.Bx,0.,0.)                                   # [T]
    
#grad_B = 0.0001                                    # [T/m]
#H13.fS = QuadEMfield()
#H13.fS.cxBz(grad_B)
#H13.fS.czBx(grad_B)


    






H13.bunch = b.getBunch(0)

"""
f = open('energies.txt','w')
for i in range(H13.bunch.getSize()):
    la = 1064.0e-9 
    n_states = 2
    E = H13.bunch.mass() + b.TK
    P = math.sqrt(E*E - H13.bunch.mass()*H13.bunch.mass())
    delta_E = 1./2. - 1./(2.*n_states*n_states)           
    la0 = 2*math.pi*5.291772108e-11/7.297352570e-3/delta_E
    te = b.TK - H13.bunch.mass()*(la/la0-1)
    kz = te/math.sqrt(P*P-te*te)
    
    px = H13.bunch.px(i)
    py = H13.bunch.py(i)
    pz = H13.bunch.pz(i)
    P = math.sqrt(px*px+py*py+pz*pz)
    E = math.sqrt(P*P+H13.bunch.mass()*H13.bunch.mass())
    cosb = (-px + pz*kz)/(math.fabs(P)*math.sqrt(1+kz*kz))
    beta = P/E
    gamma = 1./math.sqrt(1-beta*beta)
    w0 = gamma*(1 - beta*cosb)*2*math.pi*299792458*5.291772108e-11/la/2.187691263e6
    print >>f, w0
    
f.close()
"""

#H13.Nevol = 1000000
#H13.print_file = True

H13.TK = b.TK
H13.sigma_beam = b.sigma_beam

H13.n_sigma = 3
H13.n_step = 10000

H13.la = 532.0e-9                                   # [m]
H13.power = 150.0e6                               # [W]

H13.fx = -4.5                                      # [m]
H13.fy = -4.5                                      # [m]

H13.wx = 370.0e-6                               # [m]
H13.wy = 700.0e-6                               # [m]

    
H13.env_sigma = 21.0e-12                     #[m]



H13.rx = 2.2e-3                               # [m]
H13.ry = 1.5e-3                               # [m]

H13.ax = 1.0e-3                               # [rad]
H13.ay = 0.7e-3                               # [rad]

#----------------------End of the Laser field and excitation parameters of H0 atom----------------------#

#H13.bunch.dumpBunch("bunch_init1.dat")




#-------------------definition of the optimization function----------------------------#
#pf = printf("optimiz.dat","N_part","cpu_time", "W[MW]", "wx[um]", "wy[um]", "fx[cm]", "fy[cm]", "Population", "+- Err")
#pf = printf("optimiz_011_532.dat","N_step","cpu_time", "W[MW]","rx[mm]", "ry[mm]", "ax[mrad]", "ay[mrad]", "Population", "+- pop", "P_ioniz", "+- P_ioniz")
pf = printf("optimiz_2_1064.dat","N_step","cpu_time","sigma_laser[ps]","W[MW]", "rx[um]", "ry[um]", "ax[cm]", "ay[cm]", "Population", "+- pop", "P_ioniz", "+- P_ioniz")
#name_args, guess, increments = ['wx','wy'],[300.0e-6, 300.0e-6],[10e-6, 100e-6]
#name_args, guess, increments  = ['wx','wy','fx'], [333.0e-6, 333.0e-6,-2.800], [10e-6, 100e-6,1.00]
#name_args, guess, increments  = ['rx','ry','ax','ay'], [1.3e-3,7.7e-3,1.9e-3,1.9e-3], [1.0e-4,1.0e-4,1.0e-5,1.0e-5]


#name_args, guess, increments  = ['env_sigma','rx','ry','ax','ay'], [84.0e-12,2.5e-3,2.00e-3,0.7e-3,1.9e-3],[130e-13, 1.0e-4,1.0e-4,1.0e-4,1.0e-4]
name_args, guess, increments  = ['env_sigma','rx','ry','ax','ay'], [90.0e-12,2.6e-3,2.4e-3,0.3e-3,1.7e-3],[30e-13, 1.0e-4,1.0e-4,1.0e-4,1.0e-4]

#name_args, guess, increments  = ['env_sigma','rx','ry','ax','ay'], [91.0e-12,7.1e-3,2.00e-3,0.6e-3,1.3e-3],[30e-13, 1.0e-4,1.0e-4,1.0e-4,1.0e-4]
#name_args, guess, increments  = ['env_sigma','rx','ry','ax','ay'], [93.0e-12,7.1e-3,2.00e-3,0.3e-3,1.4e-3],[30e-13, 1.0e-4,1.0e-4,1.0e-4,1.0e-4]

#name_args, guess, increments  = ['rx','ry','ax','ay'], [68.45e-3,0.377e-3,0.2e-3,0.2e-3], [1.0e-4,1.0e-4,1.0e-4,1.0e-4]
#name_args, guess, increments  = ['wx','fx'], [370.0e-6, -4.500], [10e-6,1.00]
#name_args, guess, increments  = ['wx','wy','fx','fy'], [300.6e-6, 300.6e-6,-2.0,-2.0], [10e-6, 10e-6,1.00,1.00]



powers = [10.0e6, 20.0e6, 30.0e6, 40.0e6, 50.0e6, 60.0e6, 70.0e6, 80.0e6, 90.0e6, 100.0e6]
#powers = [0.5e6]
#powers = [0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.6e6, 0.7e6, 0.8e6, 0.9e6, 1.0e6, 1.1e6, 1.2e6, 1.3e6, 1.4e6, 1.5e6, 1.6e6, 1.7e6, 1.8e6, 1.9e6, 2.0e6]

def opt_func(args):

    for i in range(len(args)):
        H13.__dict__[name_args[i]] = args[i]

#    H13.ay = 0
        
    
    sum = 0
    for N_p in range(len(powers)):
        H13.power = powers[N_p]
        pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
        sum += pop
        pf.fdata("","",H13.power/1.0e+6,"","","","",pop,sigma_pop,p_ioniz, sigma_p_ioniz)
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






ne = 50
energies = []
for i in range(1,ne):
    energies.append(i*20.0e-5)

def opt_func3(args):

    for i in range(len(args)):
        H13.__dict__[name_args[i]] = args[i]
       
    
    sum = 0
    for i in range(ne-1):
        H13.power = energies[i]/(math.sqrt(2*math.pi)*H13.env_sigma)
#        H13.power = 3e6 
        pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
        sum += pop
        pf.fdata("","",H13.env_sigma*1e12,H13.power/1.0e+6,"","","","",pop,sigma_pop,p_ioniz, sigma_p_ioniz)
#        pf.fdata("","",H13.power/1.0e+6,"","","","",pop,sigma_pop)
    H13.count += 1
#    pf.fdata(H13.count,orbit_mpi.MPI_Wtime() - time_start,"aver_0.1-1.0",1.0e+6*H13.wx,1.0e+6*H13.wy,H13.fx*100,H13.fy*100,sum/len(powers))
    pf.fdata(H13.count,orbit_mpi.MPI_Wtime() - time_start,H13.env_sigma*1e12,"aver",1.0e+3*H13.rx,1.0e+3*H13.ry,1.0e+3*H13.ax,1.0e+3*H13.ay,sum/ne,sigma_pop,p_ioniz, sigma_p_ioniz)
    pf.fdata("")
    
    return -sum/(ne-1)


#-------------------definition of the optimization function----------------------------#












#-----------------------Beginning of optimization-----------------------------------------

s = Simplex(opt_func3, guess, increments)
(values, err, iter) = s.minimize(1e-20, 1000,0)

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
print pop, sigma_pop, p_ioniz, sigma_p_ioniz

#-----------------------End of Single point calculation with optimization function-----------------------------------------------


#-----------------------Emittance of stripped atoms-----------------------------------------------
"""
z1 = -0.220745063069
z2 = 0.220745063069
Nevol = H13.Nevol
sum = 0
for j in range(b.N_part):
    E = H13.bunch_target.mass() + H13.TK
    P = math.sqrt(E*E - H13.bunch_target.mass()*H13.bunch_target.mass()) 
    c_light = 2.99792458e+8
    beta = P/E
    gamma = 1./math.sqrt(1-beta*beta)
    
    
    ang = [0]*(Nevol+1)
    fd = [0]*(Nevol+1)
    dz = (z2 - z1)/Nevol
    prob = H13.bunch_target.partAttrValue("Evolution",j,Nevol)
        
    ang_aver = 0
    for i in range(Nevol-1,-1,-1):
        fd[i] = abs((H13.bunch_target.partAttrValue("Evolution",j,i+1) - H13.bunch_target.partAttrValue("Evolution",j,i))/dz/prob)
    
        ang[i] = ang[i+1] + H13.Bx*dz*c_light/(P*1.0e9)
        ang_aver += ang[i]*fd[i]*dz
              
          
    ang_aver2 = 0 
    for i in range(Nevol):
        ang_aver2 += (ang[i] - ang_aver)*(ang[i] - ang_aver)*fd[i]*dz
       
    
      
    print H13.Bx,1e3*math.sqrt(ang_aver2),ang_aver2*22.5*1e6, beta*gamma*ang_aver2*22.5*1e6,prob
    sum += beta*gamma*ang_aver2*22.5*1e6
print "aver=",sum/b.N_part

"""

    
#f = open('evol_out.txt','w')
#for i in range(Nevol):
#    print >>f, abs(H13.bunch_target.partAttrValue("Evolution",0,i))
#f.close()


#-----------------------Emittance of stripped atoms-----------------------------------------------




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


#-----------------------Slide show-----------------------------------------
"""
res_dir = os.environ["ORBIT_ROOT"]+"/ext/laserstripping/working_dir/"
data_name = "data"
pic_name = "pic"

for i in range(b.N_part):
     probi = H13.bunch_target.partAttrValue("Populations",i,0)
     PlotPopl([3,1],["P_ioniz= %1.6f"%probi],0.15,res_dir+data_name+"%i.dat"%(i*size+rank),res_dir+pic_name+"%i.png"%(i*size+rank))
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
