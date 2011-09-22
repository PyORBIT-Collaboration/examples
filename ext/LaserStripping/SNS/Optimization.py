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
from ext.las_str.emittance import Freq_spread


time_start = orbit_mpi.MPI_Wtime()
rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)
orbit_path = os.environ["ORBIT_ROOT"]





#----------------------Beginning of the beam parameters----------------------#

b = BunchGen()

b.N_part = 250

b.TK= 3.0                                    # [GeV]

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
b.dispDP = -2.6                                # [rad]

b.sigma_beam = 30e-12/(2*math.sqrt(2*math.log(2)))                          #[m]
b.cutOffZ = 3*b.sigma_beam



#----------------------End of the beam parameters----------------------#

#----------------------Beginning of the Laser field and excitation parameters of H0 atom----------------------#


H13 = SchredingerFunc(0)

H13.Bx = 0.001






#H13.Bx = 0.005e-3    
#H13.fS = ConstEMfield(0.,0.,0.,H13.Bx,0.,0.)                                   # [T]
    
#grad_B = 1.0                                    # [T/m]
#H13.fS = QuadEMfield()
#H13.fS.cxBz(grad_B)
#H13.fS.czBx(grad_B)








H13.bunch = b.getBunch(0)
sys.exit(0)
#H13.Nevol = 1000
#H13.print_file = True

H13.TK = b.TK
H13.sigma_beam = b.sigma_beam

H13.n_sigma = 3
H13.n_step = 100000

H13.la = 355.0e-9                               # [m]
H13.power = 2.2e6                               # [W]

H13.fx = -4.5                                      # [m]
H13.fy = -4.5                                      # [m]

H13.wx = 370.0e-6                               # [m]
H13.wy = 700.0e-6                               # [m]

    
H13.env_sigma = 55e-12/(2*math.sqrt(2*math.log(2)))  #[m]



H13.rx = 0.8580e-3                               # [m]
H13.ry = 0.8580e-3                               # [m]

H13.ax = 0.3442e-3                               # [rad]
H13.ay = 0.3442e-3                               # [rad]

#----------------------End of the Laser field and excitation parameters of H0 atom----------------------#

#H13.bunch.dumpBunch("bunch_init1.dat")




#-------------------list of parameters to be optimized----------------------------#

name_args,guess,increments,name_args_pr,print_factor = [],[],[],[],[]

#name_args.append('env_sigma'),  guess.append(15.0e-12), increments.append(24.0e-13),name_args_pr.append('env_sigma[ps]'),print_factor.append(1e12)
name_args.append('rx'),         guess.append(1.69e-3),   increments.append(1.0e-4),  name_args_pr.append('rx[mm]'),print_factor.append(1e3)
#name_args.append('ry'),         guess.append(0.77e-3),   increments.append(1.0e-4),  name_args_pr.append('ry[mm]'),print_factor.append(1e3)
name_args.append('ax'),         guess.append(0.36e-3),   increments.append(1.0e-4),  name_args_pr.append('ax[mrad]'),print_factor.append(1e3)
#name_args.append('ay'),         guess.append(0.24e-3),   increments.append(1.0e-4),  name_args_pr.append('ay[mrad]'),print_factor.append(1e3)
#name_args.append('Bx'),         guess.append(2.00),      increments.append(0.1),     name_args_pr.append('Bx[T]'),print_factor.append(1)

pf = printf("test.dat",["N_step","cpu_time","W[MW]"] + name_args_pr + ["Popul", "+- Pop", "P_ioniz", "+- P_ioniz"])
#-------------------list of parameters to be optimized----------------------------#


#powers = [0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.6e6, 0.7e6, 0.8e6, 0.9e6, 1.0e6]
#powers = [0.5e6]
#powers = [0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.6e6, 0.7e6, 0.8e6, 0.9e6, 1.0e6, 1.1e6, 1.2e6, 1.3e6, 1.4e6, 1.5e6, 1.6e6, 1.7e6, 1.8e6, 1.9e6, 2.0e6]

ne = 23
powers = []
for i in range(1,ne):
    powers.append(i*0.1e6)

def opt_func(args):
    
    em = []
    for i in range(len(args)):
        H13.__dict__[name_args[i]] = args[i]
        em.append("")
       
    H13.ay = H13.ax
    H13.ry = H13.rx
    
    sum_ioniz = 0
    sum_pop2 = 0
    for i in range(ne-1):
        H13.power = powers[i]
#        H13.power = 3e6 
        pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
        sum_ioniz += p_ioniz
        sum_pop2 += pop
        
        pf.fdata(["","",H13.power/1.0e+6] + em + [pop,sigma_pop,p_ioniz, sigma_p_ioniz])

    H13.count += 1


    pr_args = []
    for i in range(len(args)):
        pr_args.append(print_factor[i]*H13.__dict__[name_args[i]])
    pf.fdata([H13.count,orbit_mpi.MPI_Wtime() - time_start,"aver"] + pr_args + [sum_pop2/(ne - 1),sigma_pop,sum_ioniz/(ne - 1), sigma_p_ioniz])
    pf.fdata([""])
    
#    return -sum_ioniz/(ne-1)
    return -sum_pop2/(ne-1)

#-------------------definition of the optimization function----------------------------#



"""
ne = 20
energies = []
for i in range(1,ne):
    energies.append(i*1.0e-5)

def opt_func3(args):
    
    em = []
    for i in range(len(args)):
        H13.__dict__[name_args[i]] = args[i]
        em.append("")
       
    
    sum_ioniz = 0
    sum_pop2 = 0
    for i in range(ne-1):
        H13.power = energies[i]/(math.sqrt(2*math.pi)*H13.env_sigma)
#        H13.power = 3e6 
        pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
        sum_ioniz += p_ioniz
        sum_pop2 += pop
        
        pf.fdata(["","",H13.power/1.0e+6] + em + [pop,sigma_pop,p_ioniz, sigma_p_ioniz])

    H13.count += 1


    pr_args = []
    for i in range(len(args)):
        pr_args.append(print_factor[i]*H13.__dict__[name_args[i]])
    pf.fdata([H13.count,orbit_mpi.MPI_Wtime() - time_start,"aver"] + pr_args + [sum_pop2/(ne - 1),sigma_pop,sum_ioniz/(ne - 1), sigma_p_ioniz])
    pf.fdata([""])
    
    return -sum_ioniz/(ne-1)
#    return -sum_pop2/(ne-1)

"""
#-------------------definition of the optimization function----------------------------#




#-----------------------Beginning of optimization-----------------------------------------

#s = Simplex(opt_func, guess, increments)
#(values, err, iter) = s.minimize(1e-20, 1000,0)

#-----------------------End of Optimization-----------------------------------------------









#-----------------------Single point calculation-----------------------------------------------------
#print "start"
#print  H13.population() 


for i in range(0,50):
    H13.Bx = i*0.001+1.0e-10       

    pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
    if(rank==0):    print i*0.001,"  ",pop,"  ",sigma_pop

#-----------------------End of single point calculation-----------------------------------------------





#-----------------------Single point calculation with optimization function-----------------------------------------------------

#pop, sigma_pop, p_ioniz, sigma_p_ioniz = H13.population()
#print pop, sigma_pop, p_ioniz, sigma_p_ioniz

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
