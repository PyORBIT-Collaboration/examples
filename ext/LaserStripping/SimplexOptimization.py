#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------

import sys, math, os, orbit_mpi

from bunch import *
from orbit_mpi import *
from orbit_utils import *
from ext.las_str.TwoLevelFuncMod import TwoLevelFunc
from ext.las_str.SimplexMod import Simplex
from ext.las_str.part_generator import BunchGen
from ext.las_str.print_mod import printf











file_name = "results_1e5.dat"

#----------------------Beginning of the beam parameters----------------------#

b = BunchGen()

b.TK= 1.0                                     # [GeV]
b.N_part = 10
b.N_attr = 15
b.attr_name = "Amplitudes"

b.mass = 0.938256 + 0.000511                  # [GeV]
b.charge = 0                                  # [e]

b.alphaX = 0.                                 # [rad]
b.betaX = 26.7                                # [m]
b.emtX = 0.225e-6                             # [m*rad]
b.cutOffX = math.sqrt(b.emtX*b.betaX)*3.0     # [m]

b.alphaY = 1.979                              # [rad]
b.betaY =0.75                                 # [m]
b.emtY = 0.225e-6                             # [m*rad]
b.cutOffY = math.sqrt(b.emtY*b.betaY)*3.0     # [m]

b.relativeSpread = 0.5e-3

b.dispD = 0.                                  # [m]
b.dispDP = 2.6                                # [rad]

#----------------------End of the beam parameters----------------------#

#----------------------Beginning of the Laser field and excitation parameters of H0 atom----------------------#

H13 = TwoLevelFunc()




H13.bunch = b.getBunch(0)
H13.TK = b.TK

H13.n_sigma = 3
H13.n_step = 10000

H13.la = 355.0e-9                               # [m]
H13.power = 1.0e6                               # [W]

H13.dip_transition = math.sqrt(729./8192.)      # [a.u]
H13.delta_E = 4./9.                             # [a.u]

H13.fx = -0.2                                   # [m]
H13.fy = -0.2                                   # [m]

H13.wx = 87.8e-6                                # [m]
H13.wy = 1091.6e-6                              # [m]




#----------------------End of the Laser field and excitation parameters of H0 atom----------------------#
#H13.bunch.dumpBunch("bunch_init.dat")




#-------------------definition of the optimization function----------------------------#

#name_args, guess, increments = ['wx','wy'],[200.0e-6, 2000.0e-6],[10e-6, 100e-6]
#name_args, guess, increments  = ['wx','wy','fx'], [100.0e-6, 1000.0e-6,-1.000], [10e-6, 100e-6,1.00]
name_args, guess, increments  = ['wx','wy','fx','fy'], [323.6e-6, 884.6e-6,-6.907,-6.907], [10e-6, 100e-6,1.00,1.00]
time_start = orbit_mpi.MPI_Wtime()
pf = printf(file_name,"N_part","cpu_time", "W[MW]", "wx[um]", "wy[um]", "fx[cm]", "fy[cm]", "Population", "+- Err")
powers = [1.0e6,2.0e6,3.0e6,4.0e6,5.0e6,6.0e6,7.0e6,8.0e6,9.0e6,10.0e6]


H13.count  = 0
def opt_func(args):

    for i in range(len(args)):
        H13.__dict__[name_args[i]] = args[i]
        
#        H13.fy = H13.fx

    sum = 0
    for N_p in range(len(powers)):
        H13.power = powers[N_p]
        pop, sigma_pop = H13.population()
        sum += pop
        pf.fdata("","",H13.power/1.0e+6,"","","","",pop,sigma_pop)
    H13.count += 1
    pf.fdata(H13.count,orbit_mpi.MPI_Wtime() - time_start,"aver_1_10",1.0e+6*H13.wx,1.0e+6*H13.wy,H13.fx*100,H13.fy*100,sum/len(powers))

    
    return -sum/len(powers)
#-------------------definition of the optimization function----------------------------#








s = Simplex(opt_func, guess, increments)
(values, err, iter) = s.minimize(1e-20, 1000,0)












#pf = printf(file_name,"dispDP", "relativeSpread", "population")
#pf = printf(file_name,"power [MW]", "beta_X [m]", "population", "+- Err")
"""
for spread, disp  in [(spread,disp)
                        for spread in [0.1e-3,0.2e-3,0.3e-3,0.4e-3,0.5e-3,0.6e-3,0.7e-3,0.8e-3,0.9e-3,1.0e-3,1.1e-3,1.2e-3,1.3e-3,1.4e-3,1.5e-3,1.6e-3,1.7e-3,1.8e-3,1.9e-3,2.0e-3]
                        for disp in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0]]:
"""
"""
for beta_X, power  in [(beta_X,power)
                        for beta_X in [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
                        for power in [1.0e6,2.0e6,3.0e6,4.0e6,5.0e6,6.0e6,7.0e6,8.0e6,9.0e6,10.0e6]]:
    b.betaX = beta_X
    H13.power = power
    H13.bunch = b.getBunch(0)
    pop, sigma_pop = H13.population()
    


    pf.fdata(H13.power/1.0e6,b.betaX,angleX,pop,sigma_pop)

"""

#popul = population_average_powers()
#popul, sigma = H13.population()



