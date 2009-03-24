import sys, math, os, orbit_mpi

from bunch import *
from ext.las_str.TwoLevelFuncMod import TwoLevelFunc
from ext.las_str.part_generator import BunchGen
from ext.las_str.plot_mod import *










res_dir = os.environ["ORBIT_ROOT"]+"/ext/laserstripping/working_dir/"

#----------------------Beginning of the beam parameters----------------------#

b = BunchGen()

b.TK= 1.0                                     # [GeV]
b.N_part = 10
b.N_attr = 2*2+8
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

bunch_init = b.getBunch(0)
bunch_init.copyBunchTo(H13.bunch)

sum = 0
for i in range(b.N_part):
    H13.bunch.deleteAllParticles()
    (x0, px0, y0, py0, z0, pz0) = (bunch_init.x(i),bunch_init.px(i),bunch_init.y(i),bunch_init.py(i),bunch_init.z(i),bunch_init.pz(i))
#    (x0, px0, y0, py0, z0, pz0) = ...
    H13.bunch.addParticle(x0, px0, y0, py0, z0, pz0)
    H13.bunch.partAttrValue(b.attr_name,0,1,1.0)
    H13.data_addr_name = res_dir+"/data%i.dat"%i
    pop, sigma_pop = H13.population()
    sum += pop
    graph = PlotPopl([3,1],["%1.6f"%pop,"px= %1.6f"%px0,"y= %1.6f"%y0],0.15,res_dir+"/data%i.dat"%i,res_dir+"/pic%i.png"%i)
    os.remove(H13.data_addr_name)
    
print "population = ", sum/b.N_part













    
