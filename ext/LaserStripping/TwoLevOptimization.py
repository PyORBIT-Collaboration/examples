#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------

import sys
import math

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *

import os
import orbit_mpi
import random
from ext.las_str import part_generator,mygra

from orbit_mpi import mpi_comm,mpi_datatype,mpi_op

rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)

if(rank == 0): print "Start."

random.seed((rank+1)*12571)
file_name = "result_ls_3.dat"


z0=2.0e-3
delta_E=4./9.
dip_transition=math.sqrt(729./8192.)
TK = 1.0
#fx=fy=-0.1
la=355.0e-9

#wx=59.0e-6 
wy=850.0e-6
#power=5000000.
n_step = 500000
N_part = 1




alphaX = 0.                   # [rad]
betaX = 26.7                  # [m]
emtX = 0.225e-6               # [m*rad]

alphaY = 1.979            # [rad]
betaY = 0.75              # [m]
emtY = 0.225e-6           # [m*rad]

relativeSpread = 1.0e-4

dispD = 0.    # [m]
dispDP = 2.58 # [rad]

#range_power=[1.0e6,1.5e6,2.0e6]
#range_wx=[60.0e-6 ,90.0e-6, 120.0e-6, 150.0e-6, 180.0e-6, 210.0e-6]
#range_fxy=[-0.36,-0.40,-0.44,-0.48]
range_power=[1.0e6]
range_wx=[120.0e-6]
range_fxy=[-0.84]




cutOffX = math.sqrt(emtX*betaX)*3.0
trGenX = part_generator.TransverseCoordGen(alphaX,betaX,emtX,cutOffX)

cutOffY = math.sqrt(emtY*betaY)*3.0
trGenY = part_generator.TransverseCoordGen(alphaY,betaY,emtY,cutOffY)

gamaX = (1.0+alphaX*alphaX)/betaX
gamaY = (1.0+alphaY*alphaY)/betaY
sigmaXP = math.sqrt(emtX*gamaX)*1.0e+3
sigmaYP = math.sqrt(emtY*gamaY)*1.0e+3
if(rank == 0):
    print "cutOffX [mm]= %6.2f cutOffY [mm] = %6.2f "%(cutOffX*1.0e+3,cutOffY*1.0e+3)
    print "sigmaXP [mrad]= %7.4f sigmaYP [mrad] = %7.4f "%(sigmaXP,sigmaYP)

pGen = part_generator.EnergyGen(TK,relativeSpread)


partGen = part_generator.ParticlesGen(dispD,dispDP,trGenX,trGenY,pGen)

bunch_target = Bunch()
bunch_target.charge(0)
bunch_target.mass(0.938256 + 0.000511)
bunch_target.addPartAttr("Amplitudes",{"size":15})

E = bunch_target.mass() + TK
P = math.sqrt(E*E - bunch_target.mass()*bunch_target.mass())
vz=299792458*P/E
time_step = (2*z0/vz)/n_step




#print alpha
#print kz
#sys.exit(1)

if(rank == 0):
    f_out = open(file_name,"w")
    f_out.write("N  cpu_time   W [MW]     wx [um]   dist[cm]   Population +- Err \n")
    f_out.close()
    print "======START LOOP=========="

bunch = Bunch()
bunch.charge(0)
bunch.mass(0.938256 + 0.000511)
for i in range(N_part):
    (x,px,y,py,z,pz) = partGen.getCoords()
    #bunch.addParticle(x,px,y,py,-z0,pz)
    bunch.addParticle(0.,0.,0.,0.,-z0,pGen.getP0())

bunch.dumpBunch("bunch_init.dat")

bunch.addPartAttr("Amplitudes",{"size":15})    
for i in range(N_part):
    bunch.partAttrValue("Amplitudes",i,1,1.0)
    
time_start = orbit_mpi.MPI_Wtime()

count = 0

for power,wx,fxy in [(power,wx,fxy) 
                    for power in range_power 
                    for wx in range_wx 
                    for fxy in range_fxy]:

    count = count + 1
    bunch_target.deleteAllParticles()
    bunch.copyBunchTo(bunch_target)
    
    la0= 2*math.pi*5.291772108e-11/7.297352570e-3/delta_E

    kz=-1/math.sqrt(math.pow(P/(bunch_target.mass()*(la/la0-1)-TK),2)-1)

    #       alpha=360*math.acos((b.mass()*(la/la0-1)-TK)/P)/2/math.pi
    #       print alpha
    
    for i in range(bunch_target.getSize()):
        z = bunch_target.z(i)
        x = bunch_target.x(i)
        bunch_target.z(i,z - kz*x)
    bunch_target.dumpBunch("bunch_ini"+str(count)+".dat")
    fx=fy=fxy
    LFS=HermiteGaussianLFmode(math.sqrt(power),0,0,wx,wy,fx,fy,la) 
    LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   1.,0.,1./kz,  0.,1.,0.)
    tracker = RungeKuttaTracker(0.000000001)
    First = TwoLevelAtom(LFS,delta_E,dip_transition)
    fS=LSFieldSource(0.,0.,0.,0.,0.,0.)
    tracker.track(bunch_target,0,time_step*n_step, time_step,fS,First)
    bunch_target.dumpBunch("bunch_res"+str(count)+".dat")
    population = 0.
    population2 = 0    
    for i in range(bunch_target.getSize()):
        val = (1-bunch_target.partAttrValue("Amplitudes",0,1)*bunch_target.partAttrValue("Amplitudes",0,1)-bunch_target.partAttrValue("Amplitudes",0,2)*bunch_target.partAttrValue("Amplitudes",0,2))
        
        population += val
        population2 += val*val
    
    mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)
    op = mpi_op.MPI_SUM
    data_type = mpi_datatype.MPI_DOUBLE
    population = orbit_mpi.MPI_Allreduce(population,data_type,op,mpi_comm.MPI_COMM_WORLD)
    population2 = orbit_mpi.MPI_Allreduce(population2,data_type,op,mpi_comm.MPI_COMM_WORLD)
    population = population/(mpi_size*N_part)
    sigma_pop = 0.
    if(N_part*mpi_size > 1):
        sigma_pop = math.sqrt((population2 - N_part*population*population)/(N_part*mpi_size*(N_part*mpi_size - 1)))
    time_live = orbit_mpi.MPI_Wtime() - time_start
    res = " %6.0f  %4.1f  %4.1f  %4.1f    %6.3f  %6.3f "%(time_live,power/1.0e+6,1.0e+6*wx,fxy*100,population,sigma_pop)
    if(rank == 0): print  "W [MW]= %4.1f  wx [um] = %4.1f  dist[cm]= %4.1f  Population: %7.3f +-  %7.3f"%(power/1.0e+6,1.0e+6*wx,fxy*100,population,sigma_pop)
    if(rank == 0):
        f_out = open(file_name,"a")
        f_out.write(str(count)+" "+res + "\n")
        f_out.close()
        
        
#graph = mygra.PlotPopl(population)
#os.system('gthumb /home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/image.png')
#os.remove("/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt")
