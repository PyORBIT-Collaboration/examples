
import sys,math,os,orbit_mpi,random

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *

from orbit_mpi import mpi_comm,mpi_datatype,mpi_op
orbit_path = os.environ["ORBIT_ROOT"]

n_step = 10000

wx = 370.0e-6
wy = 700.0e-6

fx = -4.500
fy = -4.500

TK = 1


la = 355.0e-9                       
power = 0.1e6    
n_sigma = 3                          


n_states = 3
dip_transition = math.sqrt(256*math.pow(n_states,7)*math.pow(n_states-1,2*n_states-5)/3/math.pow(n_states+1,2*n_states+5))    
delta_E = 1./2. - 1./(2.*n_states*n_states)





#St = Stark(orbit_path+"/ext/laserstripping/transitions/",n_states)
fS = ConstEMfield(0.,0.,0.,0., 0.00001, 0.)

        



b = Bunch()
b.charge(0)
b.mass(0.938256 + 0.000511)
E = b.mass() + TK
P = math.sqrt(E*E - b.mass()*b.mass())

b.addParticle(0, 0, 0, 0, 0, P)
vz = 299792458*P/E   


la0 = 2*math.pi*5.291772108e-11/7.297352570e-3/delta_E
kz = -1/math.sqrt(math.pow(P/(b.mass()*(la/la0-1)-TK),2)-1.)    
#print math.atan(1./-kz)
z0 = n_sigma*wx*math.sqrt(1+math.pow(fx*la/(wx*wx*math.pi),2))*math.sqrt(1+kz*kz)

b.z(0,-z0 - kz*b.x(0))

time_step = (2*z0/vz)/n_step

LFS = HermiteGaussianLFmode(math.sqrt(power),0,0,abs(wx), abs(wy),fx,fy,la)
LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   1.,0.,1./kz,  0.,1.,0.)
tracker = RungeKuttaTracker(0)

#First = SchrodingerEquation(LFS,St,1000)#0.541633082102 0.544301583647
First = TwoLevelAtom(LFS,delta_E,dip_transition)#0.543385568089
#First = DensityMatrix(LFS,St,1000.) #0.541633082102

tracker.track(b,0,time_step*n_step, time_step,fS,First)

print 1 - b.partAttrValue("Populations",0,0) - b.partAttrValue("Populations",0,1)


