#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------

import sys,math,os

from bunch import *
from trackerrk4 import *
from laserstripping import *
from ext.las_str.plot_mod import *







method = 2

#----------------------Beginning of the tracker and laser parameters----------------------#
orbit_path = os.environ["ORBIT_ROOT"]
addr = orbit_path+"/ext/laserstripping/working_dir/"
trans = orbit_path+"/ext/laserstripping/transitions/"

n_step = 1000
par=10

Ex=4.0e-7

n_states = 3


data_name = "data_ampl_"
pic_name = "imageFS.png"

ta = 2.418884326505e-17                             # atomic unit of time
Rabi = 1e+12*ta                                      # Rabi frequency
Gamma = -math.pi*Rabi*Rabi/2/math.log(1-0.87)         # frequency sweep 
en_trans = 1./2. - 1./(2.*n_states*n_states)              # energy transition for two level atom
Omega = en_trans                                     # basis frequency
dip_trans = math.sqrt(256*math.pow(n_states,7)*math.pow(n_states-1,2*n_states-5)/3/math.pow(n_states+1,2*n_states+5))   # dipole transition between two levels
Elas=Rabi/dip_trans                                     # Amplitude of laser field

#----------------------End of the tracker and laser parameters----------------------#


LFS = FroissartStoraLF(Omega,Gamma,Elas) 
LFS.setLaserFieldPolarization(1.,1.,1.)
fS = ConstEMfield(Ex,0.,0.,0.,0.,0.)

b = Bunch()
b.charge(0)
b.mass(0.938256 + 0.000511)
b.addParticle(0.,0.,0.,0.,0.,0.)


if(method == 2 or method == 3):     Stark_ef = Stark(trans, n_states)

if (method == 1):   eff = TwoLevelAtom(LFS,en_trans,dip_trans)
if (method == 2):   eff = SchrodingerEquation(LFS,Stark_ef,1000.) 
if (method == 3):   eff = DensityMatrix(LFS,Stark_ef,1000.)  

pr = PrintExtEffects("Populations",10000,addr+data_name)
evo = RecordEvolution("Populations",1,200)

cont_eff = ExtEffectsContainer()
cont_eff.AddEffect(evo)
cont_eff.AddEffect(pr)
cont_eff.AddEffect(eff)


tracker = RungeKuttaTracker(0)
time_step=(2*math.pi*ta/Rabi)/n_step

tracker.track(b,-par*time_step*n_step,2*par*time_step*n_step, time_step,fS,cont_eff)

pop1 = b.partAttrValue("Populations",0,1)
sum = 1 - b.partAttrValue("Populations",0,0)

pop2 = sum - pop1


if (method == 1):                  ratio = [3,1]
if (method == 2 or method == 3):   ratio = [5,3]

PlotPopl(ratio,["%1.4f"%pop2],0.15,addr+data_name+"0.dat",addr+pic_name)
os.system('eog '+addr+pic_name)
os.remove(addr+data_name+"0.dat")
print "AttrValue=","sum=", pop2, sum

#b.dumpBunch("evol.dat")

#g = Bunch()
#g.readBunch("evol.dat")



