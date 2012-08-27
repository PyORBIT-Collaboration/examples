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


n_step = 2000
par = 6

if (method == 1):   Ez = 6.0e0
if (method == 2):   Ez = 6.0e7
if (method == 3):   Ez = 6.0e7
if (method == 4):   Ez = 5.14220642e11*0.0045
if (method == 5):   Ez = 5.14220642e11*0.0045


n_states = 3

#AttrValue= sum= 0.439180474444 0.439547597955
data_name = "data_pop_"
pic_name = "imageFS.png"

ta = 2.418884326505e-17                             # atomic unit of time
Rabi = 1e+12*ta                                      # Rabi frequency
Gamma = -math.pi*Rabi*Rabi/2/math.log(1-0.87)*1.0        # frequency sweep 
en_trans = 1./2. - 1./(2.*n_states*n_states)              # energy transition for two level atom
#en_trans = 0.5 - 0.12827056290270467212385977121
#en_trans = 0.3 - 0.12827032524816808378540833481

Omega = en_trans                                     # basis frequency
dip_trans = math.sqrt(256*math.pow(n_states,7)*math.pow(n_states-1,2*n_states-5)/3/math.pow(n_states+1,2*n_states+5))   # dipole transition between two levels
Elas=Rabi/dip_trans                                     # Amplitude of laser field

#----------------------End of the tracker and laser parameters----------------------#



fS = ConstEMfield(0.,0.,Ez,0.,0.,0.)

b = Bunch()
b.charge(0)
b.mass(0.938256 + 0.000511)
b.addParticle(0.,0.,0.,0.,0.,0.)



if (method == 2 or method == 3):     
    Stark_ef = Stark(orbit_path+"/ext/laserstripping/Hydrogen_data/StarkEG/", n_states)


if (method == 4):     
    Stark_ef = StarkStrongField(orbit_path+"/ext/laserstripping/Hydrogen_data/StarkEG/",0,0,1)
    Omega = Stark_ef.deltaE(b.mass(),0.,0.,Ez,0.,0.,0.,0.,0.,0.)
#    Omega = 0.5 - 0.1275
#    Omega += 2.0e-4

if (method == 5):     
    continuum_spectr = TDMcontinuum(orbit_path + "/ext/laserstripping/Hydrogen_data/TDMcontinuum/001/")
    Omega = continuum_spectr.setField_returndE(b.mass(),0.,0.,Ez,0.,0.,0.,0.,0.,0.)
#    Omega = 0.5 - 0.12775
#    Omega += 2.0e-4




LFS = FroissartStoraLF(Omega,Gamma,Elas)

if (method == 2 or method == 3):    LFS.setLaserFieldPolarization(1.,1.,1.)
if (method == 4):    LFS.setLaserFieldPolarization(1.,1.,0.)
if (method == 5):    LFS.setLaserFieldPolarization(1.,1.,1.) #    It is supposed that polarization of the laser field has optimal direction for this method


if (method == 1):   eff = TwoLevelAtom(LFS,en_trans,dip_trans)
if (method == 2):   eff = SchrodingerEquation(LFS,Stark_ef,200.) 
if (method == 3):   eff = DensityMatrix(LFS,Stark_ef,200.)
if (method == 4):   eff = TwoLevelStrongField(LFS, Stark_ef)  
if (method == 5):   eff = ContinuumSS(LFS,continuum_spectr)

nevo = 2000
pr = PrintExtEffects("Populations",nevo,addr+data_name)
evo = RecordEvolution("Populations",0,nevo)

cont_eff = ExtEffectsContainer()
cont_eff.AddEffect(evo)
cont_eff.AddEffect(eff)
if(method != 5):
    cont_eff.AddEffect(pr)



tracker = RungeKuttaTracker(0)
time_step=(2*math.pi*ta/Rabi)/n_step

#b.addPartAttr("Amplitudes",{"size":5})
#b.partAttrValue("Amplitudes",0,1,1.0)

tracker.track(b,-par*time_step*n_step,2*par*time_step*n_step, time_step,fS,cont_eff)



pop1 = b.partAttrValue("Populations",0,1)
sum = 1 - b.partAttrValue("Populations",0,0)

pop2 = sum - pop1


if (method == 1):   ratio = [3,1]
if (method == 2):   ratio = [5,3]
if (method == 3):   ratio = [5,3]
if (method == 4):   ratio = [3,1]
if (method == 5):   ratio = [1,1]

if (method == 5):
    
    size = b.getPartAttrSize("Populations") - 1
    f = open(addr + "spectrum.txt",'a')
    for i in range(1,size):
        print >>f, b.partAttrValue("Populations",0,i)
    f.close()
    
    size = b.getPartAttrSize("Evolution")
    f = open(addr + data_name+"0.dat",'a')
    for i in range(0,size -4):
        print >>f, -par*time_step*n_step + i*2*par*time_step*n_step/nevo, b.partAttrValue("Evolution",0,i)
    f.close()
    
    
    

PlotPopl(ratio,["%1.4f"%pop2],0.15,addr+data_name+"0.dat",addr+pic_name)

os.system('eog '+addr+pic_name)
os.remove(addr + data_name+"0.dat")

if (method == 5):
    os.remove(addr + "spectrum.txt")
    
print "AttrValue=","sum=", pop2, sum

n = n_states

#for n1 in range(n):
#    n2 = n-n1-1
#    print 1+n1+(n*n*n-n)/3,"  ", b.partAttrValue("Populations",0,1+n1+(n*n*n-n)/3)/0.87
    
#for n1 in range(n-1):
#    n2 = n-n1-2
#    print n
#    print 1+n+n1+(n*n*n-n)/3,"  ", b.partAttrValue("Populations",0,1+n+n1+(n*n*n-n)/3)/0.87,"  ",3.*(n1+1)*(n2+1)/(n*(n*n-1))
#    print 2-n+n1+(n*n*n-n)/3,"  ", b.partAttrValue("Populations",0,2-n+n1+(n*n*n-n)/3)/0.87,"  ",3.*(n1+1)*(n2+1)/(n*(n*n-1))


#b.dumpBunch("evol.dat")

#g = Bunch()
#g.readBunch("evol.dat")

#AttrValue= sum= 0.882746852119 1.00000000001
#AttrValue= sum= 0.88274684111  1.00000000001
#AttrValue= sum= 0.878141409529 1.0

