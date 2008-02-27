from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch

print "Start."

b = Bunch()
b.addParticle(1.0e-3,0.0,0.0,0.0,0.0,0.0)
b.addParticle(0.0,1.0e-3,0.0,0.0,0.0,0.0)
b.addParticle(0.0,0.0,1.0e-3,0.0,0.0,0.0)
b.addParticle(0.0,0.0,0.0,1.0e-3,0.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,1.0e-3,0.0)
b.addParticle(0.0,0.0,0.0,0.0,0.0,1.0e-3)
b.compress()

syncPart = b.getSyncParticle()
#energy in GeV
energy = 1.0                          
syncPart.kinEnergy(energy)

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("sext_623_620_00.mad","RNG")
print "Track Bunch."
teapot_latt.trackBunch(b)
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()
print "Stop."
b.dumpBunch()

#---------------Results to compare -----------------------------
#Read MAD.
#Track Bunch.
#Lattice= RNG  length [m] = 248.0
#Stop.
#% PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES 
#% BUNCH_ATTRIBUTE_DOUBLE charge   1 
#% BUNCH_ATTRIBUTE_DOUBLE classical_radius   1.5347e-18 
#% BUNCH_ATTRIBUTE_DOUBLE macro_size   0 
#% BUNCH_ATTRIBUTE_DOUBLE mass   0.938272 
#%  SYNC_PART_COORDS 0 0 0  x, y, z positions in [m]
#%  SYNC_PART_MOMENTUM 0 0 1.696037912  px, py, pz momentum component in GeV/c
#%  info only: energy of the synchronous particle [GeV] = 1 
#%  info only: momentum of the synchronous particle [GeV/c] = 1.696037912 
#%  info only: beta=v/c of the synchronous particle = 0.8750256155 
#%  info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle = 2.065788684 
#%  SYNC_PART_TIME 9.453882737e-07  time in [sec]
#% x[m] px[rad] y[m] py[rad] z[m]  (pz or dE [GeV]) 
#0.00096526633 -5.4566491e-05 0 0 2.0623749e-05 0 
#0.00015595394 0.0010274038 0 0 5.6234025e-05 0 
#4.2191444e-10 4.9129255e-10 -0.0011433267 0.00019585874 2.7585565e-06 0 
#-3.5851218e-08 1.203286e-07 -0.02593513 0.0035800767 0.00074366546 0 
#0 0 0 0 0.001 0 
#2.3315706e-07 1.3078068e-07 0 0 -0.03300023 0.001 


