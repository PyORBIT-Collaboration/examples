import sys
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement, MAD_LattLine
from orbit.time_dep import time_dep
from orbit.time_dep.waveform import KickerWaveform, MagnetWaveform, LinearMagnetWaveform
from bunch import Bunch

b = Bunch()
b.addParticle(1.0e-3,0.0,0.0,0.0,0.0,0.0)
b.addParticle(0.0,1.0e-3,0.0,0.0,0.0,0.0)
b.addParticle(0.0,0.0,1.0e-3,0.0,0.0,0.0)
b.addParticle(0.0,0.0,0.0,1.0e-3,0.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,1.0e-3,0.0)
b.addParticle(0.0,0.0,0.0,0.0,0.0,1.0e-3)
b.compress()
syncPart = b.getSyncParticle()
energy = 1.0                          
syncPart.kinEnergy(energy)

time_dep_latt = time_dep.TIME_DEP_Lattice()
time_dep_latt.readMAD("./LATTICES/rcs1124.dat","RCS")
#time_dep_latt.readMAD("./LATTICES/RTBT_00.LAT","TARG1")
time_dep_latt.setLatticeOrder()

WaveForm01 = LinearMagnetWaveform()
WaveForm01.initialize(0,1.0,1.0,0.97)

time_dep_latt.setTimeDepNode("R1SF02_1",WaveForm01)
time_dep_latt.setTurns(10)
time_dep_latt.trackBunchTurns(b)

