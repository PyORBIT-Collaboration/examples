#--------------------------------------------------------------------
# Translate ORBIT_MPI bunch to pyORBIT bunch.
# pyORBIT bunch needs the energy, mass and charge 
# of the synchronous particle, but the ORBIT_MPI does not have it. 
# So, this information is specified here in this script.
# The ring length should be defined in the input.
# ORBIT_MPI file has lines: x[mm] xp[mrad] y[mm] yp[mrad]   phi[rad]  dE[MeV]
# pyORBIT: x[m] xp[rad] y[m] yp[rad]  z[m]  dE[GeV]
#--------------------------------------------------------------------

kineticEnergy = 1.0 # kinetic energy in GeV

import sys
import math

from bunch import Bunch

if len(sys.argv) != 6 or sys.argv[1] != "L" or sys.argv[2] != "=" or float(sys.argv[3]) == None:
	print "Usage: ",sys.argv[0]," L = RingLength[m] <name of input ORBIT_MPI file>  <name of generated pyORBIT file>"
	sys.exit(1)

L = float(sys.argv[3]) 
print "Kinetic Energy [GeV] =",kineticEnergy," Ring Length[m]=",L

name_f_in = sys.argv[4]
name_f_out = sys.argv[5]

file_in = open(name_f_in,"r")

b = Bunch()
b.getSyncParticle().kinEnergy(kineticEnergy)

ln = file_in.readline()
count = 0
while ln:
	count += 1
	res_arr = ln.strip().split()
	val_arr = []
	for s in res_arr:
		val_arr.append(float(s))
	for i in range(4):
		val_arr[i] /= 1000.
	val_arr[4] = val_arr[4]*L/(2*math.pi)
	val_arr[5] = val_arr[5]/1000.0
	b.addParticle(val_arr[0],val_arr[1],val_arr[2],val_arr[3],val_arr[4],val_arr[5])
	ln = file_in.readline()
	
file_in.close()
print "Total number of macro-particles = ",count

b.dumpBunch(name_f_out)

print "Stop."

