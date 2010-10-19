#--------------------------------------------------------------------
# Translate pyORBIT bunch to ORBIT_MPI bunch.
# The ring length should be defined in the input.
# ORBIT_MPI file has lines: x[mm] xp[mrad] y[mm] yp[mrad]   phi[rad]  dE[MeV]
# pyORBIT: x[m] xp[rad] y[m] yp[rad]  z[m]  dE[GeV]
#--------------------------------------------------------------------

import sys
import math

from bunch import Bunch

if len(sys.argv) != 6 or sys.argv[1] != "L" or sys.argv[2] != "=" or float(sys.argv[3]) == None:
	print "Usage: ",sys.argv[0]," L = RingLength[m]  <name of input pyORBIT file> <name of output ORBIT_MPI file>"
	sys.exit(1)

L = float(sys.argv[3]) 
print "Ring Length[m]=",L

name_f_in = sys.argv[4]
name_f_out = sys.argv[5]

b = Bunch()
b.readBunch(name_f_in)

file_out = open(name_f_out,"w")

for i in range(b.getSize()):
	x = b.x(i)*1000.
	px = b.px(i)*1000.
	y = b.y(i)*1000.
	py = b.py(i)*1000.
	z = b.z(i)*2.0*math.pi/L
	dE = b.dE(i)*1000.
	file_out.write(str(x) + " " + str(px) + " " + str(y) + " " + str(py) + " "+ str(z) + " " + str(dE) + "\n")

file_out.close()
print "Total number of macro-particles = ",b.getSize()

print "Stop."
