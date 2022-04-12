#------------------------------------------------------
# This is an example of Loop Current Field Source 
# of magnetic field for 3D tracking using RK4_Tracker 
#-------------------------------------------------------

import sys
import math
import time

from orbit_utils import field_sources
from field_sources import LoopFieldSource

from orbit_utils import Matrix, PhaseVector
from spacecharge import Grid3D

loop_fs = LoopFieldSource()

#---- set radius of the loop
radius = 2.0
print "Radius  [m]= ",loop_fs.radius(radius)

#---- set current in the loop
current = 3.0
print "Current [A]= ",loop_fs.current(current)

#---- maximal z in [m]. Outside this limit the field will be 0
z_max = 1.001
print "Z max [m] =",loop_fs.maxZ(z_max)

#---- maximal r in [m]. Outside this limit the field will be 0
r_max = 1.5
print "R max [m] =",loop_fs.maxR(r_max)

#--- analytical field in the center of the loop in CI 
mu0 = 4*math.pi*1.0e-7
Bc = mu0*current/(2*radius)
print "B center [T] = ",Bc


def B_theory(z):
	mu0 = 4*math.pi*1.0e-7
	if(z == 0.): return mu0*current/(2*radius)
	return (mu0/2)*current*radius**2/(z**2 + radius**2)**1.5

#---------------------------------------
# Let's calculate field on the z-axis
#---------------------------------------
x = 0.
y = 0.
z_start = -1.2
z_end   = +1.2
z_step = 0.1
n_steps = int((z_end - z_start)/z_step) + 2
print "======================================================"
print "               Field along the z-axis"
print " z[m]      Bz[T]        BzTheory[T]    Relative Diff. "
for ind in range(n_steps):
	z = z_start + z_step*ind
	(Ex,Ey,Ez,Bx,By,Bz) = loop_fs.getFields(x,y,z)
	Bz_th = B_theory(z)
	st  = " %+5.3f "%z
	st += " %+12.5g  %+12.5g "%(Bz,Bz_th)
	st += "    %+15.8g "%((Bz_th-Bz)/Bz_th)
	print st
print "============================================="
print "         Let's check rot(B) = 0"
print "         Our results will not be perfect,"
print "         but we will see cancellations of "
print "         partial derivatives:"
print "         rot(B)x = dBz/dy - dBy/dz"
print "         rot(B)y = dBx/dz - dBz/dx"
print "         rot(B)z = dBy/dx - dBx/dy"
x =  0.2
y = -0.1

delta = 0.0001
delta2 = 2*delta
print " z[m]         (dBx/dx , dBy/dy , dBz/dz)   (dBxdz,dBxdy,dBydz,dBydx,dBzdx,dBzdy)   (rot(B)x,rot(B)y,rot(B)z) "
for ind in range(n_steps):
	z = z_start + z_step*ind
	(Ex,Ey,Ez,Bx1p,By1p,Bz1p) = loop_fs.getFields(x+delta,y,z)
	(Ex,Ey,Ez,Bx2p,By2p,Bz2p) = loop_fs.getFields(x,y+delta,z)
	(Ex,Ey,Ez,Bx3p,By3p,Bz3p) = loop_fs.getFields(x,y,z+delta)
	(Ex,Ey,Ez,Bx1m,By1m,Bz1m) = loop_fs.getFields(x-delta,y,z)
	(Ex,Ey,Ez,Bx2m,By2m,Bz2m) = loop_fs.getFields(x,y-delta,z)
	(Ex,Ey,Ez,Bx3m,By3m,Bz3m) = loop_fs.getFields(x,y,z-delta)
	rotX =  (Bz2p - Bz2m)/delta2  - (By3p - By3m)/delta2
	rotY = -((Bz1p - Bz1m)/delta2 - (Bx3p - Bx3m)/delta2)
	rotZ =  (By1p - By1m)/delta2  - (Bx2p - Bx2m)/delta2
	dBx = (Bx1p - Bx1m)/delta2
	dBy = (By2p - By2m)/delta2
	dBz = (Bz3p - Bz3m)/delta2
	dBxdz = (Bx3p - Bx3m)/delta2
	dBxdy = (Bx2p - Bx2m)/delta2
	dBydz = (By3p - By3m)/delta2
	dBydx = (By1p - By1m)/delta2
	dBzdx = (Bz1p - Bz1m)/delta2
	dBzdy = (Bz2p - Bz2m)/delta2	
	st = " %+5.3f "%z
	st += " (%+10.3g %+10.3g %+10.3g) "%(dBx,dBy,dBz)
	st += " (%+10.3g %+10.3g  %+10.3g %+10.3g  %+10.3g %+10.3g)"%(dBxdz,dBxdy,dBydz,dBydx,dBzdx,dBzdy)
	st += " (%+10.3g %+10.3g %+10.3g) "%(rotX,rotY,rotZ)
	print st
print "============================================="
print "Stop."
