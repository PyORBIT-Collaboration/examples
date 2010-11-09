#-----------------------------------------------------
#UniformEllipsoidFieldCalculator test
#-----------------------------------------------------

import sys
import math
import time

from spacecharge import UniformEllipsoidFieldCalculator


print "Start."

#ellipsoid parameters
a = 1.0
b = 1.0
c = 1.0
r_max = 10.
uniformEllips = UniformEllipsoidFieldCalculator()
print "Start."
uniformEllips.setEllipsoid(a,b,c,r_max)

(dirX,dirY,dirZ) = (1.5,2.5,0.3)
mod_dir = math.sqrt(dirX**2 + dirY**2 + dirZ**2)
(dirX,dirY,dirZ) =(dirX/mod_dir,dirY/mod_dir,dirZ/mod_dir)

nPoints = 100
r_step = 11./(nPoints-1)

max_diff = 0.
max_diff_point = (0.,0.,0.)
for ir in range(nPoints):
	r = r_step*(ir+1)
	(x,y,z) = (dirX*r,dirY*r,dirZ*r)
	(ex,ey,ez) = uniformEllips.calcField(x,y,z)
	if(r < 1.0):
		(ex_th,ey_th,ez_th) = (x,y,z)
	else:
		(ex_th,ey_th,ez_th) = (dirX/r**2,dirY/r**2,dirZ/r**2)
	diff = (ex-ex_th)**2 + (ey-ey_th)**2 + (ez-ez_th)**2 
	mod_e = ex*ex + ey*ey + ez*ez
	diff = 100*math.sqrt(diff/mod_e)
	if(diff > max_diff): 
		max_diff=diff
		max_diff_point = (x,y,z)

print "max difference between theory and calculations %=",max_diff
print "at point (x,y,z)=",max_diff_point

print "===================point analysis==========="
#calculate field
(x,y,z) = max_diff_point
r = math.sqrt(x**2+y**2+z**2)
(ex,ey,ez) = uniformEllips.calcField(x,y,z)
print "(x,y,z)=",(x,y,z)
print "(ex,ey,ez)=",(ex,ey,ez)
print "abs(ex,ey,ez)=",math.sqrt(ex**2+ey**2+ez**2)
if(r > 1.0):
	(ex_th,ey_th,ez_th) = (x/r**3,y/r**3,z/r**3)
	print "(ex_th,ey_th,ez_th)=",(ex_th,ey_th,ez_th)
	print "abs(ex_th,ey_th,ez_th)=",math.sqrt(ex_th**2+ey_th**2+ez_th**2)
else:
	(ex_th,ey_th,ez_th) = (x,y,z)
	print "(ex_th,ey_th,ez_th)=",(ex_th,ey_th,ez_th)
	print "abs(ex_th,ey_th,ez_th)=",math.sqrt(ex_th**2+ey_th**2+ez_th**2)

"""
#--------------Speed test -----------
count = 0
count_step = 1000
time_start = time.clock()
while(1 < 2):
	uniformEllips.setEllipsoid(a,b,c,r_max)
	count += 1
	if(count % count_step == 0):
		time_stop = (time.clock() - time_start)/count_step
		print "i=",count," freq=",time_stop
		time_start = time.clock()
"""

print "Stop."

