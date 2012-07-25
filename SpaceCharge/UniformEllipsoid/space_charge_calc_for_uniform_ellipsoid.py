#-----------------------------------------------------
#Space Charge Calculator for Uniform Ellipsoid test 
#If a=b=c=1 and Q=1 it is an uniformly charged sphere 
#Inside the field is E=(x,y,z), outside E(r)=(x,y,z)/r^3
#-----------------------------------------------------

import sys
import math
import time
import random

from bunch import Bunch

from spacecharge import UniformEllipsoidFieldCalculator
from spacecharge import SpaceChargeCalcUnifEllipse

print "Start."

nEllipses = 2
spaceChargeCalc = SpaceChargeCalcUnifEllipse(2)

b = Bunch()
syncPart = b.getSyncParticle()
syncPart.kinEnergy(0.0025)

gamma = syncPart.gamma()
print "gamma=",gamma

nParticles = 100000
macrosize = 1./nParticles
b.macroSize(macrosize)
a_ellipse  = 1.0
b_ellipse  = 1.0
c_ellipse  = 1.0
r_max = max(a_ellipse,b_ellipse,c_ellipse) 
x = 0.
y = 0.
z = 0.
for i in range(nParticles):
	while(1 < 2):
		x = 2.*a_ellipse*(random.random() - 0.5)
		y = 2.*b_ellipse*(random.random() - 0.5)
		z = 2.*c_ellipse*(random.random() - 0.5)
		if((x**2/a_ellipse**2 + y**2/b_ellipse**2 + z**2/b_ellipse**2) < 1.):
			break
	b.addParticle(x,0.,y,0.,z,0.)

time_start = time.clock()

spaceChargeCalc.trackBunch(b, 1.0)

print "time =",(time.clock() - time_start)

#-------------------------------------------------	
#this is the example of using the Gnuplot package
#Comment this part if you do not have this package installed
#-------------------------------------------------
r_max_graph = 2*r_max
nGraphPoints = 1000


filed_x_arr = []
for i in range(nGraphPoints):
	r = i*2*r_max_graph/(nGraphPoints) - r_max_graph
	(ex,ey,ez) = spaceChargeCalc.calculateField(r,0.,0.)
	filed_x_arr.append((r,ex))
	
filed_y_arr = []
for i in range(nGraphPoints):
	r = i*2*r_max_graph/(nGraphPoints) - r_max_graph
	(ex,ey,ez) = spaceChargeCalc.calculateField(0.,r,0.)
	filed_y_arr.append((r,ey))	
	
filed_z_arr = []
for i in range(nGraphPoints):
	r = i*2*r_max_graph/(nGraphPoints) - r_max_graph
	(ex,ey,ez) = spaceChargeCalc.calculateField(0.,0.,r)
	filed_z_arr.append((r,ez))		
	

import Gnuplot

gX = Gnuplot.Gnuplot()
gX.title('Ex as function of x')
gX('set style data lines')
gX.plot(filed_x_arr)

time.sleep(1.0)

gY = Gnuplot.Gnuplot()
gY.title('Ey as function of y')
gY('set style data lines')
gY.plot(filed_y_arr)

time.sleep(1.0)

gZ = Gnuplot.Gnuplot()
gZ.title('Ez as function of z')
gZ('set style data lines')
gZ.plot(filed_z_arr)


print "Stop."

raw_input('Please press return to stop:\n')
#-------------------------------------------------	
