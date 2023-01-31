##############################################################
# This script test the BaseAperture and ConvexApertureShape
# class which is used to combine several [x,y] points in
# one aperture that we call "convex".
##############################################################

import math
import sys

from bunch import Bunch
from aperture import BaseAperture
from aperture import PyBaseApertureShape
from aperture import PrimitiveApertureShape
from aperture import CompositeApertureShape
from aperture import ConvexApertureShape

#---- points should create a clockwise convex shape
points_arr = [[-1.,+1.],[+1.,+1.],[+1.,-1.],[-1.,-1.]]

#---- after reverse these points will give anti-clockwise shape
#---- and script will complain about this.
#points_arr.reverse()

print "points_arr = ",points_arr

apertureShape = ConvexApertureShape()
apertureShape.setPoints(points_arr)
print "Points = ",apertureShape.getPoints()

apertureShape.name("TestConvexShape")
print "Shape name=",apertureShape.name()
print "Shape type=",apertureShape.typeName()

baseAperture = BaseAperture()
baseAperture.setApertureShape(apertureShape)
baseAperture.position(11.0)
#---- by default is True (1), but you can switch it to False (0
baseAperture.onOff(True)

bunch = Bunch()
nParts = 15
for ind in xrange(nParts):
	bunch.addParticle(0.1*ind ,0.2*ind ,0.1*ind ,0.4+ind ,0.5+ind ,0.6+ind )
	
lostBunch = Bunch()	
	
baseAperture.checkBunch(bunch,lostBunch)

#---- should work also, but lostBunch will be empty
#baseAperture.checkBunch(bunch)

print "=============bunch================="
bunch.dumpBunch()

print "=============lost bunch============"
lostBunch.dumpBunch()

print "=============done=================="

sys.exit(0)

#------------------------------------------------
# Below is a memeory leak check for all classes
# Run the script and see ">top" monitor
#------------------------------------------------

count = 0
while( 1 < 2):
	
	bunch = Bunch()
	nParts = 5
	for ind in xrange(nParts):
		bunch.addParticle(0.1+ind ,0.2+ind ,0.3+ind ,0.4+ind ,0.5+ind ,0.6+ind )
		
	lostBunch = Bunch()
		
	baseAperture.checkBunch(bunch,lostBunch)
	
	apertureShape = ConvexApertureShape()
	points_arr = [[-1.,+1.],[+1.,+1.],[+1.,-1.],[-1.,-1.]]
	apertureShape.setPoints(points_arr)
	
	points_arr = [[-1.,+1.],[+1.,+1.],[+1.,-1.],[-1.,-1.]]
	apertureShape.setPoints(points_arr)
	
	apertureShape.name("TestConvexShape")
	baseAperture = BaseAperture()
	baseAperture.setApertureShape(apertureShape)
	baseAperture.position(11.0)
	baseAperture.onOff(True)
	apertureShape = baseAperture.getApertureShape()
	
	count += 1
	if(count % 1000):
		print "count=",count

print "Done."