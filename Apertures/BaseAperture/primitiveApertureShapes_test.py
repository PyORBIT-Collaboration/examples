##############################################################
# This script test the BaseAperture and PrimitiveApertureShape
# PrimitiveApertureShape class that is used for circle, ellipse,
# and rectangular aperture shapes.
##############################################################

import math
import sys

from bunch import Bunch
from aperture import BaseAperture
from aperture import PyBaseApertureShape
from aperture import PrimitiveApertureShape

radius = 0.5
x_half_size = 0.6
y_half_size = 0.8
apertureShape = PrimitiveApertureShape("circle",radius)
#apertureShape.setParams(0.3)
print "paramsDict =",apertureShape.getParamsDict()
apertureShape = PrimitiveApertureShape("ellipse",x_half_size,y_half_size)
#apertureShape.setParams(0.2,0.4)
print "paramsDict =",apertureShape.getParamsDict()
apertureShape = PrimitiveApertureShape("rectangular",x_half_size,y_half_size)
#apertureShape.setParams(0.5,0.6)
print "paramsDict =",apertureShape.getParamsDict()
apertureShape.name("TestPyShape")
print "Shape name=",apertureShape.name()
print "Shape type=",apertureShape.typeName()

baseAperture = BaseAperture()
baseAperture.setApertureShape(apertureShape)

bunch = Bunch()
nParts = 5
for ind in xrange(nParts):
	bunch.addParticle(0.1+ind ,0.2+ind ,0.3+ind ,0.4+ind ,0.5+ind ,0.6+ind )
	
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
	
	apertureShape = PrimitiveApertureShape("circle",radius)
	print "paramsDict =",apertureShape. getParamsDict()
	apertureShape = PrimitiveApertureShape("ellipse",x_half_size,y_half_size)
	print "paramsDict =",apertureShape. getParamsDict()
	apertureShape = PrimitiveApertureShape("rectangular",x_half_size,y_half_size)
	print "paramsDict =",apertureShape. getParamsDict()
	
	apertureShape.name("TestPyShape")
	par_dict = apertureShape.getParamsDict()
	apertureShape.name()
	baseAperture = BaseAperture()
	baseAperture.setApertureShape(apertureShape)
	apertureShape = baseAperture.getApertureShape()
	
	count += 1
	if(count % 1000):
		print "count=",count

print "Done."