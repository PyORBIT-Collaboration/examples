##############################################################
# This script test the BaseAperture and CompositeApertureShape
# class which is used to combine several BaseApertureShapes in
# one aperture that we call "composite".
##############################################################

import math
import sys

from bunch import Bunch
from aperture import BaseAperture
from aperture import PyBaseApertureShape
from aperture import PrimitiveApertureShape
from aperture import CompositeApertureShape

radius = 0.6
x_half_size = 0.6
y_half_size = 0.6
apertureShape1 = PrimitiveApertureShape("circle",radius)
print "Shape 1 paramsDict =",apertureShape1.getParamsDict()
apertureShape2 = PrimitiveApertureShape("ellipse",x_half_size,0.5*y_half_size)
print "Shape 2 paramsDict =",apertureShape2.getParamsDict()
apertureShape3 = PrimitiveApertureShape("rectangular",x_half_size,0.5*y_half_size)
print "Shape 3 paramsDict =",apertureShape2.getParamsDict()

compositeApertureShape = CompositeApertureShape()
compositeApertureShape.addApertureShape(apertureShape1)
compositeApertureShape.addApertureShape(apertureShape2)
compositeApertureShape.addApertureShape(apertureShape3)


apertureShape_arr = compositeApertureShape.getApertureShapes()
print "===============Composite Shape==================="
for apertureShape in apertureShape_arr:
	print "debug apertureShape=",apertureShape.name()," type=",apertureShape.typeName()
print "================================================="

compositeApertureShape.name("TestPyCompositeShape")
print "Shape name=",compositeApertureShape.name()
print "Shape type=",compositeApertureShape.typeName()

baseAperture = BaseAperture()
baseAperture.setApertureShape(compositeApertureShape)

bunch = Bunch()
nParts = 7
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

#sys.exit(0)

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
	
	apertureShape1 = PrimitiveApertureShape("circle",radius)
	apertureShape2 = PrimitiveApertureShape("ellipse",x_half_size,0.5*y_half_size)
	apertureShape3 = PrimitiveApertureShape("rectangular",x_half_size,0.5*y_half_size)
	
	compositeApertureShape = CompositeApertureShape()
	compositeApertureShape.addApertureShape(apertureShape1)
	compositeApertureShape.addApertureShape(apertureShape2)
	compositeApertureShape.addApertureShape(apertureShape3)
	
	compositeApertureShape.name("TestCompositeShape")
	baseAperture = BaseAperture()
	baseAperture.setApertureShape(compositeApertureShape)
	apertureShape = baseAperture.getApertureShape()
	
	count += 1
	if(count % 1000):
		print "count=",count

print "Done."