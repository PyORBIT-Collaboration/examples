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

radius = 0.5
x_half_size = 0.5
y_half_size = 0.5
apertureShape1 = PrimitiveApertureShape("circle",radius)
print "Shape 1 paramsDict =",apertureShape1.getParamsDict()
apertureShape2 = PrimitiveApertureShape("ellipse",x_half_size,y_half_size)
print "Shape 2 paramsDict =",apertureShape2.getParamsDict()

apertureShape2.centerX(1.0)

compositeApertureShape = CompositeApertureShape()
compositeApertureShape.addApertureShape(apertureShape1)
compositeApertureShape.addApertureShape(apertureShape2)


apertureShape_arr = compositeApertureShape.getApertureShapes()
print "===============Composite Shape==================="
for apertureShape in apertureShape_arr:
	print "debug apertureShape=",apertureShape.name()," type=",apertureShape.typeName()
print "================================================="


baseAperture = BaseAperture()
baseAperture.setApertureShape(compositeApertureShape)
baseAperture.position(11.0)
#---- by default is True (1), but you can switch it to False (0
baseAperture.onOff(True)

bunch = Bunch()
bunch.addParticle(0.1 , 0.0, 0.1, 0.0, 0.0 , 0.0)
bunch.addParticle(0.2 , 0.0, 0.2, 0.0, 0.0 , 0.0)

bunch.addParticle(1.1 , 0.0, 0.1, 0.0, 0.0 , 0.0)
bunch.addParticle(1.2 , 0.0, 0.2, 0.0, 0.0 , 0.0)

lostBunch = Bunch()	
	
baseAperture.checkBunch(bunch,lostBunch)

#---- should work also, but lostBunch will be empty
#baseAperture.checkBunch(bunch)

print "=============bunch================="
bunch.dumpBunch()

print "=============lost bunch============"
lostBunch.dumpBunch()

print "=============done=================="
