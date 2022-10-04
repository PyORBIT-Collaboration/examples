##############################################################
# This script test the BaseAperture and PyBaseApertureShape
# PyBaseApertureShape class does not work by itself. It is an
# abstract class. User should implement the subclass with
# "inside(bunch,index)" method. 
##############################################################

import math
import sys

from bunch import Bunch
from aperture import BaseAperture, PyBaseApertureShape

class PyTestApertureShape(PyBaseApertureShape):
	
	def __init__(self):
		"""
		Constructor. The __init__ is necessary for PyBaseApertureShape 
		"""
		PyBaseApertureShape.__init__(self)
		self.radius = 2.7
		#---- it will be by default - here it is just to show the example
		self.centerX(0.0)
		self.centerY(0.0)
		
	def inside(self,bunch, index):
		"""
		Implementation of the inside method from PyBaseApertureShape base class
		"""
		centerX = self.centerX()
		centerY = self.centerY()
		x = bunch.x(index) - centerX
		y = bunch.y(index) - centerY
		if(math.sqrt(x**2 + y**2) < self.radius): return 1
		return 0
		
apertureShape = PyTestApertureShape()
apertureShape.name("TestPyShape")
print "Shape name=",apertureShape.name()
print "Shape type=",apertureShape.typeName()
print "Shape parameters=",apertureShape.getParamsDict()

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
	apertureShape = PyTestApertureShape()
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