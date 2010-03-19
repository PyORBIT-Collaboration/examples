import sys
from orbit_utils import Function

f = Function()
for i in range(10):
	f.add(1.0*i,0.5*i)
	
#f.add(5.1,1000.)

f.dump()
print "const step=",f.setConstStep(1)

print "y(5.5)=",f.getY(5.5)
print "norm=",f.normalize()
f.dump()
