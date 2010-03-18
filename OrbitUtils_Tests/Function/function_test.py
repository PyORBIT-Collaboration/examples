import sys
from orbit_utils import Function

f = Function()
f.add(0.,0.5)
f.add(1.,1.0)
print "size=",f.getSize()
for i in range(f.getSize()):
	x = f.x(i)
	y = f.y(i)
	print "==== i=",i
	print "x=",x," y=",y
	print "(x,y)=",f.xy(i) 
	
f.dump()
f.normalize()
f.dump()
