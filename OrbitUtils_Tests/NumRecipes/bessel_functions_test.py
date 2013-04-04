#------------------------------------------------------
#This is an example of Bessel function J(n,x) and I(n,x)
#-------------------------------------------------------

import sys
import math

from orbit_utils.numrecipes import bessi0,bessi1,bessi
from orbit_utils.numrecipes import bessj0,bessj1,bessj

x = 1.0

print "bessi(0,",x,")=",bessi(0,x)
print "bessi(1,",x,")=",bessi(1,x)
print "bessi0(",x,")=",bessi0(x)
print "bessi1(",x,")=",bessi1(x)
print "============================"
print "bessj(0,",x,")=",bessj(0,x)
print "bessj(1,",x,")=",bessj(1,x)
print "bessj0(",x,")=",bessj0(x)
print "bessj1(",x,")=",bessj1(x)
print "============================"
step = 0.1
x_start = 0.
n_steps = 30
print "x  bessj0 bessj1 bessj2   bessi0 bessi1 bessi2  "
for i in range(n_steps):
	x = x_start + step*i
	res = (x,bessj(0,x),bessj(1,x),bessj(2,x),bessi(0,x),bessi(1,x),bessi(2,x))
	st = " %4.1f %6.4f %6.4f %6.4f    %6.4f %6.4f %6.4f "%res
	print st
