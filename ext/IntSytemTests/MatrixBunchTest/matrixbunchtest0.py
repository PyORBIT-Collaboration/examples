#---------------------------------------------------
# This script was prepared by Kirill Danilov during the
# summer of 2010. It will track one particle through
# one transport matrix and one non-linear element and 
# plot the phase trajectory for 2000 turns
#----------------------------------------------------


import os, sys, math, tempfile, random, numpy
from bunch import Bunch
from orbit_utils import Matrix
import intsystem
from intsystem import Nll

try:
    import Gnuplot, Gnuplot.PlotItems, Gnuplot.funcutils
except ImportError:
    # kludge in case Gnuplot hasn't been installed as a module yet:
    import __init__
    Gnuplot = __init__
    import PlotItems
    Gnuplot.PlotItems = PlotItems
    import funcutils
    Gnuplot.funcutils = funcutils

#-----------------------------------------------------
#Test of Bunch and matrix
#-----------------------------------------------------

gp = Gnuplot.Gnuplot(persist = 1)

plotx = []
ploty = []
plotxy = []

def printM(m):
    print "----matrix--- size=",m.size()
    for i in xrange(m.size()[0]):
        for j in xrange(m.size()[1]):
            print ("m(" + str(i) + "," + str(j)+")="+str(m.get(i,j)) + " "),
        print ""    


print "Start."

n = Nll()

k = .01

b = Bunch()

b.addParticle(1.,2.,3.,4.,5.,6.)

b.dumpBunch()

m = Matrix(6,6)
m.unit()
m.set(0,0,math.cos(k))
m.set(0,1,math.sin(k))
m.set(1,0,-math.sin(k))
m.set(1,1,math.cos(k))
m.set(2,2,math.cos(k))
m.set(2,3,math.sin(k))
m.set(3,2,-math.sin(k))
m.set(3,3,math.cos(k))

printM(m)

for i in range(20000):
    m.track(b); n.TRACK_EXT(b); plotx.append([b.x(0), b.y(0)]); #ploty.append([b.y(0), b.yp(0)]); ploty.append([b.x(0), b.y(0)])

gp.plot(plotx)
#gp.plot(ploty)
#gp.plot(plotxy)

print "Stop." 
